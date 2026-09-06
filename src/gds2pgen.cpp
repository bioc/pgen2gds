// ===========================================================
//
// gds2pgen.cpp: Format Conversion from GDS to PLINK2 PGEN
//
// Copyright (C) 2026    Xiuwen Zheng
//
// pgen2gds is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License Version 3 as
// published by the Free Software Foundation.
//
// pgen2gds is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along
// with pgen2gds.
// If not, see <http://www.gnu.org/licenses/>.
//
// The PGEN layout implemented here follows the PLINK2 specification
// (plink-ng/2.0/include/pgenlib_misc.h, "The actual format"):
//   * storage mode 0x10, one index block per 65536 variants;
//   * per-variant record types 0 (dense 2-bit), 1 (1-bit + difflist) and
//     4/6/7 (difflist), chosen with the same thresholds as plink2;
//   * multiallelic hardcall patches (vrtype bit 0x08) and hardcall phase
//     (vrtype bit 0x10);
//   * no LD compression (vrtype 2/3) and no dosage tracks.

#include <cstdio>
#include <cstring>
#include <stdint.h>
#include <string>
#include <vector>
#include <stdexcept>
#include <algorithm>

#define R_NO_REMAP
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#ifndef COREARRAY_DLL_EXPORT
#   define COREARRAY_DLL_EXPORT
#endif

using namespace std;


namespace {

static const uint32_t kVblockSize = 65536u;  ///< # of variants per index block
static const uint32_t kGroupSize  = 64u;     ///< difflist group size
static const uint8_t  kMissingRaw = 0xFF;    ///< SeqArray raw missing genotype

/// # of bytes needed to store the nonzero value x
static inline uint32_t bytes_to_represent(uint32_t x)
{
	uint32_t n = 1;
	while (x >= 256u) { x >>= 8; n++; }
	return n;
}

static inline uint32_t div_up(uint32_t a, uint32_t b) { return (a + b - 1) / b; }

typedef vector<uint8_t> ByteBuf;

static inline void put_le(ByteBuf &b, uint64_t v, uint32_t nbytes)
{
	for (uint32_t i=0; i < nbytes; i++)
	{
		b.push_back((uint8_t)(v & 0xFF));
		v >>= 8;
	}
}

/// plink2 Vint32Append: 7 bits per byte, low bits first
static inline void put_vint(ByteBuf &b, uint32_t v)
{
	while (v > 127u)
	{
		b.push_back((uint8_t)((v & 127u) | 128u));
		v >>= 7;
	}
	b.push_back((uint8_t)v);
}

/// Append bits to a byte buffer, least significant bit first
class BitWriter
{
public:
	BitWriter(ByteBuf &buf): fBuf(buf), fCur(0), fNBit(0) {}
	inline void put(uint32_t bit)
	{
		if (bit) fCur |= (1u << fNBit);
		if (++fNBit == 8) flush_byte();
	}
	inline void put_bits(uint32_t v, uint32_t nbit)
	{
		for (uint32_t i=0; i < nbit; i++) put((v >> i) & 1u);
	}
	inline void finish() { if (fNBit) flush_byte(); }
private:
	ByteBuf &fBuf;
	uint32_t fCur, fNBit;
	inline void flush_byte() { fBuf.push_back((uint8_t)fCur); fCur = 0; fNBit = 0; }
};

/// Append a difflist (raregeno != NULL) or a deltalist (raregeno == NULL):
///   <vint length> <group start sample ids> <extra byte counts>
///   [<2-bit genotypes>] <vint sample-id deltas, per group>
static void put_difflist(ByteBuf &rec, const uint32_t *ids,
	const uint8_t *raregeno, uint32_t n, uint32_t sid_bytes)
{
	put_vint(rec, n);
	if (n == 0) return;
	const uint32_t group_ct = div_up(n, kGroupSize);
	for (uint32_t g=0; g < group_ct; g++)
		put_le(rec, ids[g*kGroupSize], sid_bytes);
	const size_t extra_pos = rec.size();
	rec.resize(extra_pos + group_ct - 1, 0);
	if (raregeno)
	{
		BitWriter bw(rec);
		for (uint32_t i=0; i < n; i++) bw.put_bits(raregeno[i] & 3u, 2);
		bw.finish();
	}
	for (uint32_t g=0; g < group_ct; g++)
	{
		const uint32_t st = g * kGroupSize;
		const uint32_t ed = min(st + kGroupSize, n);
		const size_t seg_start = rec.size();
		for (uint32_t i=st+1; i < ed; i++)
			put_vint(rec, ids[i] - ids[i-1]);
		if (g + 1 < group_ct)
			rec[extra_pos + g] = (uint8_t)(rec.size() - seg_start - (kGroupSize - 1));
	}
}


// ===================================================================== //

/// PGEN writer (single pass, header and index filled in at the end)
class CPgenWriter
{
public:
	enum { cDense=0, cOnebit=1, cDifflist=2, cMultiallelic=3, cPhased=4, cNumCnt=5 };

	CPgenWriter(const char *fn, uint32_t nsamp, uint32_t nvar, bool use_phase,
		uint32_t max_allele_ct);
	~CPgenWriter();

	/// append nv variants; geno: 2*nsamp*nv raw allele indices (0xFF missing),
	/// phase: nsamp*nv raw flags or NULL, allele_ct: nv integers
	void AppendBlock(const uint8_t *geno, const uint8_t *phase,
		const int *allele_ct, uint32_t nv, int half_call);
	/// write header and index, close the file
	void Finish();
	/// close and remove the incomplete file
	void Abort();

	const uint32_t *Counts() const { return fCnt; }
	uint32_t NumVariant() const { return fVidx; }
	uint32_t LenBytes() const { return fLenBytes; }

private:
	FILE *fFile;
	string fFileName;
	uint32_t fN, fM, fVidx;
	bool fUsePhase;
	uint32_t fMaxAlleleCt, fSidBytes, fLenBytes;
	uint64_t fHeaderSize, fFpos, fMaxVrecLen;
	vector<uint8_t> fVrtype;        ///< one byte per variant
	vector<uint8_t> fVrecLen;       ///< fLenBytes per variant
	vector<uint64_t> fVblockFpos;   ///< file offset of the first record in each block
	uint32_t fCnt[cNumCnt];

	// per-variant scratch
	vector<uint8_t> fG;             ///< main 2-bit codes, one byte per sample
	vector<uint8_t> fHet, fPhased, fSwap;
	vector<uint32_t> fP01, fP10;    ///< sample indices of patched genotypes
	vector<uint8_t> fP01v, fP10a, fP10b;  ///< allele codes of patched genotypes
	vector<uint32_t> fIds;
	vector<uint8_t> fRG;
	ByteBuf fRec, fTmp;

	uint64_t max_vrec_len() const;
	void write_zeros(uint64_t n);
	uint8_t encode_main();
	void encode_aux_multiallelic(uint32_t allele_ct);
	bool encode_aux_phase();
};


CPgenWriter::CPgenWriter(const char *fn, uint32_t nsamp, uint32_t nvar,
	bool use_phase, uint32_t max_allele_ct)
{
	fFile = NULL;
	fFileName = fn;
	fN = nsamp; fM = nvar; fVidx = 0;
	fUsePhase = use_phase;
	fMaxAlleleCt = max_allele_ct;
	memset(fCnt, 0, sizeof(fCnt));
	if (fN == 0) throw runtime_error("No sample.");
	if (fM == 0) throw runtime_error("No variant.");
	if (fMaxAlleleCt < 2) fMaxAlleleCt = 2;
	if (fMaxAlleleCt > 255)
		throw runtime_error("PLINK2 supports at most 255 alleles per variant.");

	fSidBytes = bytes_to_represent(fN);
	fMaxVrecLen = max_vrec_len();
	if (fMaxVrecLen >= 0xFFFFFFFFull)
		throw runtime_error("Too many samples for a PGEN variant record.");
	fLenBytes = bytes_to_represent((uint32_t)fMaxVrecLen);

	const uint32_t vblock_ct = div_up(fM, kVblockSize);
	fHeaderSize = 12u + 8ull*vblock_ct +
		(fUsePhase ? (uint64_t)fM : (uint64_t)div_up(fM, 2)) +
		(uint64_t)fM * fLenBytes;

	fVrtype.reserve(fM);
	fVrecLen.reserve((size_t)fM * fLenBytes);
	fVblockFpos.reserve(vblock_ct);
	fG.resize(fN); fHet.resize(fN); fPhased.resize(fN); fSwap.resize(fN);

	fFile = fopen(fn, "wb");
	if (!fFile) throw runtime_error(string("Fail to create '") + fn + "'.");
	write_zeros(fHeaderSize);
	fFpos = fHeaderSize;
}

CPgenWriter::~CPgenWriter()
{
	if (fFile) { fclose(fFile); fFile = NULL; }
}

/// an upper bound of the record length for the encodings used here
uint64_t CPgenWriter::max_vrec_len() const
{
	const uint64_t N = fN;
	uint64_t n = (N + 3) / 4;  // dense main track (alternatives never exceed it)
	if (fMaxAlleleCt > 2)
	{
		// patch sets: bitarray (N/8) or deltalist (< N/9 entries)
		const uint64_t nlist = N/9 + 1;
		const uint64_t set_bound = 5 + (nlist/kGroupSize + 1)*(fSidBytes + 1) +
			nlist + N/32 + (N + 7)/8;
		// patch values
		uint64_t v01, v10;
		if (fMaxAlleleCt <= 3) v01 = 0;
		else if (fMaxAlleleCt == 4) v01 = (N + 7)/8;
		else if (fMaxAlleleCt <= 6) v01 = (N + 3)/4;
		else if (fMaxAlleleCt <= 18) v01 = (N + 1)/2;
		else v01 = N;
		if (fMaxAlleleCt == 3) v10 = (N + 7)/8;
		else if (fMaxAlleleCt <= 5) v10 = (N + 1)/2;
		else if (fMaxAlleleCt <= 17) v10 = N;
		else v10 = 2*N;
		n += 1 + 2*set_bound + v01 + v10;
	}
	if (fUsePhase)
		n += 2 + 2*((N + 7)/8);
	return n;
}

void CPgenWriter::write_zeros(uint64_t n)
{
	static const uint8_t zeros[65536] = { 0 };
	while (n > 0)
	{
		size_t m = (n > sizeof(zeros)) ? sizeof(zeros) : (size_t)n;
		if (fwrite(zeros, 1, m, fFile) != m)
			throw runtime_error("Fail to write the PGEN file.");
		n -= m;
	}
}

void CPgenWriter::Abort()
{
	if (fFile) { fclose(fFile); fFile = NULL; }
	::remove(fFileName.c_str());
}


/// encode the main data track into fRec, return the low 3 bits of vrtype
uint8_t CPgenWriter::encode_main()
{
	const uint32_t N = fN;
	const uint8_t *g = &fG[0];
	uint32_t cnt[4] = { 0, 0, 0, 0 };
	for (uint32_t i=0; i < N; i++) cnt[g[i]]++;

	// the two most common genotypes (plink2 PwcAppendBiallelicGenovecMain)
	uint32_t most = (cnt[1] > cnt[0]) ? 1 : 0;
	uint32_t second = 1 - most;
	uint32_t largest = cnt[most], second_largest = cnt[second];
	for (uint32_t c=2; c < 4; c++)
	{
		if (cnt[c] > second_largest)
		{
			if (cnt[c] > largest)
			{
				second_largest = largest; second = most;
				largest = cnt[c]; most = c;
			} else {
				second_largest = cnt[c]; second = c;
			}
		}
	}
	const uint32_t difflist_len = N - largest;
	const uint32_t rare2 = difflist_len - second_largest;
	const uint32_t n_d8 = N / 8, n_d64 = N / 64;
	uint32_t max_difflist_len = n_d8 - 2*n_d64 + rare2;
	if (max_difflist_len > n_d8) max_difflist_len = n_d8;
	const bool difflist_viable = (most != 1) && (difflist_len <= max_difflist_len);
	const size_t dense_size = (N + 3) / 4;

	fRec.clear();
	if (!difflist_viable && (rare2 < N/16))
	{
		// 1-bit + difflist representation
		const uint32_t larger = max(most, second), smaller = min(most, second);
		fRec.push_back((uint8_t)(larger + smaller*3));
		{
			BitWriter bw(fRec);
			for (uint32_t i=0; i < N; i++) bw.put(g[i] == larger);
			bw.finish();
		}
		fIds.clear(); fRG.clear();
		for (uint32_t i=0; i < N; i++)
		{
			if (g[i] != larger && g[i] != smaller)
				{ fIds.push_back(i); fRG.push_back(g[i]); }
		}
		put_difflist(fRec, fIds.empty() ? NULL : &fIds[0],
			fRG.empty() ? NULL : &fRG[0], (uint32_t)fIds.size(), fSidBytes);
		if (fRec.size() <= dense_size)
			{ fCnt[cOnebit]++; return 1; }
	} else if (difflist_viable)
	{
		fIds.clear(); fRG.clear();
		for (uint32_t i=0; i < N; i++)
		{
			if (g[i] != most)
				{ fIds.push_back(i); fRG.push_back(g[i]); }
		}
		put_difflist(fRec, fIds.empty() ? NULL : &fIds[0],
			fRG.empty() ? NULL : &fRG[0], (uint32_t)fIds.size(), fSidBytes);
		if (fRec.size() <= dense_size)
			{ fCnt[cDifflist]++; return (uint8_t)(4 + most); }
	}

	// dense 2-bit encoding
	fRec.clear();
	BitWriter bw(fRec);
	for (uint32_t i=0; i < N; i++) bw.put_bits(g[i], 2);
	bw.finish();
	fCnt[cDense]++;
	return 0;
}


/// append auxiliary data track #1 (multiallelic hardcall patches)
void CPgenWriter::encode_aux_multiallelic(uint32_t allele_ct)
{
	const uint32_t N = fN;
	const uint8_t *g = &fG[0];
	const size_t fmt_pos = fRec.size();
	fRec.push_back(0);
	uint8_t format_byte = 0;

	// ref/altx patches
	const uint32_t p01_ct = (uint32_t)fP01.size();
	if (p01_ct == 0)
	{
		format_byte = 15;
	} else {
		uint32_t het_ct = 0;
		for (uint32_t i=0; i < N; i++) het_ct += (g[i] == 1);
		if (p01_ct < het_ct / 9)
		{
			put_difflist(fRec, &fP01[0], NULL, p01_ct, fSidBytes);
			format_byte = 1;
		} else {
			BitWriter bw(fRec);
			uint32_t k = 0;
			for (uint32_t i=0; i < N; i++)
			{
				if (g[i] == 1)
				{
					const bool hit = (k < p01_ct) && (fP01[k] == i);
					if (hit) k++;
					bw.put(hit);
				}
			}
			bw.finish();
		}
		if (allele_ct > 3)
		{
			if (allele_ct <= 18)
			{
				const uint32_t width = (allele_ct == 4) ? 1 : ((allele_ct <= 6) ? 2 : 4);
				BitWriter bw(fRec);
				for (uint32_t k=0; k < p01_ct; k++) bw.put_bits(fP01v[k] - 2, width);
				bw.finish();
			} else {
				for (uint32_t k=0; k < p01_ct; k++) fRec.push_back(fP01v[k] - 2);
			}
		}
	}

	// altx/alty patches
	const uint32_t p10_ct = (uint32_t)fP10.size();
	if (p10_ct == 0)
	{
		format_byte |= 0xF0;
	} else {
		uint32_t altxy_ct = 0;
		for (uint32_t i=0; i < N; i++) altxy_ct += (g[i] == 2);
		if (p10_ct < altxy_ct / 9)
		{
			put_difflist(fRec, &fP10[0], NULL, p10_ct, fSidBytes);
			format_byte |= 0x10;
		} else {
			BitWriter bw(fRec);
			uint32_t k = 0;
			for (uint32_t i=0; i < N; i++)
			{
				if (g[i] == 2)
				{
					const bool hit = (k < p10_ct) && (fP10[k] == i);
					if (hit) k++;
					bw.put(hit);
				}
			}
			bw.finish();
		}
		if (allele_ct == 3)
		{
			BitWriter bw(fRec);
			for (uint32_t k=0; k < p10_ct; k++) bw.put(fP10a[k] - 1);
			bw.finish();
		} else if (allele_ct <= 5)
		{
			BitWriter bw(fRec);
			for (uint32_t k=0; k < p10_ct; k++)
			{
				bw.put_bits(fP10a[k] - 1, 2);
				bw.put_bits(fP10b[k] - 1, 2);
			}
			bw.finish();
		} else if (allele_ct <= 17)
		{
			for (uint32_t k=0; k < p10_ct; k++)
				fRec.push_back((uint8_t)((fP10a[k] - 1) | ((fP10b[k] - 1) << 4)));
		} else {
			for (uint32_t k=0; k < p10_ct; k++)
			{
				fRec.push_back(fP10a[k] - 1);
				fRec.push_back(fP10b[k] - 1);
			}
		}
	}

	fRec[fmt_pos] = format_byte;
	fCnt[cMultiallelic]++;
}


/// append auxiliary data track #2 (hardcall phase), return false if no phased het
bool CPgenWriter::encode_aux_phase()
{
	const uint32_t N = fN;
	uint32_t het_ct = 0, phased_ct = 0;
	for (uint32_t i=0; i < N; i++)
	{
		if (fHet[i]) { het_ct++; if (fPhased[i]) phased_ct++; }
	}
	if (phased_ct == 0) return false;
	BitWriter bw(fRec);
	if (phased_ct == het_ct)
	{
		// no explicit phasepresent bitarray
		bw.put(0);
		for (uint32_t i=0; i < N; i++)
			if (fHet[i]) bw.put(fSwap[i]);
		bw.finish();
	} else {
		bw.put(1);
		for (uint32_t i=0; i < N; i++)
			if (fHet[i]) bw.put(fPhased[i]);
		bw.finish();
		BitWriter bw2(fRec);
		for (uint32_t i=0; i < N; i++)
			if (fHet[i] && fPhased[i]) bw2.put(fSwap[i]);
		bw2.finish();
	}
	fCnt[cPhased]++;
	return true;
}


void CPgenWriter::AppendBlock(const uint8_t *geno, const uint8_t *phase,
	const int *allele_ct, uint32_t nv, int half_call)
{
	if (!fFile) throw runtime_error("The PGEN writer has been closed.");
	const uint32_t N = fN;
	for (uint32_t v=0; v < nv; v++)
	{
		if (fVidx >= fM)
			throw runtime_error("More variants than expected.");
		const int nallele = allele_ct[v];
		if (nallele < 1 || nallele > (int)fMaxAlleleCt)
			throw runtime_error("Invalid number of alleles.");
		const uint8_t *pg = geno + (size_t)2 * N * v;
		const uint8_t *pp = phase ? (phase + (size_t)N * v) : NULL;

		// per-sample genotype codes and patches
		fP01.clear(); fP10.clear(); fP01v.clear(); fP10a.clear(); fP10b.clear();
		for (uint32_t i=0; i < N; i++, pg+=2)
		{
			uint32_t a1 = pg[0], a2 = pg[1];
			uint8_t code, het = 0, swap = 0;
			if (a1 == kMissingRaw || a2 == kMissingRaw)
			{
				if ((a1 != a2) && (half_call == 1))
				{
					// haploid call: treat as homozygous
					if (a1 == kMissingRaw) a1 = a2; else a2 = a1;
				} else {
					a1 = a2 = kMissingRaw;
				}
			}
			if (a1 == kMissingRaw)
			{
				code = 3;
			} else {
				if ((int)a1 >= nallele || (int)a2 >= nallele)
					throw runtime_error("Invalid allele index in the genotypes.");
				if (a1 == 0 && a2 == 0)
				{
					code = 0;
				} else if (a1 == 0 || a2 == 0)
				{
					code = 1; het = 1;
					const uint32_t k = a1 ? a1 : a2;
					swap = (a1 != 0);
					if (k >= 2)
						{ fP01.push_back(i); fP01v.push_back((uint8_t)k); }
				} else {
					code = 2;
					const uint32_t j = min(a1, a2), k = max(a1, a2);
					if (j != k) { het = 1; swap = (a1 > a2); }
					if (!(j == 1 && k == 1))
					{
						fP10.push_back(i);
						fP10a.push_back((uint8_t)j); fP10b.push_back((uint8_t)k);
					}
				}
			}
			fG[i] = code; fHet[i] = het; fSwap[i] = swap;
			fPhased[i] = (pp && het) ? (pp[i] != 0) : 0;
		}

		// encode
		uint8_t vrtype = encode_main();
		if (!fP01.empty() || !fP10.empty())
		{
			encode_aux_multiallelic((uint32_t)nallele);
			vrtype |= 0x08;
		}
		if (fUsePhase && encode_aux_phase())
			vrtype |= 0x10;
		if (fRec.size() > fMaxVrecLen)
			throw runtime_error("Internal error: variant record is too long.");

		// write the record
		if ((fVidx % kVblockSize) == 0) fVblockFpos.push_back(fFpos);
		if (!fRec.empty())
		{
			if (fwrite(&fRec[0], 1, fRec.size(), fFile) != fRec.size())
				throw runtime_error("Fail to write the PGEN file.");
		}
		fFpos += fRec.size();
		fVrtype.push_back(vrtype);
		put_le(fVrecLen, fRec.size(), fLenBytes);
		fVidx++;
	}
}


void CPgenWriter::Finish()
{
	if (!fFile) throw runtime_error("The PGEN writer has been closed.");
	if (fVidx != fM)
		throw runtime_error("Fewer variants than expected.");
	if (fseek(fFile, 0, SEEK_SET) != 0)
		throw runtime_error("Fail to seek in the PGEN file.");

	ByteBuf hd;
	hd.reserve(fHeaderSize);
	hd.push_back(0x6C); hd.push_back(0x1B); hd.push_back(0x10);
	put_le(hd, fM, 4);
	put_le(hd, fN, 4);
	// control byte: record length width, 8-bit vrtypes if phase present,
	// nonref flags storage 1 (all REF alleles are trusted)
	hd.push_back((uint8_t)((fLenBytes - 1) | (fUsePhase ? 4 : 0) | (1 << 6)));
	for (size_t i=0; i < fVblockFpos.size(); i++) put_le(hd, fVblockFpos[i], 8);
	for (uint32_t st=0; st < fM; st+=kVblockSize)
	{
		const uint32_t bs = min(kVblockSize, fM - st);
		if (fUsePhase)
		{
			hd.insert(hd.end(), fVrtype.begin() + st, fVrtype.begin() + st + bs);
		} else {
			for (uint32_t i=0; i < bs; i+=2)
			{
				uint8_t b = fVrtype[st + i] & 0x0F;
				if (i + 1 < bs) b |= (fVrtype[st + i + 1] & 0x0F) << 4;
				hd.push_back(b);
			}
		}
		hd.insert(hd.end(), fVrecLen.begin() + (size_t)st*fLenBytes,
			fVrecLen.begin() + (size_t)(st + bs)*fLenBytes);
	}
	if (hd.size() != fHeaderSize)
		throw runtime_error("Internal error: header size mismatch.");
	if (fwrite(&hd[0], 1, hd.size(), fFile) != hd.size())
		throw runtime_error("Fail to write the PGEN file.");
	if (fclose(fFile) != 0)
		{ fFile = NULL; throw runtime_error("Fail to close the PGEN file."); }
	fFile = NULL;
}


// ===================================================================== //

static char err_msg[1024];

static CPgenWriter *get_writer(SEXP ptr, bool allow_null=false)
{
	if (TYPEOF(ptr) != EXTPTRSXP)
		Rf_error("Invalid PGEN writer object.");
	CPgenWriter *p = (CPgenWriter*)R_ExternalPtrAddr(ptr);
	if (!p && !allow_null)
		Rf_error("The PGEN writer has been closed.");
	return p;
}

static void pgen_writer_finalizer(SEXP ptr)
{
	CPgenWriter *p = (CPgenWriter*)R_ExternalPtrAddr(ptr);
	if (p)
	{
		delete p;
		R_ClearExternalPtr(ptr);
	}
}

} // namespace


extern "C"
{

/// create a PGEN writer: SEQ_PGEN_Writer_Open(fn, nsamp, nvar, use.phase, max.allele)
COREARRAY_DLL_EXPORT SEXP SEQ_PGEN_Writer_Open(SEXP fn, SEXP nsamp, SEXP nvar,
	SEXP use_phase, SEXP max_allele)
{
	const char *s = CHAR(STRING_ELT(fn, 0));
	CPgenWriter *p = NULL;
	err_msg[0] = 0;
	try {
		p = new CPgenWriter(s, (uint32_t)Rf_asInteger(nsamp),
			(uint32_t)Rf_asInteger(nvar), Rf_asLogical(use_phase)==TRUE,
			(uint32_t)Rf_asInteger(max_allele));
	} catch (exception &e) {
		strncpy(err_msg, e.what(), sizeof(err_msg)-1);
	}
	if (err_msg[0]) Rf_error("%s", err_msg);
	SEXP ptr = PROTECT(R_MakeExternalPtr(p, R_NilValue, R_NilValue));
	R_RegisterCFinalizerEx(ptr, pgen_writer_finalizer, TRUE);
	UNPROTECT(1);
	return ptr;
}

/// append a block of variants:
/// SEQ_PGEN_Writer_Block(ptr, geno[2,N,b] raw, phase[N,b] raw or NULL,
///     allele_ct[b] int, half_call int)
COREARRAY_DLL_EXPORT SEXP SEQ_PGEN_Writer_Block(SEXP ptr, SEXP geno, SEXP phase,
	SEXP allele_ct, SEXP half_call)
{
	CPgenWriter *p = get_writer(ptr);
	if (TYPEOF(geno) != RAWSXP)
		Rf_error("'genotype' should be a raw array.");
	if (TYPEOF(allele_ct) != INTSXP)
		Rf_error("The number of alleles should be an integer vector.");
	const uint32_t nv = (uint32_t)Rf_length(allele_ct);
	const uint64_t N = (uint64_t)Rf_length(geno) / 2;
	if ((uint64_t)Rf_length(geno) != 2ull * N || (nv > 0 && N % nv != 0))
		Rf_error("Invalid dimension of the genotype array.");
	const uint8_t *pp = NULL;
	if (!Rf_isNull(phase))
	{
		if (TYPEOF(phase) != RAWSXP || (uint64_t)Rf_length(phase) != N)
			Rf_error("Invalid dimension of the phase array.");
		pp = RAW(phase);
	}
	err_msg[0] = 0;
	try {
		p->AppendBlock(RAW(geno), pp, INTEGER(allele_ct), nv,
			Rf_asInteger(half_call));
	} catch (exception &e) {
		strncpy(err_msg, e.what(), sizeof(err_msg)-1);
	}
	if (err_msg[0]) Rf_error("%s", err_msg);
	return R_NilValue;
}

/// finish (finalize=TRUE) or abort (finalize=FALSE) the writer;
/// return c(# of variants, dense, onebit, difflist, multiallelic, phased,
///     vrec_len bytes) when finalizing
COREARRAY_DLL_EXPORT SEXP SEQ_PGEN_Writer_Close(SEXP ptr, SEXP finalize)
{
	CPgenWriter *p = get_writer(ptr, true);
	if (!p) return R_NilValue;
	SEXP rv = R_NilValue;
	err_msg[0] = 0;
	if (Rf_asLogical(finalize) == TRUE)
	{
		try {
			p->Finish();
		} catch (exception &e) {
			strncpy(err_msg, e.what(), sizeof(err_msg)-1);
		}
		if (!err_msg[0])
		{
			rv = PROTECT(Rf_allocVector(INTSXP, 7));
			INTEGER(rv)[0] = p->NumVariant();
			for (int i=0; i < CPgenWriter::cNumCnt; i++)
				INTEGER(rv)[i+1] = p->Counts()[i];
			INTEGER(rv)[6] = p->LenBytes();
			UNPROTECT(1);
		}
	} else {
		p->Abort();
	}
	delete p;
	R_ClearExternalPtr(ptr);
	if (err_msg[0]) Rf_error("%s", err_msg);
	return rv;
}

} // extern "C"
