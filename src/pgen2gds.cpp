// ===========================================================
//
// pgen2gds.cpp: Format Conversion from PLINK2 PGEN to GDS
//
// Copyright (C) 2025    Xiuwen Zheng
//
// pgen2gds is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License Version 3 as
// published by the Free Software Foundation.
//
// pgen2gds is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU General Public License along
// with pgen2gds.
// If not, see <http://www.gnu.org/licenses/>.


#include <R_GDS_CPP.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#define class xclass
#define private xprivate
#include <R_ext/Connections.h>
#undef class
#undef private

#include <vector>
#include <ctime>

using namespace std;
using namespace CoreArray;


// ===================================================================== //

/// Progress object
class COREARRAY_DLL_LOCAL CProgress
{
public:
	CProgress();
	CProgress(C_Int64 count, SEXP conn);

	void Reset(C_Int64 count);
	void Forward(C_Int64 val);
	void ShowProgress();

protected:
	C_Int64 fTotalCount;  ///< the total number
	C_Int64 fCounter;  ///< the current counter
	double _start, _step;
	C_Int64 _hit;
	vector< pair<double, time_t> > _timer;
	time_t _start_time, _last_time, _last_check_time;
	Rconnection progress_conn;
};

static const int PROGRESS_BAR_CHAR_NUM = 50;
static const double S_MIN  =  60;
static const double S_HOUR =  60 * S_MIN;
static const double S_DAY  =  24 * S_HOUR;
static const double S_YEAR = 365 * S_DAY;

static int job_index = 0;
static int job_count = 0;

static const char *time_str(double s)
{
	if (R_FINITE(s))
	{
		static char buffer[64];
		if (s < S_MIN)
			snprintf(buffer, sizeof(buffer), "%.0fs", s);
		else if (s < S_HOUR)
			snprintf(buffer, sizeof(buffer), "%.1fm", s/S_MIN);
		else if (s < S_DAY)
			snprintf(buffer, sizeof(buffer), "%.1fh", s/S_HOUR);
		else if (s < S_YEAR)
			snprintf(buffer, sizeof(buffer), "%.1fd", s/S_DAY);
		else
			snprintf(buffer, sizeof(buffer), "%.1f years", s/S_YEAR);
		return buffer;
	} else
		return "---";
}

extern "C"
{
	static void chkIntFn(void *dummy) { R_CheckUserInterrupt(); }
}

static bool CheckInterrupt()
{
	// this will call the above in a top-level context so it won't longjmp-out of your context
	return (R_ToplevelExec(chkIntFn, NULL) == FALSE);
}

inline static void put_text(Rconnection conn, const char *fmt, ...)
{
	va_list args;
	va_start(args, fmt);
	(*conn->vfprintf)(conn, fmt, args);
	va_end(args);
}



CProgress::CProgress()
{
	fTotalCount = 0;
	fCounter = 0;
	progress_conn = NULL;
}

CProgress::CProgress(C_Int64 count, SEXP conn)
{
	fTotalCount = 0;
	fCounter = 0;
	progress_conn = (!conn || Rf_isNull(conn)) ? NULL : R_GetConnection(conn);
	Reset(count);
}

void CProgress::Reset(C_Int64 count)
{
	bool flag = (fTotalCount==0) || (fCounter > 0);
	fTotalCount = count;
	fCounter = 0;
	if (count > 0)
	{
		int n = 100;
		if (n > count) n = count;
		if (n < 1) n = 1;
		_start = _step = (double)count / n;
		_hit = (C_Int64)(_start);
		double percent = (double)fCounter / count;
		time_t s; time(&s);
		_start_time = _last_time = _last_check_time = s;
		_timer.clear();
		_timer.reserve(128);
		_timer.push_back(pair<double, time_t>(percent, s));
		if (flag) ShowProgress();
	}
}

void CProgress::Forward(C_Int64 val)
{
	if (fTotalCount > 0)
	{
		fCounter += val;
		if (fCounter >= _hit)
		{
			do {
				_start += _step;
				_hit = (C_Int64)(_start);
			} while (fCounter >= _hit);
			if (_hit > fTotalCount) _hit = fTotalCount;
			ShowProgress();
		}
		// check whether user interrupts
		time_t now; time(&now);
		if (difftime(now, _last_check_time) > 0.25)
		{
			_last_check_time = now;
			if (CheckInterrupt())
				throw ErrCoreArray("User interrupts!");
		}
	}
}

void CProgress::ShowProgress()
{
	if (fTotalCount > 0)
	{
		char bar[PROGRESS_BAR_CHAR_NUM + 1];
		double p = (double)fCounter / fTotalCount;
		int n = (int)round(p * PROGRESS_BAR_CHAR_NUM);
		memset(bar, '.', sizeof(bar));
		memset(bar, '=', n);
		if ((fCounter > 0) && (n < PROGRESS_BAR_CHAR_NUM))
			bar[n] = '>';
		bar[PROGRESS_BAR_CHAR_NUM] = 0;

		// ETC: estimated time to complete
		n = (int)_timer.size() - 20;  // 20% as a sliding window size
		if (n < 0) n = 0;

		time_t now; time(&now);
		_timer.push_back(pair<double, time_t>(p, now));

		// in seconds
		double interval = difftime(now, _last_time);
		double s = difftime(now, _timer[n].second);
		double diff = p - _timer[n].first;
		if (diff > 0)
			s = s / diff * (1 - p);
		else
			s = R_NaN;
		p *= 100;

		// show
		_last_time = now;
		if (fCounter >= fTotalCount)
		{
			s = difftime(_last_time, _start_time);
			if (!progress_conn)
			{
				if (job_count > 0)
				{
					Rprintf("\r[%s] 100%% (%d/%d), completed, %s\n", bar,
						job_index, job_count, time_str(s));
				} else {
					Rprintf("\r[%s] 100%%, completed, %s\n", bar, time_str(s));
				}
			} else {
				if (job_count > 0)
				{
					put_text(progress_conn, "[%s] 100%% (%d/%d), completed, %s\n",
						bar, job_index, job_count, time_str(s));
				} else {
					put_text(progress_conn, "[%s] 100%%, completed, %s\n", bar,
						time_str(s));
				}
				(*progress_conn->fflush)(progress_conn);
			}
		} else if ((interval >= 5) || (fCounter <= 0))
		{
			if (!progress_conn)
			{
				if (job_count > 0)
				{
					Rprintf("\r[%s] %2.0f%% (%d/%d), ETC: %s        ", bar,
						p, job_index, job_count, time_str(s));
				} else {
					Rprintf("\r[%s] %2.0f%%, ETC: %s        ", bar, p,
						time_str(s));
				}
			} else {
				if (job_count > 0)
				{
					put_text(progress_conn, "[%s] %2.0f%% (%d/%d), ETC: %s\n",
						bar, p, job_index, job_count, time_str(s));
				} else {
					put_text(progress_conn, "[%s] %2.0f%%, ETC: %s\n", bar, p,
						time_str(s));
				}
				(*progress_conn->fflush)(progress_conn);
			}
			// fflush(stdout);
		}
	}
}



extern "C"
{

/// Import a pgen file
COREARRAY_DLL_EXPORT SEXP SEQ_SetJobStatus(SEXP index, SEXP count)
{
	job_index = Rf_asInteger(index);
	job_count = Rf_asInteger(count);
	return R_NilValue;
}


/// Import a pgen file
COREARRAY_DLL_EXPORT SEXP SEQ_PGEN_Geno_Import(
	SEXP R_read_gt_fc, SEXP R_allele_num_fc, SEXP R_buf, SEXP R_ii, SEXP R_ia,
	SEXP R_variant_sel, SEXP R_env, SEXP gds_root, SEXP Start, SEXP Count,
	SEXP progfile, SEXP Verbose)
{
	int start = Rf_asInteger(Start);
	if (start < 1) start = 1;
	int count = Rf_asInteger(Count);
	const bool verbose = Rf_asLogical(Verbose)==TRUE;
	// variant selection
	const int *p_var_sel =
		Rf_isNull(R_variant_sel) ? NULL : INTEGER(R_variant_sel);

	COREARRAY_TRY

		// GDS nodes
		PdAbstractArray Root = GDS_R_SEXP2Obj(gds_root, FALSE);
		PdAbstractArray varGeno = GDS_Node_Path(Root, "genotype/data", FALSE);
		PdAbstractArray varGenoLen = GDS_Node_Path(Root, "genotype/@data", FALSE);
		PdAbstractArray varPhase = GDS_Node_Path(Root, "phase/data", FALSE);

		const size_t num_sample = Rf_length(R_buf);
		int *ptr_gt_buf = INTEGER(R_buf);
		vector<C_UInt8> zeros(num_sample);
		vector<C_Int8> gt(2*num_sample);
		const size_t gt_sz = gt.size();
		size_t nbit = 0;
		const C_Int8 NA = 0xFF;

		// progress information
		CProgress prog((verbose || !Rf_isNull(progfile)) ? count : -1, progfile);

		// for-loop
		for (size_t idx=start; count > 0; idx++, count--)
		{
			// call GetAlleleCt(pvar, ii)
			INTEGER(R_ii)[0] = p_var_sel ? p_var_sel[idx-1] : idx;
			int allele_cnt = Rf_asInteger(Rf_eval(R_allele_num_fc, R_env));
			if (allele_cnt > 127)
				throw "PLINK2 does not support > 127 alleles at a site.";
			// initialize all zeros
			memset(&gt[0], 0, gt.size());
			for (int allele_i=1; allele_i<allele_cnt; allele_i++)
			{
				// call ReadHardcalls(pgen, buf, ii, ia)
				INTEGER(R_ia)[0] = allele_i + 1;
				Rf_eval(R_read_gt_fc, R_env);
				// find non-zero
				for (size_t i=0; i < num_sample; i++)
				{
					switch(ptr_gt_buf[i])
					{
						case 0: break;
						case 1:
							{
								C_Int8 &g1 = gt[2*i+0], &g2 = gt[2*i+1];
								if (g1==0)
									g1 = allele_i;
								else if (g1!=NA && g2==0)
									g2 = allele_i;
								break;
							}
						case 2:
							gt[2*i+0] = gt[2*i+1] = allele_i; break;
						default:
							// NA, missing or invalid genotype
							gt[2*i+0] = gt[2*i+1] = NA;
					}
				}
			}

			// append genotypes
			C_UInt8 num_bits = 1;
			// using bit2 array
			bool flag = true;
			while(flag)
			{
				GDS_Array_AppendData(varGeno, gt_sz, &gt[0], svUInt8);
				nbit += gt_sz * 2;
				allele_cnt -= 4;
				flag = (allele_cnt>=0) || (count==1 && (nbit & 0x7));
				// the last one should be padded to a byte
				if (flag)
				{
					num_bits ++;
					for (size_t i=0; i < gt_sz; i++) gt[i] >>= 2;
				}
			}

			// append the number of bit2 rows
			GDS_Array_AppendData(varGenoLen, 1, &num_bits, svUInt8);
			// append phase data
			GDS_Array_AppendData(varPhase, num_sample, &zeros[0], svUInt8);

			// update progress
			prog.Forward(1);
		}

	COREARRAY_CATCH
}


// gds2pgen.cpp
extern SEXP SEQ_PGEN_Writer_Open(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP SEQ_PGEN_Writer_Block(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP SEQ_PGEN_Writer_Close(SEXP, SEXP);

/// initialize the package
COREARRAY_DLL_EXPORT void R_init_pgen2gds(DllInfo *info)
{
	#define CALL(name, num)	   { #name, (DL_FUNC)&name, num }
	static R_CallMethodDef callMethods[] =
	{
		CALL(SEQ_PGEN_Geno_Import, 12),
		CALL(SEQ_SetJobStatus, 2),
		CALL(SEQ_PGEN_Writer_Open, 5),
		CALL(SEQ_PGEN_Writer_Block, 5),
		CALL(SEQ_PGEN_Writer_Close, 2),
		{ NULL, NULL, 0 }
	};
	R_registerRoutines(info, NULL, callMethods, NULL, NULL);
	R_useDynamicSymbols(info, FALSE);
	Init_GDS_Routines();
}

} // extern "C"
