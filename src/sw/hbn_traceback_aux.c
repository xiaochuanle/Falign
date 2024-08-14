#include "hbn_traceback_aux.h"

#include <stdlib.h>

void extract_sw_scoring_params(int* match_reward, int* mismatch_penalty, int* gap_open, int* gap_extend, int* gap_open1, int* gap_extend1)
{
#if 0
    data->reward = 2;
    data->penalty = 4;
    data->go = 4;
    data->ge = 2;
    data->go1 = 24;
    data->ge1 = 1;
#elif 0
    //-o 5 -O 56 -e 4 -E 1 -A 2 -B 5
    data->reward = 2;
    data->penalty = 5;
    data->go = 5;
    data->ge = 4;
    data->go1 = 56;
    data->ge1 = 1;
#endif 
    if (match_reward) *match_reward = 2;
    if (mismatch_penalty) *mismatch_penalty = 4;
    if (gap_open) *gap_open = 4;
    if (gap_extend) *gap_extend = 2;
    if (gap_open1) *gap_open1 = 24;
    if (gap_extend1) *gap_extend1 = 1;
}

double
calc_ident_perc(const char* query_mapped_string, 
				const char* target_mapped_string,
			    const int align_size,
                int* dist,
				int* score)
{
	if (align_size == 0) return 0.0;
	
	int n = 0;
	for (int i = 0; i < align_size; ++i) {
		if (query_mapped_string[i] == target_mapped_string[i]) ++n;
	}
    if (dist) *dist = align_size - n;

	if (score) {
		int a, b, go, ge;
		extract_sw_scoring_params(&a, &b, &go, &ge, NULL, NULL);
		int i = 0;
		int S = 0;
		while (i < align_size) {
			int qc = query_mapped_string[i];
			int tc = target_mapped_string[i];
			if (qc != GAP_CHAR && tc != GAP_CHAR) {
				if (qc == tc) {
					S += a;
				} else {
					S -= b;
				}
				++i;
				continue;
			}
			if (qc == GAP_CHAR) {
				int j = i + 1;
				while (j < align_size && query_mapped_string[j] == GAP_CHAR) ++j;
				S += -(ge * (j - i));
				i = j;
				continue;
			}
			if (tc == GAP_CHAR) {
				int j = i + 1;
				while (j < align_size && target_mapped_string[j] == GAP_CHAR) ++j;
				S += -(ge * (j - i));
				i = j;
				continue;
			}
		}
		*score = S;
	}

	return 100.0 * n / align_size;
}

double
calc_effective_ident_perc(const char* query_mapped_string, 
				const char* target_mapped_string,
			    const int align_size)
{
	const int E = 20;
	int eff_len = 0;
	int eff_mat = 0;
	int i = 0;
	while (i < align_size) {
		int qc = query_mapped_string[i];
		int tc = target_mapped_string[i];

		if (qc != GAP_CHAR && tc != GAP_CHAR) {
			if (qc == tc) ++eff_mat;
			++eff_len;
			++i;
			continue;
		}

		if (qc == GAP_CHAR && tc == GAP_CHAR) {
			++i;
			continue;
		}

		if (qc == GAP_CHAR) {
			int j = i + 1;
			while (j < align_size) {
				if (query_mapped_string[j] == GAP_CHAR && target_mapped_string[j] == GAP_CHAR) {
					++j;
					continue;
				}
				if (query_mapped_string[j] != GAP_CHAR) break;
				++j;
			}
			if (j - i < E) {
				for (int k = i; k < j; ++k) {
					int qc = query_mapped_string[k];
					int tc = target_mapped_string[k];
					if (qc == GAP_CHAR && tc == GAP_CHAR) continue;
					if (qc == tc) ++eff_mat;
					++eff_len;
				}
			}
			i = j;
			continue;
		}

		hbn_assert(tc == GAP_CHAR);
		if (tc == GAP_CHAR) {
			int j = i + 1;
			while (j < align_size) {
				if (query_mapped_string[j] == GAP_CHAR && target_mapped_string[j] == GAP_CHAR) {
					++j;
					continue;
				}
				if (target_mapped_string[j] != GAP_CHAR) break;
				++j;
			}
			if (j - i < E) {
				for (int k = i; k < j; ++k) {
					qc = query_mapped_string[k];
					tc = target_mapped_string[k];
					if (qc == GAP_CHAR && tc == GAP_CHAR) continue;
					if (qc == tc) ++eff_mat;
					++eff_len;
				}
			}
			i = j;
			continue;
		}
	}
	if (eff_len == 0) return 0.0;
	return 100.0 * eff_mat / eff_len;
}

void
validate_aligned_string(const char* source_file,
						const char* source_func,
						const int source_line,
						int qid,
						const u8* query,
						const int qoff,
						const int qend,
						const char* query_mapped_string,
						int tid,
						const u8* target,
						const int toff,
						const int tend,
						const char* target_mapped_string,
						const size_t align_size,
					    const BOOL right_extend)
{
	//return;
	int x = qoff, y = toff;
	for (size_t i = 0; i != align_size; ++i) {
		const char qc = query_mapped_string[i];
		if (qc != GAP_CHAR) {
            int c = right_extend ? query[x] : query[-x];
			const char qc1 = DECODE_RESIDUE(c);
            if (qc != qc1) {
			    fprintf(stderr, "[%s, %s, %d] qid = %d, tid = %d, right_extend = %d, i = %lu, x = %d, y = %d, qc = %c, qc1 = %c, qoff = %d, qend = %d, toff = %d, tend = %d, align_size = %lu\n",
					  source_file,
					  source_func,
					  source_line,
					  qid,
					  tid,
					  right_extend,
					  i,
					  x,
					  y,
					  qc,
					  qc1,
					  qoff,
					  qend,
					  toff,
					  tend,
					  align_size);
                abort();
            }		  
			++x;
		}
		const char tc = target_mapped_string[i];
		if (tc != GAP_CHAR) {
            int c = right_extend ? target[y] : target[-y];
			const char tc1 = DECODE_RESIDUE(c);
            if (tc != tc1) {
			    fprintf(stderr, "[%s, %s, %d] qid = %d, tid = %d, right_extend = %d, i = %lu, x = %d, y = %d, tc = %c, tc1 = %c, qoff = %d, qend = %d, toff = %d, tend = %d\n",
						  source_func,
						  source_func,
						  source_line,
						  qid,
						  tid,
						  right_extend,
						  i,
						  x,
						  y,
						  tc,
						  tc1,
						  qoff,
						  qend,
						  toff,
						  tend);
                dump_align_string(query_mapped_string, target_mapped_string, align_size, stderr);
                abort();
            }
			++y;
		}
	}
}
