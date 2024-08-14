#ifndef __DALIGN_HPP
#define __DALIGN_HPP

#include "align.h"

#include <string>
#include <vector>

struct DalignData {
    Align_Spec*     align_spec;
    Work_Data*      work;
    float           basis_freq[4];
    Path            path;
    Alignment       align;
    std::vector<char>       aseq;
    std::vector<char>       bseq;
    double          error;
    double          ident_perc;

    DalignData(double _error) {
	    basis_freq[0] = .25;
	    basis_freq[1] = .25;
	    basis_freq[2] = .25;
	    basis_freq[3] = .25;
	    align_spec = New_Align_Spec(1.0 - _error, 100, basis_freq, 1);
        error = _error;
	    work = New_Work_Data();
    }

    ~DalignData() {
        Free_Align_Spec(align_spec);
	    Free_Work_Data(work);
    }
} ;

#define ocda_ident_perc(ocda)     ((ocda).ident_perc)
#define ocda_query_start(ocda)    ((ocda).path.abpos)
#define ocda_query_end(ocda)      ((ocda).path.aepos)
#define ocda_target_start(ocda)   ((ocda).path.bbpos)
#define ocda_target_end(ocda)     ((ocda).path.bepos)
#define ocda_distance(ocda)		  ((ocda).path.diffs)

int
dalign_align(DalignData* ocda,
	const u8* query,
    int qoff,
    int qsize,
	const u8* target,
    int toff,
    int tsize,
	const int min_align_size,
	const double min_ident_perc,
	int* qbeg,
	int* qend,
	int* tbeg,
	int* tend,
	double* ident_perc);

#endif // __DALIGN_HPP