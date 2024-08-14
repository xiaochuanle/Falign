#include "edlib_wrapper.hpp"

#include "edlib.h"
#include "hbn_traceback_aux.h"

using namespace std;

int
edlib_nw(EdlibAlignData* data,
	const u8* query,
    int query_size,
	const u8* target,
    int target_size,
    int tolerance,
    std::string& qaln,
    std::string& taln)
{
    qaln.clear();
    taln.clear();
    if (query_size == 0 || target_size == 0) return 1;

    data->query.clear();
    for (int i = 0; i < query_size; ++i) {
        int c = query[i];
        hbn_assert(c >= 0 && c < 4);
        data->query += DECODE_RESIDUE(c);
    }

    data->target.clear();
    for (int i = 0; i < target_size; ++i) {
        int c = target[i];
        hbn_assert(c >= 0 && c < 4);
        data->target += DECODE_RESIDUE(c);
    }
    EdlibAlignTask task = EDLIB_TASK_PATH;
    if (tolerance < 0) tolerance = hbn_max(query_size, target_size) * 0.5;
    EdlibAlignResult align = edlibAlign(data->query.c_str(),
                                query_size,
                                data->target.c_str(),
                                target_size,
                                edlibNewAlignConfig(tolerance, EDLIB_MODE_NW, task, NULL, 0));

    if (align.numLocations == 0) {
        tolerance = hbn_max(query_size, target_size);
        align = edlibAlign(data->query.c_str(),
                           query_size,
                           data->target.c_str(),
                           target_size,
                           edlibNewAlignConfig(tolerance, EDLIB_MODE_NW, task, NULL, 0));        
    }
    hbn_assert(align.numLocations);
    if (align.numLocations == 0) return 0;

    int qBgn = 0;
    int qEnd = query_size;
    int tBgn = align.startLocations[0];
    int tEnd = align.endLocations[0] + 1;
    data->vqaln.resize(align.alignmentLength+1);
    data->vtaln.resize(align.alignmentLength+1);
    edlibAlignmentToStrings(align.alignment,
        align.alignmentLength,
        tBgn,
        tEnd,
        qBgn,
        qEnd,
        data->target.c_str(),
        data->query.c_str(),
        data->vtaln.data(),
        data->vqaln.data());
    edlibFreeAlignResult(align);
    data->vqaln.pop_back();
    data->vtaln.pop_back();
    qaln.assign(data->vqaln.begin(), data->vqaln.end());
    taln.assign(data->vtaln.begin(), data->vtaln.end());
    hbn_assert(qaln.size() == taln.size());
    validate_aligned_string(HBN_LOG_ARGS_DEFAULT,
        0,
        query,
        qBgn,
        qEnd,
        qaln.c_str(),
        0,
        target,
        tBgn,
        tEnd,
        taln.c_str(),
        qaln.size(),
        TRUE);

    return 1;
}

int
edlib_nw_dist(EdlibAlignData* data,
	const u8* query,
    int query_size,
	const u8* target,
    int target_size,
    int tolerance,
    int* score)
{
    if (query_size == 0 || target_size == 0) return 1;

    data->query.clear();
    for (int i = 0; i < query_size; ++i) {
        int c = query[i];
        hbn_assert(c >= 0 && c < 4);
        data->query += DECODE_RESIDUE(c);
    }
    data->target.clear();
    for (int i = 0; i < target_size; ++i) {
        int c = target[i];
        hbn_assert(c >= 0 && c < 4);
        data->target += DECODE_RESIDUE(c);
    }
    EdlibAlignTask task = EDLIB_TASK_DISTANCE;
    if (tolerance < 0) tolerance = hbn_max(query_size, target_size) * 0.5;
    EdlibAlignResult align = edlibAlign(data->query.c_str(),
                                query_size,
                                data->target.c_str(),
                                target_size,
                                edlibNewAlignConfig(tolerance, EDLIB_MODE_NW, task, NULL, 0));

    if (align.numLocations == 0) {
        tolerance = hbn_max(query_size, target_size);
        align = edlibAlign(data->query.c_str(),
                           query_size,
                           data->target.c_str(),
                           target_size,
                           edlibNewAlignConfig(tolerance, EDLIB_MODE_NW, task, NULL, 0));        
    }
    hbn_assert(align.numLocations);
    edlibFreeAlignResult(align);
    if (align.numLocations == 0) return -1;
    if (score) {
        *score = ((query_size + target_size) / 2 - align.editDistance) * 2 - align.editDistance * 5;
    }
    return align.editDistance;
}

int
edlib_shw(EdlibAlignData* data,
    const u8* _query,
    int query_size,
    const u8* _target,
    int target_size,
    int* qend,
    int* tend,
    std::string& qaln,
    std::string& taln)
{
    qaln.clear();
    taln.clear();
    *qend = *tend = 0;
    if (query_size == 0 || target_size == 0) return 0;

    data->query.clear();
    for (int i = 0; i < query_size; ++i) {
        int c = _query[i];
        hbn_assert(c >= 0 && c < 4, "qid = %d, sid = %d, i = %d, c = %d, query_size = %d, target_size = %d", data->qid, data->sid, i, c, query_size, target_size);
        data->query += DECODE_RESIDUE(c);
    }
    data->target.clear();
    for (int i = 0; i < target_size; ++i) {
        int c = _target[i];
        hbn_assert(c >= 0 && c < 4);
        data->target += DECODE_RESIDUE(c);
    }
    EdlibAlignTask task = EDLIB_TASK_PATH;
    int tolerance = hbn_max(query_size, target_size) * 0.35;
    EdlibAlignResult align = edlibAlign(data->query.c_str(),
                                query_size,
                                data->target.c_str(),
                                target_size,
                                edlibNewAlignConfig(tolerance, EDLIB_MODE_SHW, task, NULL, 0));
    if (align.numLocations == 0) return 0;

    int qBgn = 0;
    int qEnd = query_size;
    int tBgn = align.startLocations[0];
    int tEnd = align.endLocations[0] + 1;
    hbn_assert(tBgn == 0);
    data->vqaln.resize(align.alignmentLength+1);
    data->vtaln.resize(align.alignmentLength+1);
    edlibAlignmentToStrings(align.alignment,
        align.alignmentLength,
        tBgn,
        tEnd,
        qBgn,
        qEnd,
        data->target.c_str(),
        data->query.c_str(),
        data->vtaln.data(),
        data->vqaln.data());
    edlibFreeAlignResult(align);
    data->vqaln.pop_back();
    data->vtaln.pop_back();
    qaln.assign(data->vqaln.begin(), data->vqaln.end());
    taln.assign(data->vtaln.begin(), data->vtaln.end());
    hbn_assert(qaln.size() == taln.size());
    *qend = qEnd;
    *tend = tEnd;

    validate_aligned_string(HBN_LOG_ARGS_DEFAULT,
        0,
        _query,
        qBgn,
        qEnd,
        qaln.c_str(),
        0,
        _target,
        tBgn,
        tEnd,
        taln.c_str(),
        qaln.size(),
        TRUE);

    return 1;    
}

static BOOL
get_next_sequence_block(const u8* query,
						int qidx,
						const int qsize,
						const u8* target,
						int tidx,
						const int tsize,
						const int desired_block_size,
						const BOOL right_extend,
						vector<u8>& qfrag,
						vector<u8>& tfrag)
{
	BOOL last_block = FALSE;
	int qleft = qsize - qidx;
	int tleft = tsize - tidx;
	int qblk;
	int tblk;
	if (qleft < desired_block_size + 100 || tleft < desired_block_size + 100) {
		qblk = tleft * 1.3;
		qblk = hbn_min(qblk, qleft);
		tblk = qleft * 1.3;
		tblk = hbn_min(tblk, tleft);
		last_block = TRUE;
	} else {
		qblk = desired_block_size;
		tblk = desired_block_size;
		last_block = FALSE;
	}
	
    qfrag.clear();
    tfrag.clear();
	if (right_extend) {
		const u8* Q = query + qidx;
		for (int i = 0; i < qblk; ++i) qfrag.push_back(Q[i]);
		const u8* R = target + tidx;
		for (int i = 0; i < tblk; ++i) tfrag.push_back(R[i]);
	} else {
		const u8* Q = query - qidx;
		for (int i = 0; i < qblk; ++i) qfrag.push_back(Q[-i]);
		const u8* R = target - tidx;
		for (int i = 0; i < tblk; ++i) tfrag.push_back(R[-i]);
	}
	
	return last_block;
}

void
edlib_extend(EdlibAlignData* data,
    const u8* query,
    const int query_size,
    const u8* target,
    const int target_size,
    const int block_size,
    const BOOL right_extend,
    std::vector<u8>& qfrag,
    std::vector<u8>& tfrag,
    int* qend,
    int* tend,
    std::string& qaln,
    std::string& taln)
{
    qaln.clear();
    taln.clear();
    int qidx = 0, tidx = 0;
    while (1) {
        int qfae, tfae, qfrag_size, tfrag_size;
        BOOL last_block = get_next_sequence_block(query,
                            qidx,
                            query_size,
                            target,
                            tidx,
                            target_size,
                            block_size,
                            right_extend,
                            qfrag,
                            tfrag);
        qfrag_size = qfrag.size();
        tfrag_size = tfrag.size();
        //fprintf(stderr, "qsize = %d, ssize = %d\n", qfrag_size, tfrag_size);
        if (qfrag_size == 0 || tfrag_size == 0) break;

        edlib_shw(data,
            qfrag.data(),
            qfrag_size,
            tfrag.data(),
            tfrag_size,
            &qfae,
            &tfae,
            data->sqaln,
            data->staln);

        BOOL done = last_block;
        if (qfrag_size - qfae > 100 || tfrag_size - tfae > 100) done = TRUE;
        int acnt = 0, qcnt = 0, tcnt = 0;
        const int M = 8;
        int align_size = data->sqaln.size();
        int k = align_size - 1, m = 0;
        while (k >= 0) {
            const char qc = data->sqaln[k];
            const char tc = data->staln[k];
            if (qc != GAP_CHAR) ++qcnt;
            if (tc != GAP_CHAR) ++tcnt;
            m = (qc == tc) ? (m+1) : 0;
            ++acnt;
            if (m == M) break;
            --k;
        }

        if (m != M || k < 1) {
            align_size = 0;
            data->sqaln.clear();
            data->staln.clear();
            for (int i = 0; i < qfrag_size && i < tfrag_size; ++i) {
                int qc = qfrag[i];
                int tc = tfrag[i];
                if (qc != tc) break;
                qc = DECODE_RESIDUE(qc);
                tc = DECODE_RESIDUE(tc);
                data->sqaln += qc;
                data->staln += tc;
                ++align_size;
            }
            done = TRUE;
        } else {
            align_size -= acnt;
            qidx += (qfae - qcnt);
            tidx += (tfae - tcnt);
            if (done) align_size += M;
        }
       // HBN_LOG("qidx = %d/%d, tidx = %d/%d", qidx, query_size, tidx, target_size);

        qaln.append(data->sqaln.c_str(), align_size);
        taln.append(data->staln.c_str(), align_size);
        if (done) break;
    }

	int qe = 0, te = 0;
	for (size_t i = 0; i != qaln.size(); ++i) {
		if (qaln[i] != GAP_CHAR) ++qe;
		if (taln[i] != GAP_CHAR) ++te;
	}
    validate_aligned_string(HBN_LOG_ARGS_DEFAULT,
        0,
        query,
        0,
        qe,
        qaln.c_str(),
        0,
        target,
        0,
        te,
        taln.c_str(),
        qaln.size(),
        right_extend);
    *qend = qe;
    *tend = te;
}