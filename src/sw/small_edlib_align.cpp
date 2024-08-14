#include "small_edlib_align.hpp"

#include "../corelib/hbn_aux.h"

#include <math.h>
#include <algorithm>

using namespace std;

typedef u64             Word;
#define WORD_SIZE       64
#define WORD_ONE        ((Word)1)
#define HIGH_BIT_MASK   (WORD_ONE << (WORD_SIZE - 1))

//#define MaxSeqSize		4096
#define MaxSeqSize      SMALL_EDLIB_MAX_SEQ_SIZE
#define MaxNumBlocks	64
#define AlphabetSize    4

#define EDLIB_OP_SUB 	('M')
#define EDLIB_OP_INS 	('I')
#define EDLIB_OP_DEL 	('D')

/**
 * Describes cigar format.
 * @see http://samtools.github.io/hts-specs/SAMv1.pdf
 * @see http://drive5.com/usearch/manual/cigar.html
 */
typedef enum {
	EEDLIB_CIGAR_STANDARD,  //!< Match: 'M', Insertion: 'I', Deletion: 'D', Mismatch: 'M'.
	EEDLIB_CIGAR_EXTENDED   //!< Match: '=', Insertion: 'I', Deletion: 'D', Mismatch: 'X'.
} SmallEdlibCigarFormat;

typedef enum {
    EEDLIB_MODE_NW,
    EEDLIB_MODE_SHW,
    EEDLIB_MODE_HW
} SmallEdlibAlignMode;

typedef struct {
	Word* Ps;
	Word* Ms;
	int* scores;
	int* first_blocks;
	int* last_blocks;
} EdlibAlignMatrix;

EdlibAlignMatrix*
EdlibAlignMatrixNew(int maxNumBlocks, int targetLength)
{
    EdlibAlignMatrix* data = (EdlibAlignMatrix*)malloc(sizeof(EdlibAlignMatrix));
    data->Ps = (Word*)malloc(sizeof(Word) * maxNumBlocks * targetLength);
    data->Ms = (Word*)malloc(sizeof(Word) * maxNumBlocks * targetLength);
    data->scores = (int*)malloc(sizeof(int) * maxNumBlocks * targetLength);
    data->first_blocks = (int*)malloc(sizeof(int) * targetLength);
    data->last_blocks = (int*)malloc(sizeof(int) * targetLength);

    return data;
}

EdlibAlignMatrix*
EdlibAlignMatrixFree(EdlibAlignMatrix* data)
{
    free(data->Ps);
    free(data->Ms);
    free(data->scores);
    free(data->first_blocks);
    free(data->last_blocks);
    free(data);
    return 0;
}

typedef struct {
	Word P;
	Word M;
	int score;
} EdlibBlock;

struct SmallEdlibAlignResult {
	int 		edit_distance;
    vector<int> end_locations;
    vector<int> start_locations;
    vector<u8> alignment;

    void clear() {
        edit_distance = -1;
        end_locations.clear();
        start_locations.clear();
        alignment.clear();
    }
} ;

typedef struct {
	int num;
	int op;
} EdlibGapAlignOp;

struct SmallEdlibAlignData {
    EdlibAlignMatrix*   align_matrix;
    SmallEdlibAlignResult   result;
    EdlibBlock*         blocks;
    Word*               peq;
    vector<EdlibGapAlignOp>      cigar;
	int					maxNumBlocks;
	int					maxTargetLength;

    SmallEdlibAlignData(int maxNumBlocks, int targetLength) {
        align_matrix = EdlibAlignMatrixNew(maxNumBlocks, targetLength);
        blocks = (EdlibBlock*)malloc(sizeof(EdlibBlock) * maxNumBlocks);
        peq = (Word*)malloc(sizeof(Word) * (AlphabetSize + 1) * maxNumBlocks);
	    this->maxNumBlocks = maxNumBlocks;
	    this->maxTargetLength = targetLength;
    }

    ~SmallEdlibAlignData() {
        EdlibAlignMatrixFree(align_matrix);
        free(blocks);
	    free(peq);
    }
};

void
SmallEdlibAlignDataAllocSpace(SmallEdlibAlignData* data, int query_length, int target_length)
{
	int newMaxNumBlocks = (query_length + WORD_SIZE - 1) / WORD_SIZE;
	int newMaxTargetLength = target_length;
	int newMaxSpace = sizeof(Word) * newMaxNumBlocks * newMaxTargetLength;
	
	if (newMaxNumBlocks > data->maxNumBlocks || newMaxTargetLength > data->maxTargetLength) {
		data->align_matrix->Ps = (Word*)realloc(data->align_matrix->Ps, newMaxSpace);
		data->align_matrix->Ms = (Word*)realloc(data->align_matrix->Ms, newMaxSpace);
		data->align_matrix->scores = (int*)realloc(data->align_matrix->scores, sizeof(int) * newMaxNumBlocks * newMaxTargetLength);
		data->align_matrix->first_blocks = (int*)realloc(data->align_matrix->first_blocks, sizeof(int) * newMaxTargetLength);
		data->align_matrix->last_blocks = (int*)realloc(data->align_matrix->last_blocks, sizeof(int) * newMaxTargetLength);
		
		data->blocks = (EdlibBlock*)realloc(data->blocks, sizeof(EdlibBlock) * newMaxNumBlocks);
		data->peq = (Word*)realloc(data->peq, sizeof(Word) * (AlphabetSize + 1) * newMaxNumBlocks);
		
		data->maxNumBlocks = newMaxNumBlocks;
		data->maxTargetLength = newMaxTargetLength;
	}
}

/////////////////////

#define EDLIB_STATUS_OK     0
#define EDLIB_STATUS_ERROR  0

// Edit operations.
#define EDLIB_EDOP_MATCH 0    //!< Match.
#define EDLIB_EDOP_INSERT 1   //!< Insertion to target = deletion from query.
#define EDLIB_EDOP_DELETE 2   //!< Deletion from target = insertion to query.
#define EDLIB_EDOP_MISMATCH 3 //!< Mismatch.

static inline int
calc_num_blocks(const int n)
{
    return (n + WORD_SIZE - 1) / WORD_SIZE;
}

static inline void
calc_block_cell_scores(const EdlibBlock block, int scores[])
{
    int score = block.score;
    Word mask = HIGH_BIT_MASK;
    for (int i = 0; i < WORD_SIZE - 1; ++i) {
        scores[i] = score;
        if (block.P & mask) --score;
        if (block.M & mask) ++score;
        mask >>= 1;
    }
    scores[WORD_SIZE - 1] = score;
}

static void
build_peq(const u8* query, const int query_size, Word* peq)
{
    const int nblk = calc_num_blocks(query_size);
    for (int s = 0; s <= AlphabetSize; ++s) {
        for (int b = 0; b < nblk; ++b) {
            const int bid = s * MaxNumBlocks + b;
            if (s < AlphabetSize) {
                peq[bid] = 0;
                for (int r = (b + 1) * WORD_SIZE - 1; r >= b * WORD_SIZE; --r) {
                    peq[bid] <<= 1;
                    if (r >= query_size || query[r] == s) peq[bid] += 1;
                }
            } else {
                peq[bid] = (Word)-1;
            }
        }
    }
}

/**
 * Corresponds to Advance_Block function from Myers.
 * Calculates one word(block), which is part of a column.
 * Highest bit of word (one most to the left) is most bottom cell of block from column.
 * Pv[i] and Mv[i] define vin of cell[i]: vin = cell[i] - cell[i-1].
 * @param [in] Pv  Bitset, Pv[i] == 1 if vin is +1, otherwise Pv[i] == 0.
 * @param [in] Mv  Bitset, Mv[i] == 1 if vin is -1, otherwise Mv[i] == 0.
 * @param [in] Eq  Bitset, Eq[i] == 1 if match, 0 if mismatch.
 * @param [in] hin  Will be +1, 0 or -1.
 * @param [out] PvOut  Bitset, PvOut[i] == 1 if vout is +1, otherwise PvOut[i] == 0.
 * @param [out] MvOut  Bitset, MvOut[i] == 1 if vout is -1, otherwise MvOut[i] == 0.
 * @param [out] hout  Will be +1, 0 or -1.
 */
static inline int calculateBlock(Word Pv, Word Mv, Word Eq, const int hin,
                                 Word* PvOut, Word* MvOut) {
    // hin can be 1, -1 or 0.
    // 1  -> 00...01
    // 0  -> 00...00
    // -1 -> 11...11 (2-complement)

    Word hinIsNeg = (Word)(hin >> 2) & WORD_ONE; // 00...001 if hin is -1, 00...000 if 0 or 1

    Word Xv = Eq | Mv;
    // This is instruction below written using 'if': if (hin < 0) Eq |= (Word)1;
    Eq |= hinIsNeg;
    Word Xh = (((Eq & Pv) + Pv) ^ Pv) | Eq;

    Word Ph = Mv | ~(Xh | Pv);
    Word Mh = Pv & Xh;

    int hout = 0;
    // This is instruction below written using 'if': if (Ph & HIGH_BIT_MASK) hout = 1;
    hout = (Ph & HIGH_BIT_MASK) >> (WORD_SIZE - 1);
    // This is instruction below written using 'if': if (Mh & HIGH_BIT_MASK) hout = -1;
    hout -= (Mh & HIGH_BIT_MASK) >> (WORD_SIZE - 1);

    Ph <<= 1;
    Mh <<= 1;

    // This is instruction below written using 'if': if (hin < 0) Mh |= (Word)1;
    Mh |= hinIsNeg;
    // This is instruction below written using 'if': if (hin > 0) Ph |= (Word)1;
    Ph |= (Word)((hin + 1) >> 1);

    (*PvOut) = Mh | ~(Xv | Ph);
    (*MvOut) = Ph & Xv;

    return hout;
}

static void
calc_edit_distance_nw(const u8* query, const int query_size,
					  const u8* target, const int target_size,
					  const Word* peq,
					  int k,
					  const BOOL traceback,
					  EdlibAlignMatrix* align_matrix,
					  EdlibBlock* blocks,
					  const int target_stop_position,
					  int* edit_distance,
					  int* end_position)
{
    *edit_distance = -1;
    *end_position = -1;

    if (target_stop_position > -1 && traceback) {
		HBN_ERR("invalid parameters: target_stop_position = %d, traceback = %d", target_stop_position, traceback);
    }

    if (k < abs(target_size - query_size)) {
        return;
    }

    //k = min(k, max(query_size, target_size));
	{
		int X = hbn_max(query_size, target_size);
		k = hbn_min(k, X);
	}

    const int nblk = calc_num_blocks(query_size);
    const int W = nblk * WORD_SIZE - query_size;
    int fblk = 0;
    int lblk; // = min(nblk, calc_num_blocks(min(k, (k + query_size - target_size) / 2) + 1)) - 1;
	{
		int X = (k + query_size - target_size) / 2;
		int Y = hbn_min(k, X);
		int Z = calc_num_blocks(Y + 1);
		lblk = hbn_min(nblk, Z) - 1;
	}
    EdlibBlock* blk = blocks;

    for (int b = 0; b <= lblk; ++b) {
        blk->score = (b + 1) * WORD_SIZE;
        blk->P = (Word)-1;
        blk->M = (Word)0;
        ++blk;
    }

    for (int c = 0; c < target_size; ++c) {
        const int tc = target[c];
        const Word* cpeq = peq + tc * MaxNumBlocks;
        int hout = 1;
        blk = blocks + fblk;
        for (int b = fblk; b <= lblk; ++b) {
            hout = calculateBlock(blk->P, blk->M, cpeq[b], hout, &blk->P, &blk->M);
            blk->score += hout;
            ++blk;
        }
        --blk;

        //k = min(k, blk->score
        //        + max(target_size -c - 1, query_size - ((1 + lblk) * WORD_SIZE - 1) - 1)
        //        + (lblk == nblk - 1 ? W : 0));
		
		{
			int X1 = target_size - c - 1, X2 = query_size - ((1 + lblk) * WORD_SIZE - 1) - 1;
			int X = hbn_max(X1, X2);
			int Y = (lblk == nblk - 1) ? W : 0;
			int Z = X + Y + blk->score;
			k = hbn_min(k, Z);
		}

        if (lblk + 1 < nblk) {
            BOOL r = (lblk + 1) * WORD_SIZE - 1 > k - blk->score + 2 * WORD_SIZE - 2 - target_size + c + query_size;
            if (!r) {
                ++lblk;
                ++blk;
                blk->P = (Word)-1;
                blk->M = (Word)0;
                int newHout = calculateBlock(blk->P, blk->M, cpeq[lblk], hout, &blk->P, &blk->M);
                blk->score = (blk - 1)->score - hout + WORD_SIZE + newHout;
                hout = newHout;
            }
        }

        while (lblk >= fblk
				&&
				(blk->score >= k + WORD_SIZE
				 || ((lblk + 1) * WORD_SIZE -1 > 
					 k - blk->score + 2 * WORD_SIZE - 2 - target_size + c + query_size + 1))) {
            --lblk;
            --blk;
        }

        while (fblk <= lblk
				&&
				(blocks[fblk].score >= k + WORD_SIZE
				 || ((fblk + 1) * WORD_SIZE - 1 <
					 blocks[fblk].score - k - target_size + query_size + c))) {
            ++fblk;
        }

        if (lblk < fblk) {
            *edit_distance = *end_position = -1;
            return;
        }

        if (traceback) {
            blk = blocks + fblk;
            for (int b = fblk; b <= lblk; ++b) {
                align_matrix->Ps[MaxNumBlocks * c + b] = blk->P;
                align_matrix->Ms[MaxNumBlocks * c + b] = blk->M;
                align_matrix->scores[MaxNumBlocks * c + b] = blk->score;
                align_matrix->first_blocks[c] = fblk;
                align_matrix->last_blocks[c] = lblk;
                ++blk;
            }
        }

        if (c == target_stop_position) {
            for (int b = fblk; b <= lblk; ++b) {
                align_matrix->Ps[b] = blocks[b].P;
                align_matrix->Ms[b] = blocks[b].M;
                align_matrix->scores[b] = blocks[b].score;
                align_matrix->first_blocks[0] = fblk;
                align_matrix->last_blocks[0] = lblk;
            }
            *edit_distance = -1;
            *end_position = target_stop_position;
            return;
        }
    }

    if (lblk == nblk - 1) {
        int scores[WORD_SIZE];
        calc_block_cell_scores(blocks[lblk], scores);
        int col_score = scores[W];
        if (col_score <= k) {
            *edit_distance = col_score;
            *end_position = target_size - 1;
            return;
        }
    }

    *edit_distance = *end_position = -1;
}

/**
 * Finds one possible alignment that gives optimal score by moving back through the dynamic programming matrix,
 * that is stored in alignData. Consumes large amount of memory: O(queryLength * targetLength).
 * @param [in] queryLength  Normal length, without W.
 * @param [in] targetLength  Normal length, without W.
 * @param [in] bestScore  Best score.
 * @param [in] alignData  Data obtained during finding best score that is useful for finding alignment.
 * @param [out] alignment  Alignment.
 * @param [out] alignmentLength  Length of alignment.
 * @return Status code.
 */
static int obtainAlignmentTraceback(const int queryLength, 
									const int targetLength,
                                    const int bestScore, //const AlignmentData* const alignData,
									const EdlibAlignMatrix* align_matrix,
									vector<u8>& alignment) {
                                    //vector<unsigned char>& alignment) {
	const int nblk = calc_num_blocks(queryLength);
	const int W = nblk * WORD_SIZE - queryLength;
	alignment.clear();
    int c = targetLength - 1; // index of column
    int b = nblk - 1; // index of block in column
    int currScore = bestScore; // Score of current cell
    int lScore  = -1; // Score of left cell
    int uScore  = -1; // Score of upper cell
    int ulScore = -1; // Score of upper left cell
    Word currP = align_matrix->Ps[c * MaxNumBlocks + b]; // P of current block
    Word currM = align_matrix->Ms[c * MaxNumBlocks + b]; // M of current block
    // True if block to left exists and is in band
    BOOL thereIsLeftBlock = c > 0 && b >= align_matrix->first_blocks[c-1] && b <= align_matrix->last_blocks[c-1];
    // We set initial values of lP and lM to 0 only to avoid compiler warnings, they should not affect the
    // calculation as both lP and lM should be initialized at some moment later (but compiler can not
    // detect it since this initialization is guaranteed by "business" logic).
    Word lP = 0, lM = 0;
    if (thereIsLeftBlock) {
        lP = align_matrix->Ps[(c - 1) * MaxNumBlocks + b]; // P of block to the left
        lM = align_matrix->Ms[(c - 1) * MaxNumBlocks + b]; // M of block to the left
    }
    currP <<= W;
    currM <<= W;
    int blockPos = WORD_SIZE - W - 1; // 0 based index of current cell in blockPos

    // TODO(martin): refactor this whole piece of code. There are too many if-else statements,
    // it is too easy for a bug to hide and to hard to effectively cover all the edge-cases.
    // We need better separation of logic and responsibilities.
    while (1) {
        if (c == 0) {
            thereIsLeftBlock = TRUE;
            lScore = b * WORD_SIZE + blockPos + 1;
            ulScore = lScore - 1;
        }

        // TODO: improvement: calculate only those cells that are needed,
        //       for example if I calculate upper cell and can move up,
        //       there is no need to calculate left and upper left cell
        //---------- Calculate scores ---------//
        if (lScore == -1 && thereIsLeftBlock) {
            lScore = align_matrix->scores[(c - 1) * MaxNumBlocks + b]; // score of block to the left
            for (int i = 0; i < WORD_SIZE - blockPos - 1; i++) {
                if (lP & HIGH_BIT_MASK) lScore--;
                if (lM & HIGH_BIT_MASK) lScore++;
                lP <<= 1;
                lM <<= 1;
            }
        }
        if (ulScore == -1) {
            if (lScore != -1) {
                ulScore = lScore;
                if (lP & HIGH_BIT_MASK) ulScore--;
                if (lM & HIGH_BIT_MASK) ulScore++;
            }
            else if (c > 0 && b-1 >= align_matrix->first_blocks[c-1] && b-1 <= align_matrix->last_blocks[c-1]) {
                // This is the case when upper left cell is last cell in block,
                // and block to left is not in band so lScore is -1.
                ulScore = align_matrix->scores[(c - 1) * MaxNumBlocks + b - 1];
            }
        }
        if (uScore == -1) {
            uScore = currScore;
            if (currP & HIGH_BIT_MASK) uScore--;
            if (currM & HIGH_BIT_MASK) uScore++;
            currP <<= 1;
            currM <<= 1;
        }
        //-------------------------------------//

        // TODO: should I check if there is upper block?

        //-------------- Move --------------//
        // Move up - insertion to target - deletion from query
        if (uScore != -1 && uScore + 1 == currScore) {
            currScore = uScore;
            lScore = ulScore;
            uScore = ulScore = -1;
            if (blockPos == 0) { // If entering new (upper) block
                if (b == 0) { // If there are no cells above (only boundary cells)
					//alignment.push_back(EDLIB_EDOP_INSERT);
					//for (int i = 0; i < c + 1; ++i) alignment.push_back(EDLIB_EDOP_DELETE);
                    alignment.push_back(EDLIB_EDOP_INSERT);
					for (int i = 0; i < c + 1; ++i) alignment.push_back(EDLIB_EDOP_DELETE);
					break;
                } else {
                    blockPos = WORD_SIZE - 1;
                    b--;
                    currP = align_matrix->Ps[c * MaxNumBlocks + b];
                    currM = align_matrix->Ms[c * MaxNumBlocks + b];
                    if (c > 0 && b >= align_matrix->first_blocks[c-1] && b <= align_matrix->last_blocks[c-1]) {
                        thereIsLeftBlock = TRUE;
                        lP = align_matrix->Ps[(c - 1) * MaxNumBlocks + b]; // TODO: improve this, too many operations
                        lM = align_matrix->Ms[(c - 1) * MaxNumBlocks + b];
                    } else {
                        thereIsLeftBlock = FALSE;
                        // TODO(martin): There may not be left block, but there can be left boundary - do we
                        // handle this correctly then? Are l and ul score set correctly? I should check that / refactor this.
                    }
                }
            } else {
                blockPos--;
                lP <<= 1;
                lM <<= 1;
            }
            // Mark move
            ///-(*alignment)[(*alignmentLength)++] = EDLIB_EDOP_INSERT;
			alignment.push_back(EDLIB_EDOP_INSERT);
        }
        // Move left - deletion from target - insertion to query
        else if (lScore != -1 && lScore + 1 == currScore) {
            currScore = lScore;
            uScore = ulScore;
            lScore = ulScore = -1;
            c--;
            if (c == -1) { // If there are no cells to the left (only boundary cells)
                ///(*alignment)[(*alignmentLength)++] = EDLIB_EDOP_DELETE; // Move left
                ///int numUp = b * WORD_SIZE + blockPos + 1;
                ///for (int i = 0; i < numUp; i++) // Move up until end
                ///    (*alignment)[(*alignmentLength)++] = EDLIB_EDOP_INSERT;
                ///break;
				
				alignment.push_back(EDLIB_EDOP_DELETE);
				int numUp = b * WORD_SIZE + blockPos + 1;
				for (int i = 0; i < numUp; i++) alignment.push_back(EDLIB_EDOP_INSERT);
				break;
            }
            currP = lP;
            currM = lM;
            if (c > 0 && b >= align_matrix->first_blocks[c-1] && b <= align_matrix->last_blocks[c-1]) {
                thereIsLeftBlock = TRUE;
                lP = align_matrix->Ps[(c - 1) * MaxNumBlocks + b];
                lM = align_matrix->Ms[(c - 1) * MaxNumBlocks + b];
            } else {
                if (c == 0) { // If there are no cells to the left (only boundary cells)
                    thereIsLeftBlock = TRUE;
                    lScore = b * WORD_SIZE + blockPos + 1;
                    ulScore = lScore - 1;
                } else {
                    thereIsLeftBlock = FALSE;
                }
            }
            // Mark move
            ///(*alignment)[(*alignmentLength)++] = EDLIB_EDOP_DELETE;
			alignment.push_back(EDLIB_EDOP_DELETE);
        }
        // Move up left - (mis)match
        else if (ulScore != -1) {
            unsigned char moveCode = ulScore == currScore ? EDLIB_EDOP_MATCH : EDLIB_EDOP_MISMATCH;
            currScore = ulScore;
            uScore = lScore = ulScore = -1;
            c--;
            if (c == -1) { // If there are no cells to the left (only boundary cells)
                ///(*alignment)[(*alignmentLength)++] = moveCode; // Move left
                ///int numUp = b * WORD_SIZE + blockPos;
                ///for (int i = 0; i < numUp; i++) // Move up until end
                ///    (*alignment)[(*alignmentLength)++] = EDLIB_EDOP_INSERT;
                ///break;
				
				alignment.push_back(moveCode);
				int numUp = b * WORD_SIZE + blockPos;
				for (int i = 0; i < numUp; i++) alignment.push_back(EDLIB_EDOP_INSERT);
				break;
            }
            if (blockPos == 0) { // If entering upper left block
                if (b == 0) { // If there are no more cells above (only boundary cells)
                    ///(*alignment)[(*alignmentLength)++] = moveCode; // Move up left
                    ///for (int i = 0; i < c + 1; i++) // Move left until end
                    ///    (*alignment)[(*alignmentLength)++] = EDLIB_EDOP_DELETE;
                    ///break;
					
					alignment.push_back(moveCode);
					for (int i = 0; i < c + 1; i++) alignment.push_back(EDLIB_EDOP_DELETE);
					break;
                }
                blockPos = WORD_SIZE - 1;
                b--;
                currP = align_matrix->Ps[c * MaxNumBlocks + b];
                currM = align_matrix->Ms[c * MaxNumBlocks + b];
            } else { // If entering left block
                blockPos--;
                currP = lP;
                currM = lM;
                currP <<= 1;
                currM <<= 1;
            }
            // Set new left block
            if (c > 0 && b >= align_matrix->first_blocks[c-1] && b <= align_matrix->last_blocks[c-1]) {
                thereIsLeftBlock = TRUE;
                lP = align_matrix->Ps[(c - 1) * MaxNumBlocks + b];
                lM = align_matrix->Ms[(c - 1) * MaxNumBlocks + b];
            } else {
                if (c == 0) { // If there are no cells to the left (only boundary cells)
                    thereIsLeftBlock = TRUE;
                    lScore = b * WORD_SIZE + blockPos + 1;
                    ulScore = lScore - 1;
                } else {
                    thereIsLeftBlock = FALSE;
                }
            }
            // Mark move
            ///(*alignment)[(*alignmentLength)++] = moveCode;
			alignment.push_back(moveCode);
        } else {
            // Reached end - finished!
            break;
        }
        //----------------------------------//
    }

    reverse(alignment.begin(), alignment.end());
    return EDLIB_STATUS_OK;
}

static void
edlibAlignmentToCigar(const unsigned char* alignment, 
					  int alignmentLength, 
					  SmallEdlibCigarFormat cigarFormat,
					  vector<EdlibGapAlignOp>& cigar)
{
	cigar.clear();
	if (cigarFormat != EEDLIB_CIGAR_EXTENDED && cigarFormat != EEDLIB_CIGAR_STANDARD) {
        return;
    }

    // Maps move code from alignment to char in cigar.
    //                        0    1    2    3
    char moveCodeToChar[] = {'=', 'I', 'D', 'X'};
    if (cigarFormat == EEDLIB_CIGAR_STANDARD) {
        moveCodeToChar[0] = moveCodeToChar[3] = 'M';
    }

    char lastMove = 0;  // Char of last move. 0 if there was no previous move.
    int numOfSameMoves = 0;
	EdlibGapAlignOp op;
    for (int i = 0; i <= alignmentLength; i++) {
        // if new sequence of same moves started
        if (i == alignmentLength || (moveCodeToChar[alignment[i]] != lastMove && lastMove != 0)) {
			op.num = numOfSameMoves;
			op.op = lastMove;
			cigar.push_back(op);
            // If not at the end, start new sequence of moves.
            if (i < alignmentLength) {
                // Check if alignment has valid values.
                if (alignment[i] > 3) {
                    cigar.clear();
                    return;
                }
                numOfSameMoves = 0;
            }
        }
        if (i < alignmentLength) {
            lastMove = moveCodeToChar[alignment[i]];
            numOfSameMoves++;
        }
    }
}

static void
cigar2AlignedString(vector<EdlibGapAlignOp>& cigar,
	const u8* query,
	const u8* target,
    string& qaln,
    string& saln)
{
	int aidx = 0;
	int qend = 0;
	int tend = 0;
	int n_cigar = cigar.size();
	EdlibGapAlignOp* oplist = cigar.data();
	
	for (int i = 0; i < n_cigar; ++i) {
		int n = oplist[i].num;
		char op = oplist[i].op;
		switch (op) {
			case 'M':
				for (int c = 0; c < n; ++c) {
					int qc = *query++;
					int tc = *target++;
					hbn_assert(qc >= 0 && qc < 4);
					hbn_assert(qc >= 0 && qc < 4);
					qc = DECODE_RESIDUE(qc);
					tc = DECODE_RESIDUE(tc);
                    qaln += qc;
                    saln += tc;
					++aidx;
					++qend;
					++tend;
				}
				break;
			case 'I':
				for (int c = 0; c < n; ++c) {
					int qc = *query++;
					hbn_assert(qc >= 0 && qc < 4);
					qc = DECODE_RESIDUE(qc);
					int tc = GAP_CHAR;
                    qaln += qc;
                    saln += tc;
					++aidx;
					++qend;
				}
				break;
			case 'D':
				for (int c = 0; c < n; ++c) {
					int tc = *target++;
					hbn_assert(tc >= 0 && tc < 4);
					int qc = GAP_CHAR;
					tc = DECODE_RESIDUE(tc);
                    qaln += qc;
                    saln += tc;
					++aidx;
					++tend;
				}
				break;
			default:
				HBN_ERR("invalid operation: %c", op);
				break;
		}
	}
}

small_edlib_align_struct*
small_edlib_align_struct_new()
{
    return (small_edlib_align_struct*)(new SmallEdlibAlignData(MaxNumBlocks, MaxSeqSize));
}

small_edlib_align_struct*
small_edlib_align_struct_free(small_edlib_align_struct* data)
{
    SmallEdlibAlignData* x_data = (SmallEdlibAlignData*)(data);
    delete x_data;
    return NULL;
}

int
small_edlib_nw(const u8* query,
		const int query_size,
		const u8* target,
		const int target_size,
		small_edlib_align_struct* _align_data,
        std::string& qaln,
        std::string& saln)
{
    qaln.clear();
    saln.clear();
    if (query_size == 0 || target_size == 0) return 1;
    SmallEdlibAlignData* align_data = (SmallEdlibAlignData*)(_align_data);
	align_data->result.clear();
	
	build_peq(query, query_size, align_data->peq);
	
	int edit_distance = -1;
	int end_position;
    int max_dist = hbn_max(query_size, target_size);
    max_dist = max_dist * 0.5;
	calc_edit_distance_nw(query,
						  query_size,
						  target,
						  target_size,
						  align_data->peq,
						  max_dist,
						  TRUE,
						  align_data->align_matrix,
						  align_data->blocks,
						  -1,
						  &edit_distance,
						  &end_position);
    if (edit_distance == -1) {
        max_dist = hbn_max(query_size, target_size);
	    calc_edit_distance_nw(query,
			query_size,
			target,
			target_size,
			align_data->peq,
			max_dist,
			TRUE,
			align_data->align_matrix,
			align_data->blocks,
			-1,
			&edit_distance,
			&end_position);        
    }
    hbn_assert(edit_distance >= 0,
        "edist_distance = %d, end_position = %d, qsize = %d, tsize = %d",
        edit_distance, end_position, query_size, target_size);
	hbn_assert(end_position == target_size - 1, 
        "edist_distance = %d, end_position = %d, qsize = %d, tsize = %d",
        edit_distance, end_position, query_size, target_size);
	
	obtainAlignmentTraceback(query_size,
							 target_size,
							 edit_distance,
							 align_data->align_matrix,
							 align_data->result.alignment);
	edlibAlignmentToCigar(align_data->result.alignment.data(),
						  align_data->result.alignment.size(),
						  EEDLIB_CIGAR_STANDARD,
						  align_data->cigar);
	cigar2AlignedString(align_data->cigar, query, target, qaln, saln);
	
	return 1;
}

int
small_edlib_nw_dist(const u8* query,
		const int query_size,
		const u8* target,
		const int target_size,
		small_edlib_align_struct* _align_data)
{
    SmallEdlibAlignData* align_data = (SmallEdlibAlignData*)(_align_data);
	align_data->result.clear();
	
	build_peq(query, query_size, align_data->peq);
	
	int edit_distance = -1;
	int end_position;
    int max_dist = hbn_max(query_size, target_size);
    max_dist = max_dist * 0.5;
	calc_edit_distance_nw(query,
						  query_size,
						  target,
						  target_size,
						  align_data->peq,
						  max_dist,
						  TRUE,
						  align_data->align_matrix,
						  align_data->blocks,
						  -1,
						  &edit_distance,
						  &end_position);
    if (edit_distance == -1) {
        max_dist = hbn_max(query_size, target_size);
	    calc_edit_distance_nw(query,
			query_size,
			target,
			target_size,
			align_data->peq,
			max_dist,
			TRUE,
			align_data->align_matrix,
			align_data->blocks,
			-1,
			&edit_distance,
			&end_position);        
    }
    hbn_assert(edit_distance >= 0);
	hbn_assert(end_position == target_size - 1, 
        "edist_distance = %d, end_position = %d, qsize = %d, tsize = %d",
        edit_distance, end_position, query_size, target_size);
	
    return edit_distance;
}
