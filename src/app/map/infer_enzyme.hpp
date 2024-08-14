#ifndef __INFER_ENZYME_HPP
#define __INFER_ENZYME_HPP

#include "../../corelib/lookup_table.hpp"
#include "../../corelib/unpacked_seqdb.hpp"
#include "hbn_options.hpp"
#include "hbn_word_finder.hpp"

#include <string>

std::string
infer_enzyme_mt(const HbnProgramOptions* options,
    HbnUnpackedDatabase* subjects,
    HbnUnpackedDatabase* queries,
    HbnWordFinder* word_finder);

#endif // __INFER_ENZYME_HPP