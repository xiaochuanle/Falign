ifeq "$(strip ${BUILD_DIR})" ""
  BUILD_DIR    := ../$(OSTYPE)-$(MACHINETYPE)/obj
endif
ifeq "$(strip ${TARGET_DIR})" ""
  TARGET_DIR   := ../$(OSTYPE)-$(MACHINETYPE)/bin
endif

TARGET       := libhbn.a

SOURCES      := \
	./corelib/cstr_util.c \
	./corelib/db_format.c \
	./corelib/fasta.c \
	./corelib/gapped_candidate.c \
	./corelib/hbn_aux.c \
	./corelib/hbn_package_version.c \
	./corelib/kstring.c \
	./corelib/line_reader.c \
	./corelib/restrict_enzyme_loci_list.c \
	./corelib/small_object_alloc.c \
	./algo/align_format/format_aux.c \
	./algo/align_format/paf.c \
	./algo/align_format/sam.c \
	./algo/align.c \
	./algo/dalign.c \
	./algo/diff_gapalign.cpp \
	./algo/edlib.cpp \
	./algo/edlib_wrapper.c \
	./algo/hash_list_bucket_sort.c \
	./algo/hbn_lookup_table.c \
	./algo/hbn_traceback.c \
	./algo/hbn_traceback_aux.c \
	./algo/hbn_word_finder.c \
	./algo/kalloc.c \
	./algo/ksw2_extd2_sse.c \
	./algo/ksw2_extz2_sse.c \
	./algo/ksw2_wrapper.c \
	./algo/make_candidate_kmer_chain.cpp \
	./algo/seq_loader.c \
	./algo/small_edlib_align.c \
	./ncbi_blast/c_ncbi_blast_aux.c \
	./ncbi_blast/ncbi_blast_aux.cpp \
	./ncbi_blast/str_util/ncbistr_util.cpp \
	./ncbi_blast/str_util/ncbistr.cpp \
	./ncbi_blast/str_util/str_cmp.cpp \
	./ncbi_blast/str_util/numeric_str_interconv.cpp \
	./ncbi_blast/str_util/str_util.cpp \

SRC_INCDIRS  := 

SUBMAKEFILES := ./app/map/main.mk ./app/map_ngf/main.mk ./app/utility/main.mk
