ifeq "$(strip ${BUILD_DIR})" ""
  BUILD_DIR    := ../$(OSTYPE)-$(MACHINETYPE)/obj
endif
ifeq "$(strip ${TARGET_DIR})" ""
  TARGET_DIR   := ../$(OSTYPE)-$(MACHINETYPE)/bin
endif

TARGET       := libhbn.a

SOURCES      := \
	./corelib/cstr_util.c \
	./corelib/db_format.cpp \
	./corelib/fasta.cpp \
	./corelib/frag_id.cpp \
	./corelib/getMemorySize.c \
	./corelib/getRSS.c \
	./corelib/hbn_aux.c \
	./corelib/hbn_package_version.c \
	./corelib/lookup_table.cpp \
	./corelib/read_list.cpp \
	./corelib/restrict_enzyme_loci_list.cpp \
	./corelib/small_object_alloc.c \
	./ncbi_blast/c_ncbi_blast_aux.c \
	./ncbi_blast/ncbi_blast_aux.cpp \
	./ncbi_blast/str_util/ncbistr_util.cpp \
	./ncbi_blast/str_util/ncbistr.cpp \
	./ncbi_blast/str_util/str_cmp.cpp \
	./ncbi_blast/str_util/numeric_str_interconv.cpp \
	./ncbi_blast/str_util/str_util.cpp \
	./sw/align.c \
	./sw/dalign.cpp \
	./sw/edlib_wrapper.cpp \
	./sw/edlib.cpp \
	./sw/hbn_traceback_aux.c \
	./sw/hbn_traceback.cpp \
	./sw/ksw2_extd2_sse.c \
	./sw/ksw2_extz2_sse.c \
	./sw/ksw2_wrapper.cpp \
	./sw/small_edlib_align.cpp \

SRC_INCDIRS  := 

SUBMAKEFILES := ./app/map/main.mk ./app/dip3d/main.mk
