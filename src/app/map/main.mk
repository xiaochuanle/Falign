ifeq "$(strip ${BUILD_DIR})" ""
  BUILD_DIR    := ../$(OSTYPE)-$(MACHINETYPE)/obj
endif
ifeq "$(strip ${TARGET_DIR})" ""
  TARGET_DIR   := ../$(OSTYPE)-$(MACHINETYPE)/bin
endif

TARGET   := falign
SOURCES  := \
            align_one_read.cpp \
            build_repeat_lookup_table.cpp \
            build_repeat_reference_regions.cpp \
            chain_align_list.cpp \
            extend_hit_list.cpp \
            hbn_options.cpp \
            hbn_outputs.cpp \
            hbn_word_finder.cpp \
            infer_enzyme.cpp \
            main.cpp \
            make_candidate_kmer_chain.cpp \
            map_one_volume.cpp \
            necat_info.cpp \
            smooth_pca_list.cpp \
            trim_overlap_subseq.cpp \
            ../dip3d/split_string_by_char.cpp \
            window_masker.cpp

SRC_INCDIRS  := .

TGT_LDFLAGS := -L${TARGET_DIR}
TGT_LDLIBS  := -lhbn
TGT_PREREQS := libhbn.a

SUBMAKEFILES :=
