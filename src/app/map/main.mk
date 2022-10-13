ifeq "$(strip ${BUILD_DIR})" ""
  BUILD_DIR    := ../$(OSTYPE)-$(MACHINETYPE)/obj
endif
ifeq "$(strip ${TARGET_DIR})" ""
  TARGET_DIR   := ../$(OSTYPE)-$(MACHINETYPE)/bin
endif

TARGET   := falign
SOURCES  := align_one_read.cpp \
			chain_align_list.cpp \
			extend_hit_list.cpp \
			hbn_options.cpp \
			main.cpp \
			map_one_volume.cpp \
			trim_overlap_subseq.cpp \
			smooth_pca_list.cpp

SRC_INCDIRS  := .

TGT_LDFLAGS := -L${TARGET_DIR}
TGT_LDLIBS  := -lhbn
TGT_PREREQS := libhbn.a

SUBMAKEFILES :=
