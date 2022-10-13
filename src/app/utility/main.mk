ifeq "$(strip ${BUILD_DIR})" ""
  BUILD_DIR    := ../$(OSTYPE)-$(MACHINETYPE)/obj
endif
ifeq "$(strip ${TARGET_DIR})" ""
  TARGET_DIR   := ../$(OSTYPE)-$(MACHINETYPE)/bin
endif

TARGET   := u4falign
SOURCES  := frag_id.cpp \
            frag_sam_main.cpp \
            main.cpp \
            paf_reader.cpp \
            paf2pwc_main.cpp \
            sam_reader.cpp \
            sam2pwc_main.cpp

SRC_INCDIRS  := .

TGT_LDFLAGS := -L${TARGET_DIR}
TGT_LDLIBS  := -lhbn
TGT_PREREQS := libhbn.a

SUBMAKEFILES :=