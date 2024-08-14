ifeq "$(strip ${BUILD_DIR})" ""
  BUILD_DIR    := ../$(OSTYPE)-$(MACHINETYPE)/obj
endif
ifeq "$(strip ${TARGET_DIR})" ""
  TARGET_DIR   := ../$(OSTYPE)-$(MACHINETYPE)/bin
endif

TARGET   := dip3d
SOURCES  := \
	bam-identity-mapQ-stats.cpp \
	cat-fastq.cpp \
	compress-vcf.cpp \
	eval-hap-tag.cpp \
	extract-ashic-chr-size.cpp \
	extract-chr-bam.cpp \
	extract-chr.cpp \
	extract-mvp-het-snp.cpp \
	frag-hap.cpp \
	frag-to-ashic-read-pair.cpp \
	haplo-tag-stats.cpp \
	haplotype-consensus.cpp \
	index-bam.cpp \
	index-fasta.cpp \
	index-vcf.cpp \
	merge-snp-bam.cpp \
	main.cpp \
	make-pore-c-frag-pairs.cpp \
	paf-frag-stats.cpp \
	sam_map_info.cpp \
	sample-vcf.cpp \
	select-chr-snp-bam.cpp \
	select-snp-bam.cpp \
	sort-bam.cpp \
	split_string_by_char.cpp \
	split-bam.cpp \
	split-vcf.cpp \
	tag-bam.cpp \
	vcf-reader.cpp \
	x-tag-bam-haplotype-imputation.cpp \
	x-tag-bam-with-snp.cpp \
	x-tag-bam.cpp

SRC_INCDIRS  := .

TGT_LDFLAGS := -L${TARGET_DIR}
TGT_LDLIBS  := -lhbn
TGT_PREREQS := libhbn.a