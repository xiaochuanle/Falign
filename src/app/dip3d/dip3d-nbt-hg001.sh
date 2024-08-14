#!/bin/bash

### 1) Output directory of Dip3D
I_output=nbt-hg001-dip3d-output

### 2) Absolute path of reference genome
I_reference=/data1/chenying/dip3d-1/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.major.fna

### 3) Chromosome list
I_CHR_LIST=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY)

### 4) Number of CPU threads used by Dip3D
I_threads=48

### 5) Directory storing Clair Pore-C model
###   The directory should contains four data file
###   1. full_alignment.data-00000-of-00001
###   2. full_alignment.index
###   3. pileup.data-00000-of-00001
###   4. pileup.index
I_clair3_model_directory=/data1/chenying/dip3d-1/dip3d/clair3_pore_c_model

### 6) Options for calling SNP by Clair3
###   1. Coverage of Pore-C data
I_snp_frag_coverage=90
###   2. Minimum Pore-C fragment length
I_snp_frag_length=100
###   3. Minimum alignment identities
I_snp_frag_identity=90.0
###   4. Minimum mapping qualities
I_snp_frag_mapQ=5

### 7) Options for haplo-tagging Pore-C fragments with SNP
###   1. Number of flanking match residues in overlapped SNP sites
I_snp_tag_match_bases=3
###   2. Only haplo-tagging Pore-C fragments with at least this mapping qualities
I_snp_tag_mapQ=5
###   3. Only haplo-tagging Pore-C fragments with at least this aligment identites
I_snp_tag_identity=85.0

### 8) Options for haplotype imputation
### 1. Maximum genomic distance between adjacent fragment pairs
I_imputation_adj_dist=29500000
### 2. Maximum genomic distacne between non-adjacent fragment pairs
I_imputation_non_adj_dist=16500000

### 9) Options for running ASHIC algorithm
###   1. Absolute path of ASHIC
I_ASHIC_PATH=/data1/chenying/dip3d-1/dip3d/third-party/ASHIC
###   2. Runing ASHIC algorithm (1) or not (0)
I_RUN_ASHIC=0
###   3. Resulotions of contact matrices
I_ASHIC_RES=25000
###   4. Minimum mapping qualities of Pore-C fragments
I_ASHIC_mapQ=5
###   5. Which algorithm to perform by ASHIC
I_ASHIC_MODEL="ASHIC-ZIPM"

### 10) Absolute path list of Pore-C list
###     Please don't concat concat FASTQ files genergeted with different enzymes into one
I_pore_c_reads=( \
/data3/nbt-pore-c-data/fq-hg001/SRR11589389.fastq.gz
/data3/nbt-pore-c-data/fq-hg001/SRR11589390.fastq.gz
/data3/nbt-pore-c-data/fq-hg001/SRR11589391.fastq.gz
/data3/nbt-pore-c-data/fq-hg001/SRR11589392.fastq.gz
/data3/nbt-pore-c-data/fq-hg001/SRR11589393.fastq.gz
/data3/nbt-pore-c-data/fq-hg001/SRR11589394.fastq.gz
/data3/nbt-pore-c-data/fq-hg001/SRR11589395.fastq.gz
/data3/nbt-pore-c-data/fq-hg001/SRR11589396.fastq.gz
/data3/nbt-pore-c-data/fq-hg001/SRR11589397.fastq.gz
/data3/nbt-pore-c-data/fq-hg001/SRR11589398.fastq.gz
/data3/nbt-pore-c-data/fq-hg001/SRR11589399.fastq.gz
/data3/nbt-pore-c-data/fq-hg001/SRR11589400.fastq.gz
/data3/nbt-pore-c-data/fq-hg001/SRR11589401.fastq.gz
/data3/nbt-pore-c-data/fq-hg001/SRR11589402.fastq.gz
/data3/nbt-pore-c-data/fq-hg001/SRR11589403.fastq.gz
)

########################

G_HOME=$(pwd)/${I_output}

G_REFERENCE_DIR=${G_HOME}/0-reference

G_FALIGN_DIR=${G_HOME}/1-falign

G_CHR_BAM_DIR=${G_HOME}/2-chr-bams

G_SNP_DIR=${G_HOME}/3-snp

G_BAM_HAPLOTAG_DIR=${G_HOME}/4-bam-haplotag

G_ASHIC_DIR=${G_HOME}/5-ashic-1

########################

function validate_params()
{
	if [ ! -n "${I_output}" ]; then
		echo "[$(date)] ERROR: output directory (I_output) is not provided"
		exit 1
	else
		echo "Output directory:"
		echo "    ${I_output}"
	fi

	if [ ! -n "${I_reference}" ]; then
		echo "[$(date)] ERROR: reference genome (I_reference) is not provided"
		exit 1
	else
		echo "Reference genome:"
		echo "    ${I_reference}"
	fi

	if [ ${#I_CHR_LIST[@]} -eq 0 ]; then
		echo "[$(date)] ERROR: chromosome list (I_CHR_LIST) is empty"
		exit 1
	else
		echo "Chromosomes:"
		echo "    ${I_CHR_LIST[@]}"
	fi

	if [ ! -n "${I_threads}" ]; then
		echo "[$(date)] ERROR: cpu threads (I_threads) is not provided"
		exit 1
	else
		echo "CPU threads:"
		echo "    ${I_threads}"
	fi
	
	if [ ${#I_pore_c_reads[@]} -eq 0 ]; then
		echo "[$(date)] ERROR: Pore-C reads (I_pore_c_reads) are not provided"
		exit 1
	else
		echo "Pore-C reads:"
		for R in ${I_pore_c_reads[@]}
		do
			echo "    ${R}"
		done
	fi
}

function prepare_reference()
{
	local all_is_done=${G_REFERENCE_DIR}/all.done
	if [ -f ${all_is_done} ]; then
		return
	fi

	mkdir -p ${G_REFERENCE_DIR}
	if [ $? -ne 0 ]; then
		echo "[$(date)] ERROR: fail at creating directory ${G_REFERENCE_DIR}"
		eixt 1
	fi

	if [ ! -f ${I_reference} ]; then
		echo "[$(date)] ERROR: Reference genome '${I_reference}' does not exist"
		exit 1
	fi

	local reference=${G_REFERENCE_DIR}/reference.fa
	if [ ! -f ${reference} ]; then
		cp ${I_reference} ${reference}
		if [ $? -ne 0 ]; then
			exit 1
		fi
	fi
	
	if [ ! -f ${reference}.fai ]; then
		CMD="dip3d index-fasta ${reference}"
		echo "[$(date)] Running command"
		echo "  ${CMD}"
		${CMD}
		if [ $? -ne 0 ]; then
			echo "[$(date)] ERROR fail at running command"
			echo "  ${CMD}"
			exit 1
		fi
	fi

	if [ ! -f ${reference}.wm.bed ]; then 
		CMD="falign-wm ${reference} ${reference}.wm.bed"
		echo "[$(date)] Running command"
		echo "  ${CMD}"
		${CMD}
		if [ $? -ne 0 ]; then
			echo "[$(date)] ERROR fail at running command"
			echo "  ${CMD}"
			exit 1
		fi
	fi

	touch ${all_is_done}
}

function map_pore_c_reads_with_falign()
{
	local all_is_done=${G_FALIGN_DIR}/all.done
	if [ -f ${all_is_done} ]; then
		return
	fi

	mkdir -p ${G_FALIGN_DIR}
	if [ $? -ne 0 ]; then
		exit 1
	fi

	local bam=${G_FALIGN_DIR}/frag.map.bam
	local opts="-repeat_bed ${G_REFERENCE_DIR}/reference.fa.wm.bed -num_threads ${I_threads} -outfmt frag-bam -out ${bam}"
	if [ -n "${I_enzyme}" ]; then
		opts="${opts} -enzyme ${I_enzyme}"
	fi

	local CMD="falign ${opts} ${G_REFERENCE_DIR}/reference.fa ${I_pore_c_reads[@]}"
	echo "[$(date)] Running command"
	echo "  ${CMD}"
	${CMD}
	if [ $? -ne 0 ]; then
		echo "[$(date)] fail"
		exit 1
	fi

	touch ${all_is_done}
}

function split_chr_bams()
{
	local all_is_done=${G_CHR_BAM_DIR}/all.done
	
	if [ -f ${all_is_done} ]; then
		return
	fi

	mkdir -p ${G_CHR_BAM_DIR}
	if [ $? -ne 0 ]; then
		exit 1
	fi

	local bam=${G_FALIGN_DIR}/frag.map.bam
	CMD="dip3d split-bam ${bam} ${G_CHR_BAM_DIR} ${I_CHR_LIST[@]}"
	echo "[$(date)] Running command"
	echo "  ${CMD}"
	${CMD}
	if [ $? -ne 0 ]; then
		exit 1
	fi

	touch ${all_is_done}
}

function s_select_snp_bam()
{
	local j_select_snp_bam=${G_SNP_DIR}/select-snp-bam.done
	if [ -f ${j_select_snp_bam} ]; then
		return
	fi

	local snp_bam=${G_SNP_DIR}/sorted.snp.bam
	local snp_bam_opts="-l ${I_snp_frag_length} -q ${I_snp_frag_mapQ} -i ${I_snp_frag_identity}"
	local CMD="dip3d select-chr-snp-bam ${snp_bam_opts} ${G_CHR_BAM_DIR} ${I_snp_frag_coverage} ${snp_bam} ${I_CHR_LIST[@]}"
	echo "[$(date)] Running command"
	echo "  ${CMD}"
	${CMD}
	if [ $? -ne 0 ]; then
		exit 1
	fi

	touch ${j_select_snp_bam}
}

function s_call_snp_clair3_v0.1_r10()
{
	local j_call_snp=${G_SNP_DIR}/call-snp.done
	if [ -f ${j_call_snp} ]; then
		return
	fi

	local clair3_output_dir=${G_SNP_DIR}/clair3-snp
	local ref_dir=${G_REFERENCE_DIR}
	local bam_dir=${G_SNP_DIR}
	local model_dir=${I_clair3_model_directory}

	rm -rf ${clair3_output_dir}
	mkdir -p ${clair3_output_dir}
	if [ $? -ne 0 ]; then 
		exit 1
	fi

	docker run -it --user $(id -u):$(id -g) \
		-v ${ref_dir}:${ref_dir} \
		-v ${clair3_output_dir}:${clair3_output_dir} \
		-v ${bam_dir}:${bam_dir} \
		-v ${model_dir}:${model_dir} \
		hkubal/clair3:v0.1-r10 \
		/opt/bin/run_clair3.sh \
		--bam_fn=${bam_dir}/sorted.snp.bam \
		--ref_fn=${ref_dir}/reference.fa \
		--threads=${I_threads} \
		--platform="ont" \
		--model_path=${model_dir} \
		--output=${clair3_output_dir} \
		--var_pct_full=0.1 \
		--call_snp_only \
		--snp_min_af=0.25

	if [ $? -ne 0 ]; then
		exit 1
	fi

	if [ ! -f ${G_SNP_DIR}/clair3-snp/merge_output.vcf.gz ]; then
		exit 1
	fi

	if [ ! -f ${G_SNP_DIR}/clair3-snp/merge_output.vcf.gz.tbi ]; then
		exit 1
	fi

	touch ${j_call_snp}
}

function s_call_snp_clair3_v1.0.5()
{
	local j_call_snp=${G_SNP_DIR}/call-snp.done
	if [ -f ${j_call_snp} ]; then
		return
	fi

	local clair3_output_dir=${G_SNP_DIR}/clair3-snp
	local ref_dir=${G_REFERENCE_DIR}
	local bam_dir=${G_SNP_DIR}
	local model_dir=${I_clair3_model_directory}

	rm -rf ${clair3_output_dir}
	mkdir -p ${clair3_output_dir}
	if [ $? -ne 0 ]; then 
		exit 1
	fi

	docker run -it --user $(id -u):$(id -g) \
		-v ${ref_dir}:${ref_dir} \
		-v ${clair3_output_dir}:${clair3_output_dir} \
		-v ${bam_dir}:${bam_dir} \
		-v ${model_dir}:${model_dir} \
		hkubal/clair3:v1.0.5 \
		/opt/bin/run_clair3.sh \
		-b ${bam_dir}/sorted.snp.bam \
		-f ${ref_dir}/reference.fa \
		-m ${model_dir} \
		-t ${I_threads} \
		-p ont \
		--disable_c_impl \
		-o ${clair3_output_dir} \
		--threads=${I_threads} \
		--snp_min_af=0.25 \
		--call_snp_only 

	if [ $? -ne 0 ]; then
		exit 1
	fi

	if [ ! -f ${G_SNP_DIR}/clair3-snp/merge_output.vcf.gz ]; then
		exit 1
	fi

	if [ ! -f ${G_SNP_DIR}/clair3-snp/merge_output.vcf.gz.tbi ]; then
		exit 1
	fi

	touch ${j_call_snp}
}

function s_split_chr_snp()
{
	local j_split_chr_snp=${G_SNP_DIR}/split-chr-snp.done
	if [ -f ${j_split_chr_snp} ]; then
		return
	fi

	if [ ! -f ${G_SNP_DIR}/clair3-snp/merge_output.vcf.gz ]; then
		echo "[$(date)] Clair3 VCF ${G_SNP_DIR}/clair3-snp/merge_output.vcf.gz does not exist"
		exit 1
	fi

	if [ ! -f ${G_SNP_DIR}/clair3-snp/merge_output.vcf.gz.tbi ]; then
		echo "[$(date)] Clair3 VCF index ${G_SNP_DIR}/clair3-snp/merge_output.vcf.gz.tbi does not exist"
		exit 1
	fi

	for chr in ${I_CHR_LIST[@]}
	do
		local dir=${G_SNP_DIR}/${chr}
		mkdir -p ${dir}
		if [ $? -ne 0 ]; then
			echo "[$(date)] Could not create directory ${dir}"
			exit 1
		fi
	done 

	CMD="dip3d split-vcf -v -p ${G_SNP_DIR}/clair3-snp/merge_output.vcf.gz ${G_SNP_DIR} ${I_CHR_LIST[@]}"
	echo "[$(date)] Running command"
	echo "  ${CMD}"
	${CMD}
	if [ $? -ne 0 ]; then
		exit
	fi 

	touch ${j_split_chr_snp}
}

function call_snp()
{
	mkdir -p ${G_SNP_DIR}
	if [ $? -ne 0 ]; then
		exit 1
	fi 

	s_select_snp_bam

	s_call_snp_clair3_v1.0.5

	s_split_chr_snp
}

function s_make_hapcut_frag_pairs()
{
	for chr in ${I_CHR_LIST[@]}
	do
		local j_extract_hapcut_frag=${G_SNP_DIR}/${chr}/extract-hapcut-frag.done
		if [ -f ${j_extract_hapcut_frag} ]; then
			continue
		fi 
		local chr_snp_bam=${G_SNP_DIR}/clair3-snp/tmp/phase_output/phase_bam/${chr}.bam
		CMD="extractHAIRS --ont 1 --mmq 5 --ref ${G_REFERENCE_DIR}/reference.fa --VCF ${G_SNP_DIR}/${chr}/${chr}.snp.vcf --bam ${chr_snp_bam} --out ${G_SNP_DIR}/${chr}/${chr}.frag"
		echo "[$(date)] Running command"
		echo "  ${CMD}"
		${CMD}
		if [ $? -ne 0 ]; then
			exit 1
		fi
		touch ${j_extract_hapcut_frag}
	done 

	for chr in ${I_CHR_LIST[@]}
	do
		local j_make_pore_c_frag_pair=${G_SNP_DIR}/${chr}/make-pore-c-frag-pair.done
		if [ -f ${j_make_pore_c_frag_pair} ]; then
			continue
		fi
		CMD="dip3d make-pore-c-frag-pair ${G_SNP_DIR}/${chr}/${chr}.frag ${G_SNP_DIR}/${chr}/${chr}.frag-pair"
		echo "[$(date)] Running command"
		echo "  ${CMD}"
		${CMD}
		if [ $? -ne 0 ]; then
			exit 1
		fi
		touch ${j_make_pore_c_frag_pair}
	done 
}

function s_phase_snp()
{
	local j_hapcut2_phase=${G_SNP_DIR}/hapcut2-phase.done
	if [ ! -f ${j_hapcut2_phase} ]; then
		echo "[$(date)] HAPCUT2 haplotype assembly"
		parallel -j ${I_threads} \
			"HAPCUT2 --fragments ${G_SNP_DIR}/{}/{}.frag-pair \
			--VCF ${G_SNP_DIR}/{}/{}.snp.vcf \
			--output ${G_SNP_DIR}/{}/{}.hapcut \
			--hic 1 \
			> ${G_SNP_DIR}/{}/{}.hapcut2.log 2>&1" ::: ${I_CHR_LIST[@]}
		ERR=$?
		for chr in ${I_CHR_LIST[@]}
		do
			cat ${G_SNP_DIR}/${chr}/${chr}.hapcut2.log
		done 
		if [ ${ERR} -ne 0 ]; then
			echo "[$(date)] Fail at haplotype assembly with HAPCUT2"
			exit 1
		fi 
		touch ${j_hapcut2_phase}
	fi

	local j_whatshap_phase=${G_SNP_DIR}/whatshap-phase.done
	if [ ! -f ${j_whatshap_phase} ]; then
		echo "[$(date)] whatshap phase"
		parallel -j ${I_threads} \
		"whatshap phase --reference ${G_REFERENCE_DIR}/reference.fa \
		${G_SNP_DIR}/{}/{}.snp.vcf \
		${G_SNP_DIR}/{}/{}.hapcut.phased.VCF \
		${G_SNP_DIR}/clair3-snp/tmp/phase_output/phase_bam/{}.bam \
		-o ${G_SNP_DIR}/{}/{}.whatshap.vcf \
		--ignore-read-groups \
		> ${G_SNP_DIR}/{}/{}.whatshap.log 2>&1" ::: ${I_CHR_LIST[@]}
		ERR=$?
		for chr in ${I_CHR_LIST[@]}
		do
			cat ${G_SNP_DIR}/${chr}/${chr}.whatshap.log
		done 
		if [ ${ERR} -ne 0 ]; then
			echo "[$(date)] Fail at vcf phase with whatshap"
			exit 1
		fi
		touch ${j_whatshap_phase}
	fi

	local j_extract_mvp_block=${G_SNP_DIR}/extract-mvp-block.done
	if [ ! -f ${j_extract_mvp_block} ]; then
		for chr in ${I_CHR_LIST[@]}
		do
			CMD="dip3d extract-mvp-het-snp ${G_SNP_DIR}/${chr}/${chr}.whatshap.vcf ${G_SNP_DIR}/${chr}/${chr}.mvpblock.vcf"
			echo "[$(date)] Running command"
			echo "  ${CMD}"
			${CMD}
			if [ $? -ne 0 ]; then
				exit 1
			fi
		done 
		touch ${j_extract_mvp_block}
	fi
}

function phase_snp()
{
	s_make_hapcut_frag_pairs

	s_phase_snp
}

function bam_haplotag()
{
	mkdir -p ${G_BAM_HAPLOTAG_DIR}
	if [ $? -ne 0 ]; then
		exit 1
	fi

	for chr in ${I_CHR_LIST[@]}
	do
		local WDIR=${G_BAM_HAPLOTAG_DIR}/${chr}
		mkdir -p ${WDIR}
		if [ $? -ne 0 ]; then
			exit 1
		fi
		local j_tag_bam=${WDIR}/tag-bam.done
		if [ -f ${j_tag_bam} ]; then
			continue
		fi
		local opts="-q ${I_snp_tag_mapQ} -m ${I_snp_tag_match_bases} -i ${I_snp_tag_identity} -a ${I_imputation_adj_dist} -d ${I_imputation_non_adj_dist} -t ${I_threads}"
		CMD="dip3d tag-bam ${opts} ${WDIR} ${G_REFERENCE_DIR}/reference.fa ${G_SNP_DIR}/${chr}/${chr}.mvpblock.vcf ${G_CHR_BAM_DIR}/${chr}/${chr}.bam"
		echo "[$(date)] Running command"
		echo "  ${CMD}"
		${CMD}
		if [ $? -ne 0 ]; then
			exit 1
		fi
		touch ${j_tag_bam}
	done 
}

function bam_haplotag_parallel()
{
	mkdir -p ${G_BAM_HAPLOTAG_DIR}
	if [ $? -ne 0 ]; then
		exit 1
	fi

	local j_tag_bam=${G_BAM_HAPLOTAG_DIR}/tag-bam.done
	if [ -f ${j_tag_bam} ]; then
		return
	fi

	for chr in ${I_CHR_LIST[@]}
	do
		mkdir -p ${G_BAM_HAPLOTAG_DIR}/${chr}
		if [ $? -ne 0 ]; then
			exit 1
		fi
	done

	echo "[$(date)] Haplotagging BAMs"
	local opts="-q ${I_snp_tag_mapQ} -m ${I_snp_tag_match_bases} -i ${I_snp_tag_identity} -a ${I_imputation_adj_dist} -d ${I_imputation_non_adj_dist} -t ${I_threads}"
	parallel -j 8 \
		"dip3d tag-bam ${opts} \
		${G_BAM_HAPLOTAG_DIR}/{} \
		${G_REFERENCE_DIR}/reference.fa \
		${G_SNP_DIR}/{}/{}.mvpblock.vcf \
		${G_CHR_BAM_DIR}/{}/{}.bam \
		> ${G_BAM_HAPLOTAG_DIR}/{}/tag-bam.log 2>&1" ::: ${I_CHR_LIST[@]}
	ERR=$?
	for chr in ${I_CHR_LIST[@]}
	do
		cat ${G_BAM_HAPLOTAG_DIR}/${chr}/tag-bam.log
	done
	if [ ${ERR} -ne 0 ]; then
		exit 1
	fi
}

function ashic_dip_mat_construction()
{
	mkdir -p ${G_ASHIC_DIR}
	if [ $? -ne 0 ]; then
		exit 1
	fi 

	local j_all=${G_ASHIC_DIR}/all.done
	if [ -f ${j_all} ]; then
		return
	fi

	for chr in ${I_CHR_LIST[@]}
	do
		mkdir -p ${G_ASHIC_DIR}/${chr}
		if [ $? -ne 0 ]; then
			exit 1
		fi
	done

	local j_extract_ashic_chr_size=${G_ASHIC_DIR}/extract-chr-size.done
	if [ ! -f ${j_extract_ashic_chr_size} ]; then
		CMD="dip3d extract-ashic-chr-size ${G_REFERENCE_DIR}/reference.fa ${G_ASHIC_DIR}/chrom.sizes"
		echo "[$(date)] Running command"
		echo "  ${CMD}"
		${CMD}
		if [ $? -ne 0 ]; then
			exit 1
		fi
		touch ${j_extract_ashic_chr_size}
	fi

	local j_make_ashic_read_pair=${G_ASHIC_DIR}/make-ashic-read-pair.done
	if [ ! -f ${j_make_ashic_read_pair} ]; then
		for chr in ${I_CHR_LIST[@]}
		do
			CMD="dip3d frag-to-ashic-read-pair ${G_BAM_HAPLOTAG_DIR}/${chr}/imputed-frag-hap-list ${chr} ${I_ASHIC_mapQ} ${G_ASHIC_DIR}/${chr}/ashic-read-pair"
			echo "[$(date)] Running command"
			echo "  ${CMD}"
			${CMD}
			if [ $? -ne 0 ]; then
				exit 1
			fi
		done 
		touch ${j_make_ashic_read_pair}
	fi

	local j_ashic_split2chrs=${G_ASHIC_DIR}/ashic-split2chr2.done
	if [ ! -f ${j_ashic_split2chrs} ];
	then
		for chr in ${I_CHR_LIST[@]}
		do
			CMD="python ${I_dip_path}/third-party/ASHIC/ashic/cli/ashic_data.py split2chrs --chr1 1 --allele1 3 --chr2 4 --allele2 6 ${G_ASHIC_DIR}/${chr}/ashic-read-pair ${G_ASHIC_DIR}/${chr}"
			echo "[$(date)] Running command"
			echo "  ${CMD}"
			${CMD}
			if [ $? -ne 0 ]; then
				exit 1
			fi
		done 
		touch ${j_ashic_split2chrs}
	fi

	local j_ashic_binning=${G_ASHIC_DIR}/ashic-binning.done
	if [ ! -f ${j_ashic_binning} ]; then
		for chr in ${I_CHR_LIST[@]}
		do
			CMD="python ${I_dip_path}/third-party/ASHIC/ashic/cli/ashic_data.py binning --c1=1 --p1=2 --a1=3 --c2=4 --p2=5 --a2=6 --res=${I_ASHIC_RES} --chrom=${chr} --genome ${G_ASHIC_DIR}/chrom.sizes ${G_ASHIC_DIR}/${chr}/ashic-read-pair ${G_ASHIC_DIR}/${chr}"
			echo "[$(date)] Running command"
			echo "  ${CMD}"
			${CMD}
			if [ $? -ne 0 ]; then
				exit 1
			fi
		done 
		touch ${j_ashic_binning}
	fi 

	local j_ashic_pack=${G_ASHIC_DIR}/ashic-pack.done
	if [ ! -f ${j_ashic_pack} ]; then
		for chr in ${I_CHR_LIST[@]}
		do
			CMD="python ${I_dip_path}/third-party/ASHIC/ashic/cli/ashic_data.py pack ${G_ASHIC_DIR}/${chr} ${G_ASHIC_DIR}/${chr}"
			echo "[$(date)] Running command"
			echo "  ${CMD}"
			${CMD}
			if [ $? -ne 0 ]; then
				exit 1
			fi
		done
		touch ${j_ashic_pack}
	fi

	local j_run_ashic=${G_ASHIC_DIR}/run-ashic.done
	if [ ! -f ${j_run_ashic} ]; then
		for chr in ${I_CHR_LIST[@]}
		do
			CMD="python ${I_dip_path}/third-party/ASHIC/ashic/__main__.py -i ${G_ASHIC_DIR}/${chr}/ashic-read-pair_${chr}_${I_ASHIC_RES}.pickle -o ${G_ASHIC_DIR}/${chr}-${I_ASHIC_RES}-${I_ASHIC_MODEL} --model ${I_ASHIC_MODEL}"
			echo "[$(date)] Running command"
			echo "  ${CMD}"
			${CMD}
			if [ $? -ne 0 ]; then
				exit 1
			fi			
		done 
		touch ${j_run_ashic}
	fi

}

validate_params

prepare_reference

map_pore_c_reads_with_falign

split_chr_bams

call_snp

phase_snp

bam_haplotag
