#!/bin/bash
set -e

echo " **** Downloading source packages **** "
## Anaconda2
wget https://repo.anaconda.com/archive/Anaconda2-2018.12-Linux-x86_64.sh
echo ""
if [ -d "DeepSimulator" ];then rm -rf DeepSimulator;fi
git clone https://github.com/lykaust15/DeepSimulator.git
echo ""
## albacore
cd DeepSimulator/base_caller/albacore_2.3.1/
wget https://mirror.oxfordnanoportal.com/software/analysis/ont_albacore-2.3.1-cp36-cp36m-manylinux1_x86_64.whl
cd ../../../
echo ""
cd DeepSimulator/base_caller/guppy_3.1.5/
## ont-guppy-cpu
if [ ! -d "ont-guppy-cpu" ]
then
        wget https://mirror.oxfordnanoportal.com/software/analysis/ont-guppy-cpu_3.1.5_linux64.tar.gz
        tar xzf ont-guppy-cpu_3.1.5_linux64.tar.gz
        rm -f ont-guppy-cpu_3.1.5_linux64.tar.gz
fi
echo ""
## ont-guppy
if [ ! -d "ont-guppy" ]
then
        wget https://mirror.oxfordnanoportal.com/software/analysis/ont-guppy_3.1.5_linux64.tar.gz
        tar xzf ont-guppy_3.1.5_linux64.tar.gz
        rm -f ont-guppy_3.1.5_linux64.tar.gz
fi
cd ../../../
echo ""

echo " **** Checking Conda installation **** "
if ! [ -x "$(command -v conda)" ]
then
	echo "No Conda is installed."
	echo "Installing Anaconda2 with Python2..."
	bash Anaconda2-2018.12-Linux-x86_64.sh -b
	source anaconda2/etc/profile.d/conda.sh
	echo "Conda is successfully installed."
else
	echo "Conda has been installed."
fi
echo ""
echo " **** Checking Conda version      **** "
CONDA_PATH=`conda info --envs |grep -w 'base'|awk '{print $(NF)}'`
CONDA_PY=`${CONDA_PATH}/bin/python --version`
CONDA_EDI=${CONDA_PATH##*/}
CONDA_VER=`conda --version`
echo "Edition:" "${CONDA_EDI^}"
echo "Version:" "${CONDA_VER}"
echo ""

echo " **** Installing DeepSimulator    **** "
source ${CONDA_PATH}/etc/profile.d/conda.sh
TF_PATH=${CONDA_PATH}/envs/tensorflow_cdpm
BC_PATH=${CONDA_PATH}/envs/basecall
if [ -d "${TF_PATH}" ]
then
	echo "ERROR: Conda environment 'tensorflow_cdpm' has been created."
	echo "       Please check and rename existing environment."
	exit 1
else
	bash Anaconda2-2018.12-Linux-x86_64.sh -b -p ${TF_PATH}
	conda activate tensorflow_cdpm
	pip install protobuf==3.2.0
	pip install tensorflow==1.2.1
	pip install biopython==1.74
	pip install tflearn==0.3.2
	pip install tqdm==4.19.4
	pip install scipy==0.18.1
	pip install h5py==2.7.1
	pip install numpy==1.13.1
	pip install scikit-learn==0.20.3
fi
if [[ $CONDA_EDI == *conda2 ]];then source deactivate;else conda deactivate;fi
echo "Conda environment 'tensorflow_cdpm' is successfully configured."
echo ""
if [ -d "${BC_PATH}" ]
then
        echo "ERROR: Conda environment 'basecall' has been created."
        echo "       Please check and rename existing environment."
        exit 1
else
	conda create --name basecall python=3.6 -y
	conda activate basecall
	pip install DeepSimulator/base_caller/albacore_2.3.1/ont_albacore-2.3.1-cp36-cp36m-manylinux1_x86_64.whl
	conda install -c bioconda seqkit -y
	conda install -c bioconda seqtk -y
	conda install -c conda-forge parallel -y
	pip install pandas
	pip install numpy
	pip install joblib
fi
if [[ $CONDA_EDI == *conda2 ]];then source deactivate;else conda deactivate;fi
echo "Conda environment 'basecall' is successfully configured."
echo ""
echo "DeepSimulator is successfully installed."
echo ""

echo " **** Customing scripts for Pore-C **** "
if [[ $CONDA_EDI == *conda2 ]];
then
	cp deep_simulator.PoreC.sh DeepSimulator/deep_simulator.PoreC.sh
	cp Csim_matchup.py DeepSimulator/util/Csim_matchup.py
else
	## keep consistent with base environment:
	## use conda activate (conda3) instead of source active (conda2)
	sed 's/source/conda/g' deep_simulator.PoreC.sh > DeepSimulator/deep_simulator.PoreC.sh
	cp Csim_matchup.py DeepSimulator/util/Csim_matchup.py
	ENV_PATH=${CONDA_PATH}/envs/tensorflow_cdpm
	rm -f ${ENV_PATH}/conda-meta/conda-*
	rm -f ${ENV_PATH}/bin/conda* && ln -s ${CONDA_PATH}/bin/conda ${ENV_PATH}/bin/conda
	rm -f ${ENV_PATH}/bin/activate && ln -s ${CONDA_PATH}/bin/activate ${ENV_PATH}/bin/activate
	rm -f ${ENV_PATH}/bin/deactivate && ln -s ${CONDA_PATH}/bin/deactivate ${ENV_PATH}/bin/deactivate
	ENV_PATH=${CONDA_PATH}/envs/basecall
        rm -f ${ENV_PATH}/conda-meta/conda-*
        rm -f ${ENV_PATH}/bin/conda* && ln -s ${CONDA_PATH}/bin/conda ${ENV_PATH}/bin/conda
        rm -f ${ENV_PATH}/bin/activate && ln -s ${CONDA_PATH}/bin/activate ${ENV_PATH}/bin/activate
        rm -f ${ENV_PATH}/bin/deactivate && ln -s ${CONDA_PATH}/bin/deactivate ${ENV_PATH}/bin/deactivate
fi
echo "Custom scripts:"
echo "${PATH}/DeepSimulator/deep_simulator.PoreC.sh"
echo "${PATH}/DeepSimulator/util/Csim_matchup.py"
echo ""
echo "Installation finished."
