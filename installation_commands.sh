# = Path definition

# path to installation folder
export INSTALL_PATH=/user/student.aau.dk/sheide17
export INSTALL_PATH=/shared-nfs/SH/software

# =================================================================================================
# = GO
wget https://golang.org/dl/go1.17.linux-amd64.tar.gz
rm -rf ${INSTALL_PATH}/go && \
tar -C ${INSTALL_PATH} -xzf go1.17.linux-amd64.tar.gz
rm go1.17.linux-amd64.tar.gz

# Add go to PATH
echo 'export PATH=${INSTALL_PATH}/go/bin:$PATH' >> ~/.bashrc && \
source ~/.bashrc

# =================================================================================================
# = singularity 3.8.2
export VERSION=3.8.2 && \
wget https://github.com/sylabs/singularity/releases/download/v${VERSION}/singularity-ce-${VERSION}.tar.gz && \
tar -xzf singularity-ce-${VERSION}.tar.gz && \
cd singularity-ce-${VERSION}
mkdir ${INSTALL_PATH}/singularity
./mconfig --without-suid --prefix=${INSTALL_PATH}/singularity && \
    make -C ./builddir && \
    make -C ./builddir install
rm singularity-ce-${VERSION}.tar.gz

# Add singularity to PATH
echo 'export PATH=${INSTALL_PATH}/singularity/bin:$PATH' >> ~/.bashrc && \
source ~/.bashrc

# =================================================================================================
# = nanodisco
cd $INSTALL_PATH
mkdir nanodisco
cd nanodisco

singularity pull --name nanodisco.sif library://fanglab/default/nanodisco
singularity build nd_env nanodisco.sif

# =================================================================================================
# = conda
cd ${INSTALL_PATH}
# miniconda
wget https://repo.anaconda.com/miniconda/Miniconda3-py39_4.10.3-Linux-x86_64.sh
# anaconda3
wget https://repo.anaconda.com/archive/Anaconda3-2021.05-Linux-x86_64.sh
# Change which script to run to install either miniconda or anaconda3
bash Miniconda3-py39_4.10.3-Linux-x86_64.sh

# Add conda to PATH
echo 'export PATH=${INSTALL_PATH}/anaconda3/bin:$PATH' >> ~/.bashrc && \
source ~/.bashrc

# create enviroment
conda create --name py38 python=3.8

# =================================================================================================
# pycharm
mkdir $INSTALL_PATH/pycharm
# Install zipped community version into folder
tar -xzf "/shared-nfs/SH/software/pycharm/pycharm-community-anaconda-2019.3.5.tar.gz"

echo 'export PATH=/shared-nfs/SH/software/pycharm/pycharm-community-anaconda-2019.3.5/bin:$PATH' >> ~/.bashrc && \
source ~/.bashrc

# =================================================================================================
# = Poretools
cd $INSTALL_PATH
git clone https://github.com/arq5x/poretools
cd poretools
python setup.py install --user
export PATH=$PATH:/home/arq5x/.local/bin

# =================================================================================================
# = Java
cd /shared-nfs/SH/software
wget "https://javadl.oracle.com/webapps/download/AutoDL?BundleId=245050_d3c52aa6bfa54d3ca74e617f18309292"

tar zxvf jre-8u301-linux-i586.tar.gz

# =================================================================================================
# = File transfer IPs
# GPU server @172.19.8.14
# 1024b      @172.19.20.204

scp \
  "sheide17@student.aau.dk@IP:/path/on/server/one/file" \
  "sheide17@student.aau.dk@IP:/path/on/server/two/file"

scp -r \
  "sheide17@student.aau.dk@IP:/path/on/server/one/directionary" \
  "sheide17@student.aau.dk@IP:/path/on/server/two/directionary"

scp -r \
  "sheide17@student.aau.dk@172.19.8.14:/user/student.aau.dk/sheide17/projects/current_difference/zymo/megalodon_unmod_v2/signal_mappings.hdf5" \
  "/shared-nfs/SH/samples/zymo/megalodon_pcr_test"


# =================================================================================================
# = Add executables directionary to PATH
# Add conda to PATH
echo 'export PATH=/shared-nfs/SH/code/scripts/executables:$PATH' >> ~/.bashrc && \
source ~/.bashrc
