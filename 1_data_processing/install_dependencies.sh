#!/bin/bash

# Exit script on error
set -e

# Update system packages
sudo apt-get update

# Install necessary tools
sudo apt-get install -y git wget

# Install Conda (assuming Miniconda for a lightweight installation)
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh
bash ~/miniconda.sh -b -p $HOME/miniconda
source "$HOME/miniconda/etc/profile.d/conda.sh"

# Create a conda environment
conda create -n bioinfo -y
conda activate bioinfo

# Install SAMtools, WhatsHap, Nanofilt from Bioconda
conda install -c bioconda samtools whatshap nanofilt -y

# Clone and install BWA Aligner
git clone https://github.com/lh3/bwa.git
cd bwa
make
cd ..

# Clone and install minimap2
git clone https://github.com/lh3/minimap2.git
cd minimap2
make
cd ..

# Clone and install fastq-filter
git clone https://github.com/LUMC/fastq-filter.git
cd fastq-filter
make
cd ..

# Clean up
rm ~/miniconda.sh

