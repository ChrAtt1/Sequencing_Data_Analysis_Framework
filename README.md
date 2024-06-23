
# Sequencing Data Analysis Framework

This repository contains scripts and tools for processing and analyzing sequencing data from Illumina and ONT (Oxford Nanopore Technologies) platforms.

## Folder Structure

### 0_sample_data
- `1_data_processing`
  - Sample output as .vcf of the data processing pipeline for Illumina and ONT for 10 sample patients.
- `2_data_analyzer`
  - `2024-06-21_step_1_vcf_import_exons_Illumina`
  - `2024-06-21_step_1_vcf_import_exons_ONT`
  - `2024-06-21_step_2_Illumina_and_ONT_merged`
  - `2024-06-21_step_3_dataset_analyzed`
  - `2024-06-21_step_4_dataset_color_coded`
- `reference_sequence`
  - `amplicon_reference_sequence.fa`

### 1_data_processing
- `install_dependencies.sh`
- `sequ_data_framework_Illumina.sh`
- `sequ_data_framework_ONT.sh`

### 2_genotype_analyzer
- `GenotypeAnalyzer.exe`
- `GenotypeAnalyzer.py`
- `Setting_Genotype-Analyzer.xlsx`
- `step_1_build_genotype.py`
- `step_2_merge_files.py`
- `step_3_analyse_file.py`
- `step_4_color_count.py`

## Sequencing Data Framework Setup Guide

### Operating System Requirement
- OS required: Linux

### Installation Instructions

1. **Clone the Repository:**

    ```bash
    git clone https://github.com/ChrAtt1/Sequencing-Data-Analysis-Framework.git
    cd Sequencing-Data-Analysis-Framework
    ```

2. **Install Dependencies:**

    ```bash
    chmod +x ./install_dependencies.sh
    ./install_dependencies.sh
    ```

    The `install_dependencies.sh` script will install the following:
    - Conda
    - SAMtools
    - BWA Aligner
    - WhatsHap
    - Nanofilt
    - minimap2
    - fastq-filter

    *(Links for installation guides and repositories were last accessed on 23 June 2024.)*

### Initial Configuration

Before the first start, modify the permissions of the scripts to make them executable:

```bash
chmod +x ./sequ_data_framework_Illumina.sh
chmod +x ./sequ_data_framework_ONT.sh
```

Insert the variables directly into the shell scripts (`sequ_data_framework_Illumina.sh` and `sequ_data_framework_ONT.sh`). Assign your specific file paths to these variables. Example:

```bash
#!/bin/bash
# Specify the base path where the data is located
path_base_data="/path/to/your/base/data"
# Specify the path to the input data
input_data="/path/to/your/input/data"
# Provide the path to the reference sequence file (.fa-file)
path_reference_sequence="/path/to/your/reference/sequence.fa"
```

Replace the placeholder paths (/path/to/your/...) with the actual paths on your system.

### Running the Scripts

To run the scripts, use the following commands in the terminal:

```bash
./sequ_data_framework_Illumina.sh
./sequ_data_framework_ONT.sh
```

Ensure you have the necessary permissions and that the paths specified in the scripts are correct before execution.

# Genotype Analyzer

Genotype Analyzer is a Python-based application designed to analyze genotype data using various settings and methods. The application provides a graphical user interface (GUI) built with Tkinter, enabling users to input necessary files, configure settings, and perform genotype analysis.

## Features

- Load and parse settings from an Excel file
- Browse and select files and folders through the GUI
- Configure sequencing methods and amplicon settings
- Validate input paths and settings
- Perform genotype building, merging, and analysis
- Display progress with a progress bar
- Measure and print execution time for each step

## Requirements

- Python 3.8+
- Pandas
- Tkinter
- Linux
- Pycharm

## Installation

### Install Dependencies

Use `pip` to install the required libraries:

```bash
pip install pandas
pip install math
pip install collections
pip install datetime
pip install scipy
pip install numpy
pip install re
pip install shutil
pip install openpyxl
pip install tkinter
```

## Running the Application

### For Linux

```bash
python GenotypeAnalyzer.py
```

### For Windows

Run the `GenotypeAnalyzer.exe` file.

## Usage

1. **Load the Settings File:**

    The application expects an Excel file named `Setting_Genotype-Analyzer.xlsx` with two sheets:
    - `Sequencing Method Setting`
    - `Amplicon Setting`

2. **Configure Settings:**

    - Use the GUI to browse and select the necessary files and folders.
    - Configure the sequencing methods and amplicon settings as needed.
    - Choose between analyzing VCF files (Option A) or using previous Genotype Analyzer output Excel files (Option B).

3. **Run the Analysis:**

    - Click the "Analyse Genotypes" button to start the analysis process.
    - The application will validate the input paths and settings before proceeding.
    - Progress is displayed with a progress bar, and execution times for each step are printed in the console.
