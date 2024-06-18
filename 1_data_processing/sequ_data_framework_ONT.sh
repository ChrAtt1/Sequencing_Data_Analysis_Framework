#!/bin/bash

# Activate conda environment for samtools, bcftools, and whatshap
eval "$(command conda 'shell.bash' 'hook' 2> /dev/null)"

# Function to process reads: align, convert SAM to BAM, sort, index, and variant call
process_reads() {
    
    local file_base_name=$(basename "$cur_file" .fastq.gz)
    
    local nanofilt_output="$output_data_filtered/${file_base_name}_filtered.fastq.gz"
    local minimap_output="$output_data_aligned/${file_base_name}.sam"
    local sam_output="$output_data_aligned/${file_base_name}.sam"
    local bam_output_base="$output_data_converted/${file_base_name}"
    local bam_output="$bam_output_base.bam"
    local bam_output_sorted="$bam_output_base_sorted.bam"
    local bcf_output="$output_data_variant_called/${file_base_name}_called.bcf"
    local bcf_output_highqual="$output_data_variant_called/${file_base_name}_called-highqual.bcf"
    local vcf_output_highqual="$output_data_variant_called/${file_base_name}_called-highqual.vcf"
    local vcf_output_highqual_gz="$vcf_output_highqual.gz"
    local whatshap_output_vcf="$output_data_phased/${file_base_name}_phased.vcf"
    local whatshap_output_vcf_gz="$whatshap_output_vcf.gz"
    local whatshap_output_gtf="$output_data_phased/${file_base_name}_called-highqual.gtf"
    local whatshap_output_stats="$output_data_phased/${file_base_name}_statistics.ods"
    local whatshap_output_haplotagged_bam="$output_data_phased/${file_base_name}_haplotagged.bam"
    local bcf_output_haplotype1_fasta="$output_data_haplotypes/${file_base_name}_haplotype1.fasta"
    local bcf_output_haplotype2_fasta="$output_data_haplotypes/${file_base_name}_haplotype2.fasta"
    
    conda activate /home/smink/anaconda3/envs/Filtering
    zcat $cur_file|NanoFilt -q 8 -l 5000 --maxlength 15000| gzip > $nanofilt_output

    conda deactivate
    conda activate /home/smink/anaconda3/envs/medaka
    minimap2 -ax map-ont $path_reference_sequence $nanofilt_output > $minimap_output
    conda deactivate
    
    conda activate /home/smink/anaconda3/envs/samtools
    samtools view -S -b $minimap_output > "$bam_output"
    samtools sort "$bam_output" -o "$bam_output_sorted"
    samtools index "$bam_output_sorted"
    
    # Perform variant calling using bcftools
    bcftools mpileup --max-depth 1000 -f "$path_reference_sequence" "$bam_output_sorted" | bcftools call -P 1.1e-3 -mv -Ob -o "$bcf_output"
    bcftools view -i '%QUAL>=10' "$bcf_output" -o "$bcf_output_highqual"
    bcftools convert -Ov -o "$vcf_output_highqual" "$bcf_output_highqual"
    bgzip -c "$vcf_output_highqual" > "$vcf_output_highqual_gz"
    tabix "$vcf_output_highqual_gz"
    conda deactivate
    
    # Whatshap phase
    whatshap phase -o $whatshap_output_vcf --ignore-read-groups --reference=$path_reference_sequence $vcf_output_highqual $bam_output_sorted

    # Whatshap stats
    whatshap stats --gtf=$whatshap_output_gtf $whatshap_output_vcf
    
    conda activate /home/smink/anaconda3/envs/samtools

    # Bcftools convert to gz
    bcftools convert -Oz -o "$whatshap_output_vcf_gz" "$whatshap_output_vcf"

    # Index the gzipped VCF
    bcftools index "$whatshap_output_vcf_gz"
    
    conda deactivate

    # Whatshap haplotag
    whatshap haplotag --ignore-read-groups -o "$whatshap_output_haplotagged_bam" --reference="$path_reference_sequence" "$whatshap_output_vcf_gz" "$bam_output_sorted"
    
    conda activate /home/smink/anaconda3/envs/samtools
    
    # Samtools index the haplotagged BAM
    samtools index "$whatshap_output_haplotagged_bam"
    conda deactivate

    # Whatshap stats for phased VCF and generate statistics
    whatshap stats "$whatshap_output_vcf" --tsv="$whatshap_output_stats"
    
    conda activate /home/smink/anaconda3/envs/samtools
    # Bcftool calls
    bcftools consensus -H 1 -f $path_reference_sequence $whatshap_output_vcf_gz > $bcf_output_haplotype1_fasta
    bcftools consensus -H 2 -f $path_reference_sequence $whatshap_output_vcf_gz > $bcf_output_haplotype2_fasta
    conda deactivate
    

    
}

# Function to create directory if not exists
ensure_dir() {
    [ ! -d "$1" ] && mkdir -p "$1"
}

# Define base paths
# Specify the base path where the data is located
path_base_data="/media/smink/Data_Linux/Sequencing_Data_Analysis_Framework/0_sample_data/1_data_processing/ONT"

# Specify the path to the input data
input_data="$path_base_data/0_raw_data"

# Specify the path to the output data
output_data="$path_base_data"

# Provide the path to the reference sequence file (.fa-file)
path_reference_sequence="/media/smink/Data_Linux/Sequencing_Data_Analysis_Framework/0_sample_data/0_reference_sequence/amplicon_reference_sequence.fa"



# Set output data paths
output_data_filtered="$output_data/1_filtered_dataset"
output_data_aligned="$output_data/2_mapped_dataset"
output_data_converted="$output_data/3_converted_dataset"
output_data_variant_called="$output_data/4_variant_called_dataset"
output_data_phased="$output_data/5_phased_dataset"
output_data_haplotypes="$output_data/6_haplotypes_dataset"

# Export reference sequence once
export path_reference_sequence

# Loop over each R1 file for processing
ensure_dir "$output_data_filtered"
ensure_dir "$output_data_aligned"
ensure_dir "$output_data_converted"
ensure_dir "$output_data_variant_called"
ensure_dir "$output_data_phased"
ensure_dir "$output_data_haplotypes"


# Process each fastq file in the input directory
for cur_file in "$input_data"/*.fastq.gz; do
    if [ -f "$cur_file" ]; then
        process_reads "$cur_file"
    else
        echo "File does not exist: $cur_file"
    fi
done
