version 1.0

workflow WES_MuTect2_pipeline {
    input {
        String tumor_sample_name
        File tumor_bam
        File tumor_bam_idx
        String normal_sample_name
        File normal_bam
        File normal_bam_idx
        File pon_vcf
        File pon_vcf_tbi
        File ref_fasta
        File ref_fasta_idx
        File ref_fasta_dict
        File intervals_bed
        File exclude_intervals_bed
        File gnomad_vcf
        File gnomad_vcf_tbi
    }

    # 1. MuTect2
    call MuTect2 {
        input:
            tumor_sample_name = tumor_sample_name,
            tumor_bam = tumor_bam,
            tumor_bam_idx = tumor_bam_idx,
            normal_sample_name = normal_sample_name,
            normal_bam = normal_bam,
            normal_bam_idx = normal_bam_idx,
            ref_fasta = ref_fasta,
            ref_fasta_idx = ref_fasta_idx,
            ref_fasta_dict = ref_fasta_dict,
            gnomad_vcf = gnomad_vcf,
            gnomad_vcf_tbi = gnomad_vcf_tbi,
            pon_vcf = pon_vcf,
            pon_vcf_tbi = pon_vcf_tbi,
            intervals_bed = intervals_bed,
            exclude_intervals_bed = exclude_intervals_bed  
    }

    # 2. LearnReadOrientationModel (optional)
    call LearnReadOrientationModel {
        input:
            tumor_sample_name = tumor_sample_name,
            f1r2_tar_gz = MuTect2.output_f1r2_tar_gz
    }

    # 3. GetPileupSummaries (tumor)
    call GetPileupSummaries as tumor_GetPileupSummaries {
        input:
            sample_name = tumor_sample_name,
            bam = tumor_bam,
            bam_idx = tumor_bam_idx,
            gnomad_vcf = gnomad_vcf,
            gnomad_vcf_tbi = gnomad_vcf_tbi,
            intervals_bed = intervals_bed,
            exclude_intervals_bed = exclude_intervals_bed
    }

    # 4. GetPileupSummaries (normal)
    call GetPileupSummaries as normal_GetPileupSummaries {
        input:
            sample_name = normal_sample_name,
            bam = normal_bam,
            bam_idx = normal_bam_idx,
            gnomad_vcf = gnomad_vcf,
            gnomad_vcf_tbi = gnomad_vcf_tbi,
            intervals_bed = intervals_bed,
            exclude_intervals_bed = exclude_intervals_bed
    }

    # 5. CalculateContamination on tumor
    call CalculateContamination {
        input:
            tumor_sample_name = tumor_sample_name,
            tumor_pileup_summaries = tumor_GetPileupSummaries.output_pileup_summaries,
            normal_pileup_summaries = normal_GetPileupSummaries.output_pileup_summaries
    }

    # 6. FilterMutectCalls
    call FilterMutectCalls {
        input:
            tumor_sample_name = tumor_sample_name,
            tumor_bam = tumor_bam,
            tumor_bam_idx = tumor_bam_idx,
            normal_sample_name = normal_sample_name,
            normal_bam = normal_bam,
            normal_bam_idx = normal_bam_idx,
            ref_fasta = ref_fasta,
            ref_fasta_idx = ref_fasta_idx,
            ref_fasta_dict = ref_fasta_dict,
            gnomad_vcf = gnomad_vcf,
            gnomad_vcf_tbi = gnomad_vcf_tbi,
            pon_vcf = pon_vcf,
            pon_vcf_tbi = pon_vcf_tbi,
            intervals_bed = intervals_bed,
            exclude_intervals_bed = exclude_intervals_bed,
            mutect2_unfiltered_vcf_stats = MuTect2.output_vcf_stats,
            artifact_priors_tar_gz = LearnReadOrientationModel.output_artifact_priors_tar_gz,
            contamination_table = CalculateContamination.output_contamination_table,
            tumor_segments_table = CalculateContamination.output_segments_table
    }

    output {
        File mutect2_vcf = MuTect2.output_vcf
        File mutect2_vcf_idx = MuTect2.output_vcf_idx
        File mutect2_vcf_stats = MuTect2.output_vcf_stats
        File mutect2_f1r2_tar_gz = MuTect2.output_f1r2_tar_gz
        File mutect2_bamout = MuTect2.output_bamout
        File artifact_priors_tar_gz = LearnReadOrientationModel.output_artifact_priors_tar_gz
        File tumor_pileup_summaries = tumor_GetPileupSummaries.output_pileup_summaries
        File normal_pileup_summaries = normal_GetPileupSummaries.output_pileup_summaries
        File tumor_segments_table = CalculateContamination.output_segments_table
        File tumor_contamination_table = CalculateContamination.output_contamination_table
        File mutect2_filtered_vcf = FilterMutectCalls.output_filtered_vcf
        File mutect2_filtered_vcf_idx = FilterMutectCalls.output_filtered_vcf_idx
        File mutect2_filtered_vcf_stats = FilterMutectCalls.output_filtered_vcf_stats
    }
}

## TASK DEFINITIONS ###########################################################################################

# TODO: separate from tumor-only and matched-normal modes
task MuTect2 {
    input {
        String tumor_sample_name
        File tumor_bam
        File tumor_bam_idx
        String normal_sample_name
        File normal_bam
        File normal_bam_idx
        File ref_fasta
        File ref_fasta_idx
        File ref_fasta_dict
        # genome aggregation database allows distinguish between germline and somatic variants
        File gnomad_vcf
        File gnomad_vcf_tbi
        # exclude sequencing artifcacts and technical errors
        File pon_vcf
        File pon_vcf_tbi
        # Specify genomic regions of interests, restricts to targeted exome regions
        File intervals_bed
        # Blacklist regions to exclude from analysis (centromeres, telomeres, segmental duplications, high gc)
        File exclude_intervals_bed

        # Configurable
        Int memory
        Int disk_space
        Int num_threads
        Int num_preempt
    }

    Int java_memory = ceil(memory * 0.8)
    
    #--genotype-germline-sites false \
    #--genotype-pon-sites false \
    command {
        set -e
        export GATK_LOCAL_JAR="/root/gatk.jar"
        gatk --java-options "-Xmx~{java_memory}g" Mutect2 \
            -R ~{ref_fasta} \
            -I ~{tumor_bam} \
            -I ~{normal_bam} \
            --normal-sample ~{normal_sample_name} \
            --germline-resource ~{gnomad_vcf} \
            --panel-of-normals ~{pon_vcf} \
            --intervals ~{intervals_bed} \
            --exclude-intervals ~{exclude_intervals_bed} \
            --f1r2-tar-gz ~{tumor_sample_name}.f1r2.tar.gz \
            --bamout ~{tumor_sample_name}.bamout.bam \
            -O ~{tumor_sample_name}.unfiltered.vcf \
            --smith-waterman FASTEST_AVAILABLE \
            --pair-hmm-implementation FASTEST_AVAILABLE \
            --native-pair-hmm-threads ~{num_threads} \
            --dont-use-soft-clipped-bases false \
            --linked-de-bruijn-graph true \
            --recover-all-dangling-branches true \
            --pileup-detection true \
            --native-pair-hmm-use-double-precision false
    }

     output {
        File output_vcf = "~{tumor_sample_name}.unfiltered.vcf"
        File output_vcf_idx = "~{tumor_sample_name}.unfiltered.vcf.idx"
        File output_vcf_stats = "~{tumor_sample_name}.unfiltered.vcf.stats"
        File output_f1r2_tar_gz = "~{tumor_sample_name}.f1r2.tar.gz"
        File output_bamout = "~{tumor_sample_name}.bamout.bam"
    }

    runtime {
        docker: "broadinstitute/gatk:4.6.0.0"
        memory: "~{memory} GB"
        disks: "local-disk ~{disk_space} HDD"
        cpu: "~{num_threads}"
        preemptible: "~{num_preempt}"
    }
}


## OPTIONAL: for oxidative artifacts
task LearnReadOrientationModel{
    input {
        String tumor_sample_name
        File f1r2_tar_gz

        # Configurable
        Int memory
        Int disk_space
        Int num_threads
        Int num_preempt
    }

    Int java_memory = ceil(memory * 0.8)

    command {
        set -e
        export GATK_LOCAL_JAR="/root/gatk.jar"
        gatk --java-options "-Xmx~{java_memory}g" LearnReadOrientationModel \
            -I ~{f1r2.tar.gz} \
            -O ~{tumor_sample_name}_artifact_priors.tar.gz
    }

    output {
        File output_artifact_priors_tar_gz = "~{tumor_sample_name}_artifact_priors.tar.gz"
    }

    runtime {
        docker: "broadinstitute/gatk:4.6.0.0"
        memory: "~{memory} GB"
        disks: "local-disk ~{disk_space} HDD"
        cpu: "~{num_threads}"
        preemptible: "~{num_preempt}"
    }
}

# Use on both normal and tumor samples
# --read-filter FirstOfPairReadFilter?
# --read-filter PairendReadFilter?
task GetPileupSummaries{
    input {
        String sample_name
        File bam
        File bam_idx
        # genome aggregation database allows distinguish between germline and somatic variants
        File gnomad_vcf
        File gnomad_vcf_tbi
        # Specify genomic regions of interests, restricts to targeted exome regions
        File intervals_bed
        # Blacklist regions to exclude from analysis (centromeres, telomeres, segmental duplications, high gc)
        File exclude_intervals_bed

        # Configurable
        Int memory
        Int disk_space
        Int num_threads
        Int num_preempt
    }

    command {
        set -e
        export GATK_LOCAL_JAR="/root/gatk.jar"
        gatk --java-options "-Xmx~{java_memory}g" GetPileupSummaries \
            -I ~{bam} \
            --variant ~{gnomad_vcf} \
            --intervals ~{intervals_bed} \
            --exclude-intervals ~{exclude_intervals_bed} \
            -O ~{sample_name}_getpileupsummaries.table
    }

    output {
        File output_pileup_summaries = "~{sample_name}.getpileupsummaries.table"
    }
    
    runtime {
        docker: "broadinstitute/gatk:4.6.0.0"
        memory: "~{memory} GB"
        disks: "local-disk ~{disk_space} HDD"
        cpu: "~{num_threads}"
        preemptible: "~{num_preempt}"
    }
}

# TODO: separate tumor-only mode and matched-normal mode
task CalculateContamination{
    input {
        String tumor_sample_name
        File tumor_pileup_summaries
        File normal_pileup_summaries

        # Configurable
        Int memory
        Int disk_space
        Int num_threads
        Int num_preempt
    }

    Int java_memory = ceil(memory * 0.8)

    command {
        set -e
        export GATK_LOCAL_JAR="/root/gatk.jar"
        gatk --java-options "-Xmx~{java_memory}g" CalculateContamination \
            -I ~{tumor_pileup_summaries} \
            --matched-normal ~{normal_pileup_summaries} \
            --tumor-segmentation ~{tumor_sample_name}_segments.table \
            -O ~{tumor_sample_name}_calculatecontamination.table
    }

    output {
        File output_segments_table = "~{tumor_sample_name}_segments.table"
        File output_contamination_table = "~{tumor_sample_name}_calculatecontamination.table"
    }
    
    runtime {
        docker: "broadinstitute/gatk:4.6.0.0"
        memory: "~{memory} GB"
        disks: "local-disk ~{disk_space} HDD"
        cpu: "~{num_threads}"
        preemptible: "~{num_preempt}"
    }
}


# --max-median-fragment-length-difference = 10000  # default: 10000
# --min-median-base-quality = 20  # default: 20
# --min-median-mapping-quality = 20  # default: -1
# --min-median-read-position = 5  # default: 1
task FilterMutectCalls{
    input {
        String tumor_sample_name
        File tumor_bam
        File tumor_bam_idx
        String normal_sample_name
        File normal_bam
        File normal_bam_idx
        File ref_fasta
        File ref_fasta_idx
        File ref_fasta_dict
  
        # genome aggregation database allows distinguish between germline and somatic variants
        File gnomad_vcf
        File gnomad_vcf_tbi
        # exclude sequencing artifcacts and technical errors
        File pon_vcf
        File pon_vcf_tbi
        # Specify genomic regions of interests, restricts to targeted exome regions
        File intervals_bed
        # Blacklist regions to exclude from analysis (centromeres, telomeres, segmental duplications, high gc)
        File exclude_intervals_bed

        File mutect2_unfiltered_vcf_stats # From MuTect2
        File artifact_priors_tar_gz # From LearnReadOrientationModel
        File contamination_table # From CalculateContamination
        File tumor_segments_table # From CalculateContamination
        
        # Configurable
        Int memory
        Int disk_space
        Int num_threads
        Int num_preempt
    }

    Int java_memory = ceil(memory * 0.8)

    command {
        set -e
        export GATK_LOCAL_JAR="/root/gatk.jar"
        gatk --java-options "-Xmx~{java_memory}g" FilterMutectCalls \
            --reference ~{ref_fasta} \
            --variant ~{gnomad_vcf} \
            -O ~{tumor_sample_name}.filtered.vcf \
            --filtering-stats ~{tumor_sample_name}.filtered.vcf.stats \
            --stats ~{mutect2_unfiltered_vcf_stats} \
            --orientation-bias-artifact-priors ~{artifact_priors_tar_gz} \
            --contamination-table ~{contamination_table} \
            --tumor-segmentation ~{tumor_segments_table} \
            --max-median-fragment-length-difference = 10000 \
            --min-median-base-quality = 20 \
            --min-median-mapping-quality = 20 \
            --min-median-read-position = 5 
    }

    output {
        File output_filtered_vcf = "~{tumor_sample_name}.filtered.vcf"
        File output_filtered_vcf_idx = "~{tumor_sample_name}.filtered.vcf.idx"
        File output_filtered_vcf_stats = "~{tumor_sample_name}.filtered.vcf.stats"
    }

    runtime {
        docker: "broadinstitute/gatk:4.6.0.0"
        memory: "~{memory} GB"
        disks: "local-disk ~{disk_space} HDD"
        cpu: "~{num_threads}"
        preemptible: "~{num_preempt}"
    }
}