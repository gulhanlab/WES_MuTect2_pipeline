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
        File gnomad
        File gnomad_tbi
    }

    call MuTect2 {
        input:
            tumor_sample_name = tumor_sample_name,
            tumor_bam = tumor_bam,
            tumor_bam_idx = tumor_bam_idx,
            normal_sample_name = normal_sample_name,
            normal_bam = normal_bam,
            normal_bam_idx = normal_bam_idx,
            pon_vcf = pon_vcf,
            pon_vcf_tbi = pon_vcf_tbi,
            ref_fasta = ref_fasta,
            ref_fasta_idx = ref_fasta_idx,
            ref_fasta_dict = ref_fasta_dict,
            intervals_bed = intervals_bed,
            exclude_intervals_bed = exclude_intervals_bed,
            gnomad = gnomad,
            gnomad_tbi = gnomad_tbi
    }

    output {
        File mutect2_vcf = MuTect2.output_vcf
        File mutect2_vcf_idx = MuTect2.output_vcf_idx
        File mutect2_vcf_stats = MuTect2.output_vcf_stats
    }
}

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
        File gnomad
        File gnomad_tbi
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
        gatk --java-options "-Xmx~{java_memory}g" Mutect2 \
            -R ~{ref_fasta} \
            -I ~{tumor_bam} \
            -I ~{normal_bam} \
            --normal-sample ~{normal_sample_name} \
            --germline-resource ~{gnomad} \
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
        memory: "${memory} GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${num_threads}"
        preemptible: "${num_preempt}"
    }
}


## OPTIONAL: for oxidative artifacts
# task LearnReadOrientationModel{
#     command {
#         gatk LearnReadOrientationModel -I f1r2.tar.gz -O read-orientation-model.tar.gz
#     }
# }

# task GetPileupSummaries{
#     command {
#         gatk GetPileupSummaries \
#             -I tumor.bam \
#             -V chr17_small_exac_common_3_grch38.vcf.gz \
#             -L chr17_small_exac_common_3_grch38.vcf.gz \
#             -O getpileupsummaries.table
#     }
# }

# task CalculateContamination{
#     command {
#         gatk CalculateContamination \
#             -I getpileupsummaries.table \
#             -tumor-segmentation segments.table \
#             -O calculatecontamination.table
#     }

# }

# task FilterMutectCalls{
#     command {
#         gatk FilterMutectCalls -V unfiltered.vcf \
#             [--tumor-segmentation segments.table] \
#             [--contamination-table contamination.table] \
#             --ob-priors read-orientation-model.tar.gz \
#             -O filtered.vcf
#     }

# }