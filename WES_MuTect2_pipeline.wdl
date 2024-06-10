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
        File pon_vcf
        File ref_fasta
        File ref_fasta_idx
        File ref_fasta_dict
        File intervals_bed
        File exclude_intervals_bed
        File gnomad
        File gnomad_tbi

        # Configurable
        Int memory
        Int disk_space
        Int num_threads
        Int num_preempt
    }

    Int java_memory = ceil(memory * 0.8)

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
            -O ~{tumor_sample_name}.vcf
    }

     output {
        File output_vcf = "~{tumor_sample_name}.vcf"
        File output_vcf_idx = "~{tumor_sample_name}.vcf.idx"
        File output_vcf_stats = "~{tumor_sample_name}.vcf.stats"
    }

    runtime {
        docker: "broadinstitute/gatk:4.5.0.0"
        memory: "${memory} GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${num_threads}"
        preemptible: "${num_preempt}"
    }
}