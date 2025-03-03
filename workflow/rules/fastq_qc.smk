rule fastqc_pre_trim:
    input:
        lambda wildcards: config["samples"][wildcards.sample]["reads"]
    output:
        directory("results/QC/FastQC/{sample}/")
    conda:
        "../envs/qc.yaml"
    threads: 2
    log:
        "results/logs/FastQC/{sample}.log"
    shell:
        "mkdir -p {output}; "
        "fastqc {input} --outdir {output} --threads {threads} 2> {log}"


rule fastp_trim: 
    input:
        lambda wildcards: config["samples"][wildcards.sample]["reads"]
    output:
        read1 = "fastq/trimmed/{sample}_trimmed_1.fastq.gz",
        read2 = "fastq/trimmed/{sample}_trimmed_2.fastq.gz",
        html = "results/QC/fastp/{sample}_fastp_report.html",
        json = "results/QC/fastp/{sample}_fastp_report.json"
    conda:
        "../envs/qc.yaml"
    log:
        "results/logs/fastp/{sample}.log"
    threads: 4
    params:
        min_qual = config.get("minimum_quality",15),
        max_unqual = config.get("max_unqualified",40),
        average_min = config.get("average_minimum",0),
        mean_window_qual = config.get("mean_quality",20),
        cut_window = config.get("window_size",4),
        min_length = config.get("minimum_length",50),
        max_length = config.get("maximum_length",0),
        adapter = ("--detect_adapter_for_pe" if config.get("auto_detect_adapters", False) else (f"--adapter_sequence {config['adapter_read_1']} --adapter_sequence_r2 {config['adapter_read_2']}" if (config['adapter_read_1'] != "" and config['adapter_read_2'] != "") else "")),
        complexity_option = f"--low_complexity_filter --complexity_threshold {config.get('complexity_minimum',30)}" if config.get("complexity_filter",False) else "",
        disable_filtering = f"--disable_quality_filtering" if config.get("disable_filtering",False) else ""
    shell:
        "fastp -i {input[0]} -I {input[1]} "
            "-o {output.read1} "
            "-O {output.read2} "
            "{params.adapter} "
            "{params.complexity_option} "
            "--thread {threads} "
            "--cut_front "
            "--cut_tail "
            "--cut_window_size {params.cut_window} "
            "--cut_mean_quality {params.mean_window_qual} "
            "{params.disable_filtering} "
            "--qualified_quality_phred {params.min_qual} "
            "--unqualified_percent_limit {params.max_unqual} "
            "--average_qual {params.average_min} "
            "--length_required={params.min_length} "
            "--length_limit={params.max_length} "
            "--html {output.html} "
            "--json {output.json} "
            "2> {log}"

 
rule fastqc_post_trim:
    input:
        "fastq/trimmed/{sample}_trimmed_1.fastq.gz",
        "fastq/trimmed/{sample}_trimmed_2.fastq.gz"
    output:
        directory("results/QC/FastQC/{sample}_trimmed/")
    conda:
        "../envs/qc.yaml"
    threads: 2
    log:
        "results/logs/FastQC/{sample}_trimmed.log"
    shell:
        "mkdir -p {output}; "
        "fastqc {input} --outdir {output} --threads {threads} 2> {log}"


# rule fastq_screen:
#     input:
#         pass