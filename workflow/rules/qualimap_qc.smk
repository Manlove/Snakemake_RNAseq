
rule qualimap_bamqc:
    input:
        bam = "results/BAM/{sample}.sorted.bam",
        index = "results/BAM/{sample}.sorted.bam.bai"
    output:
        directory('results/QC/qualimap/{sample}_BAMQC/')
    threads: 8
    log:
        "results/logs/qualimap/{sample}.qualimap.bamqc.log"
    params:
        gtf = config["annotation_file"],
        wd = config["wd"]
    shell:
        "(docker run --rm --user root "
            "-v {params.wd}/{input.bam}:/{input.bam} "
            "-v {params.wd}/{input.index}:/{input.index} " 
            "-v {params.wd}/{output}/:/{output}/ "
            "-v {params.gtf}:/{params.gtf} "
            "-w / "
            "community.wave.seqera.io/library/qualimap:2.3--8375b60bba97a2a6 "
            "bash -c \"mkdir -p /{output} && "
            "qualimap bamqc "
            "-bam {input.bam} "
            "--feature-file {params.gtf} "
            "--paint-chromosome-limits "
            "--nt {threads} "
            "-outdir /{output}/ "
            "-outformat 'html'\" "
            ")2> {log}"

rule qualimap_rnaseqqc:
    input:
        bam = "results/BAM/{sample}.sorted.bam",
        index = "results/BAM/{sample}.sorted.bam.bai"
    output:
        directory('results/QC/qualimap/{sample}_RNAseqQC/'),
        'results/count_data/{sample}_countfile.rnaseqqc.txt'
    threads: 8
    log:
        "results/logs/qualimap/{sample}.qualimap.rnaseqqc.log"
    params:
        gtf = config["annotation_file"],
        wd = config["wd"]
    shell:
        "(docker run --rm --user root "
            "-v {params.wd}/{input.bam}:/{input.bam} "
            "-v {params.wd}/{input.index}:/{input.index} " 
            "-v {params.wd}/{output[0]}/:/{output[0]}/ "
            "-v {params.gtf}:/{params.gtf} "
            "-v {params.wd}/results/count_data/:/results/count_data/ "
            "-w / "
            "community.wave.seqera.io/library/qualimap:2.3--8375b60bba97a2a6 "
            "bash -c \"mkdir -p /{output[0]} && "
            "qualimap rnaseq "
            "-bam {input.bam} "
            "-gtf {params.gtf} "
            "--paired "
            "--sorted "
            "-oc /{output[1]} "
            "-outdir /{output[0]}/ "
            "-outformat 'html'\" "
            ")2> {log}"


rule qualimap_compute_counts:
    """This step can be done with the RNASeqQC function"""
    input:
        bam = "results/BAM/{sample}.sorted.bam",
        index = "results/BAM/{sample}.sorted.bam.bai",

    output:
        'results/count_data/{sample}_countfile.compcount.txt'
    threads: 8
    log:
        "results/logs/qualimap/{sample}.qualimap.compcount.log"
    params:
        gtf = config["annotation_file"],
        wd = config["wd"]
        
    shell:
        "(docker run --rm --user root "
            "-v {params.wd}/{input.bam}:/{input.bam} "
            "-v {params.wd}/{input.index}:/{input.index} " 
            "-v {params.gtf}:/{params.gtf} "
            "-v {params.wd}/results/count_data/:/results/count_data/ "
            "-w / "
            "community.wave.seqera.io/library/qualimap:2.3--8375b60bba97a2a6 "
            "bash -c \"mkdir -p /results/count_data/ && "
            "qualimap comp-counts "
            "-bam {input.bam} "
            "-gtf {params.gtf} "
            "-pe "
            "--sorted "
            "-out /{output}/\" "
            ")2> {log}"


rule qualimap_data_file:
    input:
        count_file = [f"results/count_data/{sample}_countfile.rnaseqqc.txt" for sample in config["samples"]]
    output:
        "metadata.txt"
    log:
        "results/logs/qualimap/qualimap.datafile.log"
    params:
        sample_names = [sample for sample in config["samples"]],
        condition = [config["samples"][sample]["condition"] for sample in config["samples"]]
    script:
        "../scripts/create_info_table.py"


rule qualimap_info_file:
    input:
        gtf_file = config["annotation_file"],
        reference = config["reference"]
    output:
        "info_file.txt"
    log:
        "results/logs/qualimap/qualimap.infofile.log"
    conda:
        "../envs/createQualimapInfoFile.yaml"
    script:
        "../scripts/createQualimapInfoFile.py"


rule qualimap_counts_QC:
    input:
        data_file = "metadata.txt",
        info_file = "info_file.txt"
    output:
        directory('results/QC/qualimap/CountsQC/')
    threads: 8
    log:
        "results/logs/qualimap/qualimap.countqc.log"
    params:
        wd = config["wd"]   
    shell:
        "(docker run --rm --user root "
            "-v {params.wd}/{input.data_file}:/{input.data_file} "
            "-v {params.wd}/{input.info_file}:/{input.info_file} " 
            "-v {params.wd}/{output}/:/{output}/ "
            "-v {params.wd}/results/count_data/:/results/count_data/ "
            "-w / "
            "community.wave.seqera.io/library/qualimap:2.3--8375b60bba97a2a6 "
                "bash -c \"mkdir -p /{output} && "
                "qualimap counts "
                    "--compare "
                    "--data {input.data_file} "
                    "--info {input.info_file} "
                    "--outdir {output} "
                    "&& ls /results/QC/qualimap/ \" "
                    ")2> {log}"