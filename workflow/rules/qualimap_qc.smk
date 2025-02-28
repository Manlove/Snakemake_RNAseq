rule qualimap_bamqc:
    input:
        bam = "results/BAM/{sample}.sorted.bam",
        index = "results/BAM/{sample}.sorted.bam.bai",
        wd = os.getcwd()
    output:
        directory('results/QC/qualimap/{sample}_BAMQC/')
    threads: 8
    log:
        "results/logs/qualimap/{sample}.qualimap.bamqc.log"
    params:
        gtf = config["annotation_file"]
    shell:
        "(docker run --rm --user root "
            "-v {input.wd}/{input.bam}:/{input.bam} "
            "-v {input.wd}/{input.index}:/{input.index} " 
            "-v {input.wd}/{output}/:/{output}/ "
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
        index = "results/BAM/{sample}.sorted.bam.bai",
        wd = os.getcwd()
    output:
        directory('results/QC/qualimap/{sample}_RNAseqQC/'),
        'results/count_data/{sample}_countfile.rnaseqqc.txt'
    threads: 8
    log:
        "results/logs/qualimap/{sample}.qualimap.rnaseqqc.log"
    params:
        gtf = config["annotation_file"]
    shell:
        "(docker run --rm --user root "
            "-v {input.wd}/{input.bam}:/{input.bam} "
            "-v {input.wd}/{input.index}:/{input.index} " 
            "-v {input.wd}/{output[0]}/:/{output[0]}/ "
            "-v {params.gtf}:/{params.gtf} "
            "-v {input.wd}/results/count_data/:/results/count_data/ "
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
        wd = os.getcwd()
    output:
        'results/count_data/{sample}_countfile.compcount.txt'
    threads: 8
    log:
        "results/logs/qualimap/{sample}.qualimap.compcount.log"
    params:
        gtf = config["annotation_file"]
        
    shell:
        "(docker run --rm --user root "
            "-v {input.wd}/{input.bam}:/{input.bam} "
            "-v {input.wd}/{input.index}:/{input.index} " 
            "-v {params.gtf}:/{params.gtf} "
            "-v {input.wd}/results/count_data/:/results/count_data/ "
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
        sample_names = [sample for sample in config["samples"]],
        condition = [config["samples"][sample]["condition"] for sample in config["samples"]],
        count_file = [f"results/count_data/{sample}_countfile.compcount.txt" for sample in config["samples"]]
    output:
        "metadata.txt"
    log:
        "results/logs/qualimap/qualimap.datafile.log"
    script:
        "scripts/create_info_table.py"


rule qualimap_info_file:
    input:
        gtf_file = config["annotation_file"],
        reference = config["reference"]
    output:
        "info_file.txt"
    log:
        "results/logs/qualimap/qualimap.infofile.log"
    script:
        "scripts/createQualimapInfoFile.py -g {input.gtf_file} -f {input.reference} -o {output} 2> {log}"
        
rule qualimap_counts_QC:
    input:
        data_file = "/metadata.txt",
        index = "results/BAM/{sample}.sorted.bam.bai",
    output:
        directory('results/QC/qualimap/{sample}_CountsQC/')
    threads: 8
    log:
        "results/logs/qualimap/{sample}.qualimap.compcount.log"
    params:
        gtf = config["annotation_file"],
        wd = os.getcwd()
    shell:
        "(docker run --rm --user root "
            "-v {input.wd}/{input.bam}:/{input.bam} "
            "-v {input.wd}/{input.index}:/{input.index} " 
            "-v {params.gtf}:/{params.gtf} "
            "-v {input.wd}/results/count_data/:/results/ "
            "-w / "
            "community.wave.seqera.io/library/qualimap:2.3--8375b60bba97a2a6 "
                "bash -c \"mkdir -p /{output} && "
                "qualimap counts "
                    "--compare "
                    "--data {input.data_file} "
                    "--info {input.info_file} "
                    "--outdir {output} \" "
                    ")2> {log}"