with open(snakemake.output[0], 'w') as out_file:
    out_file.write("sample,condition\n")
    for name, condition in zip(snakemake.params[0], snakemake.params[1]):
        out_file.write(f"results.BAM.{name}.sorted.bam,{condition}\n")