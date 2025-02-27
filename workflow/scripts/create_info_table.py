with open(snakemake.output[0], 'w') as out_file:
    for name, condition, file_path in zip(snakemake.output[0], snakemake.output[1], snakemake.output[2]):
        out_file.write(f"{name}\t{condition}\t{file_path}\t2\n")