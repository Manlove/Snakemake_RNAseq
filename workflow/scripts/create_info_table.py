with open(snakemake.output[0], 'w') as out_file:
    for name, condition, file_path in zip(snakemake.params[0], snakemake.params[1], snakemake.input):
        out_file.write(f"{name}\t{condition}\t{file_path}\t2\n")