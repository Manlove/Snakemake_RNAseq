with open(snakemake.output[0], 'w') as out_file:
    out_file.write("sample/tcondition/n")
    for name, condition in zip(snakemake.params[0], snakemake.params[1]):
        out_file.write(f"{name}\t{condition}\n")