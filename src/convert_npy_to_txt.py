import argparse
import numpy as np
import pandas as pd

def convert_npy_to_txt(npy_file, txt_file, pickle_file, lioness_edges):
    data = np.load(npy_file)

    # metadata
    lioness_df = pd.read_pickle(pickle_file)

    # optional separate edge-name file
    with open(lioness_edges, "w") as f:
        for row_name in lioness_df.index:
            f.write(f"{row_name}\n")

    # write matrix with both row and column names
    out_df = pd.DataFrame(data, index=lioness_df.index, columns=lioness_df.columns)
    out_df.to_csv(txt_file, sep="\t", index=True, header=True)

input_file = snakemake.input.lioness
output_file = snakemake.output.lioness_txt
pickle_file = snakemake.input.lioness_pickle
edges_file = snakemake.output.lioness_edges

convert_npy_to_txt(input_file, output_file, pickle_file, edges_file)
