import argparse
import numpy as np
import pandas as pd

def convert_npy_to_txt(npy_file, txt_file, pickle_file, lioness_edges):
    data = np.load(npy_file)

    # Load metadata (rownames & colnames) from pickle file
    lioness_df = pd.read_pickle(pickle_file)

    # Save edge names as text file, one name per line
    with open(lioness_edges, "w") as f:
        for row_name in lioness_df.index:
            f.write(f"{row_name}\n")

    
    np.savetxt(txt_file, data, fmt='%.6g', delimiter='\t')

input_file = snakemake.input.lioness
output_file = snakemake.output.lioness_txt
pickle_file = snakemake.input.lioness_pickle
edges_file = snakemake.output.lioness_edges

convert_npy_to_txt(input_file, output_file, pickle_file, edges_file)