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

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert npy file to txt.")
    parser.add_argument("--input", required=True, help="Path to the input npy file.")
    parser.add_argument("--output", required=True, help="Path to the output file.")
    parser.add_argument("--pickle", required=True, help="Path to the pickle file containing metadata.")
    parser.add_argument("--edges", required=True, help="Path to the output edges text file.")
    args = parser.parse_args()

    convert_npy_to_txt(args.input, args.output, args.pickle, args.edges)