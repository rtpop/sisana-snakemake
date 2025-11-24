import yaml
import argparse
import os

def generate_sisana_params(exp_file, motif_file, ppi_file, number, outdir, method,
                           pandafilepath, lionessfilepath, output_filename, ncores=10, input_format='csv'):
    # get the file name and append the preprocessed suffix... cos this is how sisana does things atm
    exp = os.path.join(outdir, exp_file.split('/')[-1].split('.')[0] + '_preprocessed.txt')

    data = {
        'preprocess': {
            'exp_file': exp_file,
            'filetype': input_format,
            'number': number,
            'outdir': outdir
        },
        'generate': {
            'exp': exp,
            'motif': motif_file,
            'ppi': ppi_file,
            'method': method,
            'pandafilepath': pandafilepath,
            'lionessfilepath': lionessfilepath,
            'ncores': ncores,
            'start': None,
            'end': None,
            'compute': "cpu"
        }
    }

    with open(output_filename, 'w') as file:
        yaml.dump(data, file, default_flow_style=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate SiSaNA params YAML file")
    parser.add_argument('--exp', required=True, help='Path to expression file')
    parser.add_argument('--motif', required=True, help='Path to motif file')
    parser.add_argument('--ppi', required=True, help='Path to PPI file')
    parser.add_argument('--number', type=int, required=True, help='Number of samples a gene must be expressed in')
    parser.add_argument('--outdir', required=True, help='Output directory for preprocess stage')
    parser.add_argument('--method', required=True, help='Method to use (panda or lioness)')
    parser.add_argument('--pandafilepath', required=True, help='Path to PANDA output file')
    parser.add_argument('--lionessfilepath', required=True, help='Path to LIONESS output file')
    parser.add_argument('--output', required=True, help='Output YAML file name')
    parser.add_argument('--input_format', required=False, default='csv', help='Input file format (csv or tsv)')

    args = parser.parse_args()

    generate_sisana_params(
        exp_file=args.exp,
        motif_file=args.motif,
        ppi_file=args.ppi,
        number=args.number,
        outdir=args.outdir,
        method=args.method,
        pandafilepath=args.pandafilepath,
        lionessfilepath=args.lionessfilepath,
        output_filename=args.output,
        input_format=args.input_format
    )