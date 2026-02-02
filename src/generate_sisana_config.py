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

# snakemake provides inputs/params/outputs
exp_file = snakemake.input.exp
output_filename = snakemake.output.config_file
motif_file = snakemake.params.motif
ppi_file = snakemake.params.ppi
number = snakemake.params.number
outdir = snakemake.params.outdir
method = snakemake.params.method
pandafilepath = snakemake.params.pandafilepath
lionessfilepath = snakemake.params.lionessfilepath
input_format = snakemake.params.input_format

generate_sisana_params(
        exp_file=exp_file,
        motif_file=motif_file,
        ppi_file=ppi_file,
        number=number,
        outdir=outdir,
        method=method,
        pandafilepath=pandafilepath,
        lionessfilepath=lionessfilepath,
        output_filename=output_filename,
        input_format=input_format
    )