
from shutil import copy

files = ["config.yaml",
         "Snakefile",
         "scripts/solve_network.py",
         "scripts/prepare_sector_network.py"]

for f in files:
    copy(f,snakemake.config['summary_dir'] + '/' + snakemake.config['run'] + '/configs/')
