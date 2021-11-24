
from shutil import copy

files = {
    "config.yaml": "config.yaml",
    "Snakefile": "Snakefile",
    "scripts/solve_network.py": "solve_network.py",
    "scripts/prepare_sector_network.py": "prepare_sector_network.py",
    "../pypsa-eur/config.yaml": "config.pypsaeur.yaml"
}

if __name__ == '__main__':
    if 'snakemake' not in globals():
        from helper import mock_snakemake
        snakemake = mock_snakemake('copy_config')

    for f, name in files.items():
        copy(f,snakemake.config['summary_dir'] + '/' + snakemake.config['run'] + '/configs/' + name)
