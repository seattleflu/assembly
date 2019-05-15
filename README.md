# Seattle Flu Assembly pipeline

## Installation:
First install Conda, then run:
```
conda env create -f envs/seattle-flu-environment.yaml
```

Running (this will be updated):
```
snakemake --configfile config/config.json -k
```
OR (if on Rhino)
```
snakemake -w 60 --configfile config/config-all-refs.json --cluster-config config/cluster.json --cluster "sbatch --nodes=1 --tasks=1 --mem={cluster.memory} --cpus-per-task={cluster.cores} --tmp={cluster.disk} --time={cluster.time} -o all_output.out" -j 20 -k
```

You can use `--use-conda` if there are rule-specific environment files.
