# Seattle Flu Assembly pipeline

## Installation:
First install Conda, then run:
```
conda env create -f envs/seattle-flu-environment.yaml
```


Running (this will be updated):
```
snakemake -k
```
OR (if on Rhino)
```
snakemake -w 60 --cluster-config config/cluster.json --cluster "sbatch --nodes=1 --tasks=1 --mem={cluster.memory} --cpus-per-task={cluster.cores} --tmp={cluster.disk} --time={cluster.time}" -j 20 -k
```
