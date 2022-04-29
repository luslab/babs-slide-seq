
# Running the pipeline

If you're from the Crick, just `ssh` to CAMP, then run

```
git clone https://github.com/bahnk/SlideSeq
cd SlideSeq
```

Configure your `params.yml` file and `design.csv` file and then run:

```
module load Nextflow/20.12.0-edge Singularity/3.6.4
export NXF_SINGULARITY_CACHEDIR=$WORK/singularity
nextflow run main.nf -params-file params.yml
```

