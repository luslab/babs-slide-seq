
# Slide-Seq

This pipeline pre-configured for people working at th [Francis Crick Institue](https://www.crick.ac.uk/?gclid=EAIaIQobChMIodDA66K59wIVF-vtCh3_SwEJEAAYAiAAEgKrkvD_BwE).
If you're not from the Crick, it should be easy to configure the pipeline for your system.
Just modify these two files:

 * `nexflow.config`
 * `conf/process.config`

Here is the documentation:

 1. [Pipeline structure](doc/structure.md)
 2. [Pipeline configuration](doc/config.md)
 3. [Running the pipeline](doc/run.md)
 4. [Pipeline output](doc/output.md)

In a nutshell, if you're from the Crick and want to test it, just `ssh` to CAMP and do:

```
git clone https://github.com/bahnk/SlideSeq
cd SlideSeq
sh run
```
If it fails it's probably because your `MODULEPATH` is missing some locations, and/or because you don't have access to the BABS reference area (`/camp/svc/reference/Genomics/babs`).

Any issues ==> `nourdine.bah@crick.ac.uk`. Cheers,

