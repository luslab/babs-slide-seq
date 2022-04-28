
# Output

The pipeline outputs are in the `results` directory.
For example, here we have two samples (`sample1` and `sample2`):

```
$ ls results 
files  qc  sample1.csv  sample1_dge  sample1.pdf  sample2.csv  sample2_dge  sample2.pdf
```

The `files` and `qc` directories contain intermediate files and raw QC metrics respectively.
So, they can both be deleted.

You are interested in the following files:

 * `sample1.pdf`: the QC plots
 * `sample1_dge`: the digital expression matrix
 * `sample1.csv`: the beads coordinates

## Digital matrix expression and bead coordinates

```
$ tree results/sample1.csv results/sample1_dge
results/sample1.csv [error opening dir]
results/sample1_dge
├── barcodes.tsv.gz
├── features.tsv.gz
└── matrix.mtx.gz
```

## QC metrics

![Page 1](example_output/pages/page-01.png)
![Page 2](example_output/pages/page-02.png)
![Page 3](example_output/pages/page-03.png)
![Page 4](example_output/pages/page-04.png)
![Page 5](example_output/pages/page-05.png)
![Page 6](example_output/pages/page-06.png)
![Page 7](example_output/pages/page-07.png)
![Page 8](example_output/pages/page-08.png)
![Page 9](example_output/pages/page-09.png)
![Page 10](example_output/pages/page-11.png)
![Page 10](example_output/pages/page-12.png)
![Page 10](example_output/pages/page-13.png)
![Page 10](example_output/pages/page-14.png)
![Page 10](example_output/pages/page-15.png)
![Page 10](example_output/pages/page-16.png)
![Page 10](example_output/pages/page-17.png)

