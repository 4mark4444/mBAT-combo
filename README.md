# mBAT-combo Portal — Progress Report


## How to Run

```r
# Install dependencies if needed
install.packages(c("shiny", "plotly", "bslib", "DT", "DBI", "RSQLite"))

# Run from the project directory
shiny::runApp("app.R")
```

## What Has Been Built

The core web application is complete and running locally. It is built in R Shiny and covers the following:

### Home Page
A landing page with a brief description of the portal and buttons to navigate to the Results and Search sections.

### Results Tab
- Interactive Manhattan plots for all three methods: **mBAT-combo**, **mBAT**, and **fastBAT**
- Genes are plotted by chromosomal position with colour-coding by chromosome
- Hovering over any point shows the gene symbol, Ensembl ID, genomic coordinates, number of SNPs, all three p-values, top SNP and its p-value, and number of eigenvalues — everything present in the mBAT-combo output
- A sortable results table of genes is shown below the plots
- A **Download CSV** button allows users to export the results for the currently selected trait

### Search Tab
- Users can search by gene symbol, gene name, or Ensembl ID and filter by all three p values. 
- Users can also search by trait name. 
 **This part uses inverted index search which I think should work even on a larger data base. **

---

## What Is Not Done Yet 
There are several things on the spec that I did not do. 
1. I did not really know what to do with the cell type, since I don't think the original mBAT-combo has cell type related functions. 
2. I am not sure what to plot for the gene-centric view of the graph. I am not sure how should I arrange the traits, since an alphabetical line up does not really make sense to me. 


## Reference
## References

- [Yang Lab UKB imputed fastGWA summary statistics](https://yanglab.westlake.edu.cn/software/gcta/)
- [mBAT hg19 GENCODE v40 gene list](https://github.com/Share-AL-work/mBAT/blob/main/glist_ensgid_hg19_v40.txt)
- [CNCR MAGMA 1000 Genomes Phase 3 EUR Build37 reference data](https://cncr.nl/research/magma/)
