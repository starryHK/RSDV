##### **The RSDV app allows users to visualize differentitally expressed genes(DEG) starting with count data.**

- *Explore* the app's features with the example data set pre-loaded by clicking on the tabs above.
- *Upload* your data manually.

## Instructions

The app is hosted on the website: https://RSDV/ 

Code can be found on github: https://github.com/starryHK/shiny-RSDV/

To run this app locally on your machine, download R or RStudio and run the following commands once to set up the environment:

```
install.packages(c("shiny","ggplot2","gplot2",DESeq2","RColorBrewer","DT")
## Try URL:https://rstudio.github.io/DT 
## You can download the R packages "DT",and learn more from
https://cran.r-project.org/web/packages/DT/DT.pdf

source("https://bioconductor.org/biocLite.R")
biocLite(c("voom","edgeR"))
```

You may now run the shiny app with just one command in R:

```
shiny::runGitHub("RSDV", "starryHK")
```

<a name="inputdata"></a> 

## Input Data

You may use this app by

1. Exploring the pre-loaded example data set.
   
   This is a pre-loaded mouse macrophages RNA-seq example for exploring the app's features.
2. Upload your own data that is Count data (or log2-expression data) which come from transcriptome sequencing.

<a name="dataformat"></a> 

### Data Format

- File must be the row counts,not normalized data,e.g.FPKM,TPKM,TPM.
- File must have a header row.
- First/Left-hand column(s) must be gene identifiers.
- First/Left-hand column(s) must be defined as row_name.

<a name="example"></a>

### Example of Data format

- Each row denotes a gene, each column denotes a sample.

![example3](example3.png)

Analysis: When raw counts are uploaded, the data is then analyzed by the app. The app uses the voom method from the ‘limma’ Bioconductor package to transform the raw counts into logged and normalized intensity values. These values are then analyzed via linear regression where gene intensity is regressed on the group factor. P-values from all pairwise regression tests for group effect are computed and Benjamini-Hochberg false discovery rate adjusted p-values are computed for each pairwise comparison. 

Example file: https://github.com/starryHK/shiny-RSDV/

<a name="degtable"></a>

### DEG Table

- Column A provide gene name.
- Column C and column G provide Fold Changes and FDR,respectively.
- We use both log2FC and FDR to file DEG.

![example4](example4.png)

Analyzed data must contain some kind of expression measure for each sample (i.e. counts, normalized intensities, CPMs), and a set of p-values with corresponding fold changes for those p-values. For instance, if you have a p-value for the comparison of control vs exp , you can upload the observed fold change or log2(fold change) between control vs exp. If you have a more complex design and do not have fold changes readily available, you may upload the test statistics or other similar measures of effect size as placeholders. The fold changes are mainly used in the volcano plots. We recommend uploading p-values that are adjusted for multiple comparisons (such as q-values from the qvalue package, or adjusted p-values from p.adjust() function in R).

Example file: : [https://github.com/shiny-RSDV/](https://github.com/???RSDVapp/)

<a name="vis"></a> 

## Visualizations

### Group Plots

<a name="pcaplots"></a>

#### PCA Plot

This plot uses Principal Component Analysis (PCA) to calculate the principal components of the expression data using data from all genes. Euclidean distances between expression values are used. Samples are projected on the first two principal components (PCs) and the percent variance explained by those PCs are displayed along the x and y axes. Ideally your samples will cluster by group identifier.

![PCA1](PCA1.png)

<a name="analysisplots"></a>

### Analysis Plots

These plots use the p-values and fold changes to visualize your data.

<a name="volcano"></a>

#### Volcano Plot

This is a scatter plot log fold changes vs –log10(p-values) so that genes with the largest fold changes and smallest p-values are shown on the extreme top left and top right of the plot. Hover over points to see which gene is represented by each point.

 (<https://en.wikipedia.org/wiki/Volcano_plot_(statistics)>)

![volcanoplot](volcanoplot.png)

<a name="scatterplots"></a>

#### Scatter Plot

This is a scatter plot of average gene expression in one group against another group. This allows the viewer to observe which genes have the largest differences between two groups. The smallest distances will be along the diagonal line, and points far away from the diagonal show the most differences. Hover over points to see which gene is represented by each point.

<a name="boxplots"></a>

### Gene Expression Boxplot

Use the search bar to look up genes in your data set. For selected gene(s) the stripchart (dotplot) and boxplots of the expression values are presented for each group. You may plot one or multiple genes along side each other. Hover over points for more information about the data.

![boxplot](boxplot.png)

<a name="heatmaps"></a>

### Heatmap

A heatmap of expression values are shown, with genes and samples arranged by unsupervised clustering. You may filter on test results as well as P-value cutoffs. By default the top 100 genes (with lowest P-values) are shown.![example](example.png)

![example5](example5.jpeg)

![example2](example2.png)
