# Circplot for R: A function for visualizing gene correlations in single-cell data with disease annotation
Function for a customizable publication-ready circular plot for integrated single-cell data analysis, gene correlation and GO Term annotation using ggplot2.

This tool is designed for visualizing the relationship between a gene of interest, its top correlated genes, their expression levels, and their functional annotations (GO terms) all in a single, intuitive plot.
It function focuses by default on neurodegenerative disease annotation (from Jensen_DISEASES_Curated_2025 and Rare_Diseases_GeneRIF_Gene_Lists (human), or Human_Disease_from_FlyBase_2017 (Fly)) and Top 5 Biological Function GO Terms for each bracket of 0.1 R_squared increments.
The plot axis is calculated by 360° but begins with a gap for visibility and labels.

It takes gene correlation data and imputed scRNA expression data from a seurat file as input. For both preparatin steps from a seurat object see https://github.com/MajorasMeow/scRNA/blob/main/scfunctions.
![blacksmall](https://github.com/user-attachments/assets/6c8fe91e-209d-4198-9f4b-4644501f4a86)


## install dependencies for circ_plot()
```r
# install dependencies
install.packages(c(
  "tidyverse",   
  "ggforce",     
  "ggnewscale",  
  "scales",     
  "enrichR",     
  "ggrepel"      
))
```
### install magic for imputation https://github.com/cran/Rmagic
# EXAMPLE

We take the pbmc3k dataset from the SeuratData Package and want to find out which are high expressing genes that correlate with the gene CD68 in the whole dataset.

### load and install the dataset
```r
devtools::install_github('satijalab/seurat-data')
library(SeuratData)
InstallData("pbmc3k")

```
Before we can calculate gene correlations we have to install some dependencies and prepare the data.
### install dependencies for regression function
```r

install.packages(c(    
  "broom",  
  "purrr",     
  "furrr",     
  "progressr"      
))
````
... and use Markov Affinity-based Graph Imputation of Cells (MAGIC) for imputation of expression values.
### apply MAGIC imputation on your seurat object
```r
pbmc3k <- FullMagic(pbmc3k)
```
### extract imputation data from seurat object into a data frame
```r
data <- ImpData(pbmc3k)
```
### perform linear regression analysis on a target gene. this code performs linreg and non-lin reg and chooses the better performant regression. calcluates and uses z-scores before regression.
We chose our target gene and perform the regression analysis.
```r
target_gene <- "CD68"
cor_data <- linreg_genes(data, target_gene, adj_r2_threshold = 0, workers = 1) # adjust adj_r2 threshold if you aim to exclude any low correlation. set workers to whichever for parallel processing
```
## 1. Plot your data 

First, load the circ_plot() function into your R session.
It has the following arguments:

```r
cor_data                # your gene correlation data output from linreg_genes() 
expr_data               # your expression csv with genes as columns. remember to clear all non-numerical data before.   
mean_cutoff = 0.1       # cut-off for low expressed gened set to 0.1
top_GO_BIO_terms = 5    # TOP Biological function GO terms that are plotted      
GOI = c()               # vector of genes of interest to be labeled
species = "human"       # species. can be human, fly 
save = FALSE            # set TRUE to save file as csv in working directory. saves time if you want to work on aes.
read = FALSE            # set TRUE to read the csv file from working directory
savename = "corplot_df.csv" # custom save name for the file if saved 
alpha_col = FALSE       # sets alpha for correlation bars to -log10 pval
pos_color = "#ff084a"   # color for positive correlation
neg_color = "#005582"   # color for negative correlation
viridis_color = "mako"  # choose viridis color for gene expression ring and for pval of GO Terms
theme_color= "black"    # Can be black or white
```

Sometimes we use non numerical or indexing columns etc for certain plots. Here we dont want them.
Make sure no non numerical columns are in your dataset! Also indexing columns, trajectory values, or other non-expressing numerical columns should be excluded before.

```r
# example of non-numerical and index columns we want to exclude
data<-select(data, -c("ident",       # cluster name from AssignClusterLabels
                      "cellnames",   # cellnames name from AssignClusterLabels
                      "X"))          # indexing column
```

We want to plot our first with save= TRUE as this function prepares our data and uses enrichR for GO database queries.
Once our function is done, we can set save = FALSE and read = TRUE to read directly from the saved csv.

Also let's take a mean threshold of 0.5 so we just get the top expressed genes.
The Plot has a default labeling for neurodegenererative disorders, but if no hits are available it labels all disease related genes found in The Jennsen and rare disease GO DB.

# perform circplot calculation and save the file
```r
circ_plot(cor_data =result, expr_data=data, read = F, save = T, mean_cutoff = 0.5)

# perform circplot from file
circ_plot(cor_data =result, expr_data=data, read = T, save = F) 
```
Great! now we have our circular plot!
Think as the x-axis is bent circular (0-360) and the y-axis is radial (Here we use low values as it is set to slope intensity). In theory you can add more layers and set y-values in their geom to n which defines the level of the object ring/label/etc. For distribution like the GO Biol. Function Terms, each Term and its value were mapped to a certain x-value across 360°. It shows the TOP 5 GO terms for each bracket of adj_r2.
If we set our threshold to 0 we have 10 brackets, 9 for 0.1, 8 for 0.2...
We can plot in a light "white" theme...
![whitesmall](https://github.com/user-attachments/assets/88713553-c377-4655-bc61-c0ef8dac2931)

...or as a black theme with ivory labels.
```r
circ_plot(cor_data =result, expr_data=data, read = T, save = F, theme_color = "black") 
```

![blacksmall](https://github.com/user-attachments/assets/ca1a6c2f-0bb7-423f-84b7-881fbfb68c4e)

Change and add rings as you like!
