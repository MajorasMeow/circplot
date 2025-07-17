# Circplot for R
Function for a customizable circular plot of gene expression, gene correlation and annotation data using ggplot2.

The plot axis is calculated by 360° but begins with a gap for visibility and labels.

It takes gene correlation data and imputed scRNA expression data from a seurat file as input. For both preparatin steps from a seurat object see https://github.com/MajorasMeow/scRNA/blob/main/scfunctions.

It has the following arguments:
```r
cor_data                # your gene correlation data output from linreg_genes() 
expr_data               # your expression csv with genes as columns. remember to clear all non-numerical data before.   
mean_cutoff = 0.1       # cut-off for low expressed gened set to 0.1
top_GO_BIO_terms = 5    # TOP Biological function GO terms that are plotted      
GOI = c(),              # vector of genes of interest to be labeled
species = "human"       # species. can be human, fly 
save = FALSE,           # set TRUE to save file as csv in working directory. saves time if you want to work on aes.
read = FALSE            # set TRUE to read the csv file from working directory
savename = "corplot_df.csv" # custom save name for the file if saved 
alpha_col = FALSE       # sets alpha for correlation bars to -log10 pval
pos_color = "#ff084a"   # color for positive correlation
neg_color = "#005582"   # color for negative correlation
viridis_color = "mako"  # choose viridis color for gene expression ring and for pval of GO Terms
```

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
Gene Actbeta in a population of C4da neurons from Drosophila peripheral nervous system.
Genes linked to neurodegenerative disorders are indicated via gene ontology analysis.
```r
# install dependencies
install.packages(c(    
  "broom",  
  "purrr",     
  "furrr",     
  "progressr"      
))

# apply MAGIC imputation
seurat_object <- FullMagic(seurat_object)

# extract imputation data from seurat object
data <- ImpData(seurat_object)

# perform linear regression analysis on a target gene. this code performs linreg and non-lin reg and choses the better performant regression. calcluates and uses z-scores
target_gene <- "Actbeta"
cor_data <- linreg_genes(data, target_gene, adj_r2_threshold = 0, workers = 1) # set workers to whichever for parallel processing

# perform circplot calculation and save the file

circ_plot(cor_data =result, expr_data=data, read = F, save = T) 
```
![img](https://github.com/user-attachments/assets/9c8b79ac-f0a0-4e44-b0dc-698d1017e463)

```r
# perform circplot from file
circ_plot(cor_data =result, expr_data=data, read = T, save = F) 
```

change and add rings as you like!
