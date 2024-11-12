
library(Seurat)

library(dplyr)
library(patchwork)

library(Signac)

require(magrittr)
require(readr)
require(Matrix)
require(tidyr)
library(DropletUtils)


# peak-bc matrix
mex_dir_path <- 'C:\\Users\\ssrikrishnan6\\scATAC\\scATAC_GEO_Dataset\\CLL_BMMC_ATAC'

mtx_path <- paste(mex_dir_path, "matrix.mtx", sep = '\\')
feature_path <- paste(mex_dir_path, "peaks.bed", sep = '\\')
barcode_path <- paste(mex_dir_path, "barcodes.tsv", sep = '\\')

features <- readr::read_tsv(feature_path, col_names = F) %>% tidyr::unite(feature)

barcodes <- readr::read_tsv(barcode_path, col_names = F) %>% tidyr::unite(barcode)

mtx <- Matrix::readMM(mtx_path) %>%
  magrittr::set_rownames(features$feature) %>%
  magrittr::set_colnames(barcodes$barcode)

## Converting mtx file to h5 files:

my.counts <- readMM(mtx_path)
my.counts <- as(my.counts, "dgCMatrix")
cell.ids <- paste0("BARCODE-", seq_len(ncol(my.counts)))

ngenes <- nrow(my.counts)
gene.ids <- paste0("ENSG0000", seq_len(ngenes))
gene.symb <- paste0("GENE", seq_len(ngenes))

# Writing this to file:
tmpdir <- 'C:\\Users\\ssrikrishnan6\\scATAC\\scATAC_GEO_Dataset\\CLL_BMMC_ATAC\\h5_converted'
write10xCounts(tmpdir, my.counts, gene.id=gene.ids, 
               gene.symbol=gene.symb, barcodes=cell.ids)
list.files(tmpdir)


# Creating a version 3 HDF5 file:
# tmph5 <- tempfile(tmpdir = "~/C:/Users/ssrikrishnan6/tmp", fileext=".h5")
write10xCounts(tmpdir, type = "HDF5", my.counts, gene.id=gene.ids, 
               gene.symbol=gene.symb, barcodes=cell.ids, version='2')


### Creating framents.tsv file with mtx:

#Create dataframe of 5 columns
#Columns names -> chrom, chromStart, chromEnd, barcode, readSupport

for 


##### Reading converted H5 file
counts <- Read10X_h5(filename = "C:\\Users\\ssrikrishnan6\\scATAC\\scATAC_GEO_Dataset\\CLL_BMMC_ATAC\\h5_converted.h5")

# features <- readr::read_tsv(feature_path, col_names = F) %>% tidyr::unite(feature)


temp_matrix <- readMM(mtx_path)
tempmatrix<-as.data.frame(temp_matrix)
tempmatrix

chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  genome = 'hg19',
  fragments = '../vignette_data/atac_v1_pbmc_10k_fragments.tsv.gz',
  min.cells = 10,
  min.features = 200
)

bmmc <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = metadata
)


bmmc



