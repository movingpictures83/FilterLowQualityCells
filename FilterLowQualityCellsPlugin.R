# This is demonstration of using SINCERA to analysis human IPF and normal epithelial single cell RNA-seq data (Xu et al., JCI Insight 2016)
# Author: Minzhe Guo (minzhe.guo@cchmc.org)

library(SINCERA)

dyn.load(paste("RPluMA", .Platform$dynlib.ext, sep=""))
source("RPluMA.R")

input <- function(inputfile) {
  parameters <<- read.table(inputfile, as.is=T);
  rownames(parameters) <<- parameters[,1];
    pfix = prefix()
  if (length(pfix) != 0) {
     pfix <<- paste(pfix, "/", sep="")
  }
}

run <- function() {}

output <- function(outputfile) {

exprmatrix <- read.csv(paste(pfix, parameters["expression", 2], sep="/"))
dataid <- readLines(paste(pfix, parameters["data", 2], sep="/"))
diagnosis <- readLines(paste(pfix, parameters["diagnosis", 2], sep="/"))
# The analysis starts with running the construct function to create an R S4 object, which will hold all the data and analysis results.
# The function takes expression matrix and cell sample information as input
sc <- construct(exprmatrix=exprmatrix, samplevector=dataid)

# After contruction, you can use setCellMeta to add more cell information to sincera
# such as adding the condition information (CONTROL or IPF) of cells into sincera object
sc <- setCellMeta(sc, name="CONDITION", value=diagnosis)

# use getCellMeta functoin to assess a specific meta data
# In most of the SINCERA functions, cell grouping will be based on the GROUP meta data
# The GROUP meta was initialized to sample information during the construction
#table(getCellMeta(sc, name="GROUP"))

# Identify and remove low quality cells.
# The key parameters of running this function include: “min.expression”,
# which specifies the minimum expression value for a gene to be considered as expressed,
# and “min.genes”, which specifies the lower bound of the number of expressed genes in a cell.
sc <- filterLowQualityCells(sc, min.expression=1, min.genes=1000, do.plot = T)
saveRDS(sc, outputfile)
}
