
# Specify the path to your ichorCNA installation
ichorCNA_path <- "/sci/labs/aviadz/nirdagan/projects/R-3.6.0/library"

# Load each function individually
library("ichorCNA", lib.loc = ichorCNA_path)

# Define arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if there is exactly one argument
if (length(args) != 2) {
  stop("Usage: Rscript your_script.R <your_argument>")
}

sample_id <- args[1]
sample_path <- args[2]

# # Check if the directory already exists
# if (!file.exists(sample_id)) {
#   # Create the directory
#   dir.create(sample_id)
# }

# Run runIchorCNA.R script with specified parameters
system(paste("Rscript", "/sci/labs/aviadz/nirdagan/projects/ichorCNA/scripts/runIchorCNA.R",
              "--id", sample_id,
              "--WIG", sample_path, "--ploidy", "'c(2,3)'",
              "--normal 'c(0.5,0.6,0.7,0.8,0.9)'", "--maxCN 5",
              "--gcWig", "/sci/home/shario/icore-home/R/x86_64-pc-linux-gnu-library/4.0/ichorCNA/extdata/gc_hg19_1000kb.wig", "--mapWig", "/sci/home/shario/icore-home/R/x86_64-pc-linux-gnu-library/4.0/ichorCNA/extdata/map_hg19_1000kb.wig",
              "--centromere", "/sci/home/shario/icore-home/R/x86_64-pc-linux-gnu-library/4.0/ichorCNA/extdata/GRCh37.p13_centromere_UCSC-gapTable.txt",
              "--normalPanel", "/sci/labs/aviadz/nirdagan/projects/ichorCNA/inst/extdata/HD_ULP_PoN_1Mb_median_normAutosome_mapScoreFiltered_median.rds",
              "--includeHOMD", "False", "--chrs 'c(1:22)'",
              "--chrTrain 'c(1:22)'", "--estimateNormal", "True",
              "--estimatePloidy", "True", "--estimateScPrevalence", "True",
              "--scStates 'c(1,3)'", "--txnE", "0.9999",
              "--txnStrength", "10000", "--outDir", sample_id), intern = TRUE)

# Done!

# Done!
message("Successfully ran IchorCNA using runIchorCNA.R script and generated report!")

