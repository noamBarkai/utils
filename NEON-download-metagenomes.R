# taken with modifications from: https://rpubs.com/zoeywerbin/download_metagenomes

# install from CRAN if necessary
# install.packages('neonUtilities')
# install.packages('future.apply')
# install.packages('progressr')


# view available sites at: https://www.neonscience.org/field-sites/field-sites-map

library(neonUtilities) # to download NEON files 
library(future.apply) # to run concurrently
library(progressr) # to update progress bar

download.dir <- file.path(getwd(), "raw_metagenomes") # change to some other location if needed
dir.create(file.path(download.dir))

# dpID's taken from https://data.neonscience.org/data-products/explore 

# DP1.10107.001   Soil microbe metagenome sequences
# DP1.20279.001	  Benthic microbe metagenome sequences
# DP1.20281.001	  Surface water microbe metagenome sequences

# sites taken from: https://www.neonscience.org/field-sites/explore-field-sites
# e.g. 'ABBY'
# for all sites set site = 'all'

metadata <- loadByProduct(dpID = 'DP1.10107.001',
                          site = 'ABBY', 
                          #startdate = "2016-09", enddate = "2016-09",
                          check.size = FALSE, package = 'expanded')
rawFileInfo <- metadata$mms_rawDataFiles 

cat ('number of files found in NEON: ' ,nrow(rawFileInfo))

# check if some files already exist - in which case they are not downloaded again
already.have <- list.files(download.dir, pattern = ".fastq", recursive = T)
already.have <- gsub("_R[12].fastq|_R[12].fastq.gz", "", already.have)
already.have <- basename(already.have)

cat ('number of files already presant that will not be downloaded', length(already.have))

not.have <- rawFileInfo[!rawFileInfo$dnaSampleID %in% already.have,]
not.have <- not.have[!not.have$internalLabID %in% already.have,]

url_list <- unique(not.have$rawDataFilePath)
cat ('number of files downloading: ', length(url_list))

# increase timeout (default: 60s) to 10m to be able to download large files on slower network
options(timeout = max(60 * 10, getOption("timeout")))

# do the downloads concurrently
cl <- parallel::makeCluster(4)
plan(cluster, workers=cl)
handlers(handler_progress(format="[:bar] :percent :eta :message"))

with_progress({
  p <- progressor(along=1:length(url_list))
  results <- future_lapply(1:length(url_list), FUN=function(i) {
    # use wget instead of R download.file as the latter occasionally crashes / gets stuck
    # need to verify wget is installed locally, or replace below with some other download command-line tool of choice
    # can also put this in async execution queue and wait for all to finish
    
    # download.file(url_list[i],  destfile = paste0(download.dir, basename(url_list[i])))
    cmd <- paste('wget', '--no-verbose', '-O', file.path(download.dir, basename(url_list[i])), url_list[i])
    system(cmd)
  }, future.chunk.size=1)
})


had.errors <- FALSE
for (i in 1:length(results)) {
  if (results[i] != 0) {
    cat('failed to download:', url_list[i])
    had.errors <- TRUE
  }
}
if (had.errors) {
  stop('aborting - see errors above')
}

# “untar” each of the downloaded files, and remove the original tar files.
cmd <- paste0('cd ', download.dir, ';', ' for file in *.tar.gz; do tar xzvf "${file}" && rm "${file}"; done')
system(cmd)

# Some of the resulting files are nested within many directories - let’s bring everything to the main directory, and then remove any empty directories.
cmd <- paste0('find ', download.dir, ' -mindepth 2 -type f -exec mv -t ', download.dir, ' -i "{}" +')
system(cmd)
cmd <- paste0('find ', download.dir, '/. -type d -empty -delete')
system(cmd)

# rename the samples according to their NEON sampleIDs, since the current files were named by the sequencing facility. Some files have not been renamed, if they weren’t in the metadata file, since the tar files are bundled samples from multiple sites.
setwd(download.dir)
suffix <- ifelse(grepl("R2_", rawFileInfo$rawDataFileName), "_R2", "_R1")

# Add "BMI" prefix to metagenomes if necessary
old.names <- ifelse(!grepl("BMI", rawFileInfo$internalLabID), paste0("BMI_", rawFileInfo$internalLabID, suffix, ".fastq"), paste0(rawFileInfo$internalLabID, suffix, ".fastq"))
new.names <- paste0(rawFileInfo$dnaSampleID, suffix, ".fastq")
file.rename(from = old.names, to = new.names)

# put the renamed files into their own folder
cleaned <- list.files(download.dir, pattern = "COMP-DNA")
dir.create("clean")
file.rename(cleaned, paste0("clean/", cleaned))

# Remove unpaired files
actual.files <- list.files(file.path(download.dir, "clean"), full.names = T)
basenames <- sapply(strsplit(actual.files, "_R"), `[`, 1)
to.rm <- actual.files[!basenames %in% basenames[duplicated(basenames)]]
file.remove(to.rm)
cat(paste(length(to.rm), "files do not have any read counterpart - omitting"))

# gzip all remaining files
cores <- parallel::detectCores()
clean.dir <- file.path(download.dir, "clean")
files.to.zip <- list.files(clean.dir, "*.fastq")
with_progress({
  p <- progressor(along=1:length(files.to.zip))
  results <- future_lapply(1:length(files.to.zip), FUN=function(f) {
    cmd <- paste0("gzip ", file.path(clean.dir, files.to.zip[f]))
    system(cmd)
  }, future.chunk.size=1)
})


# Get list of downloaded files
actual.files <- list.files(file.path(download.dir), full.names = T, pattern = ".fastq",recursive = F)

# optionally filter files, e.g. remove aquatic files
# to.rm <- actual.files[grepl("Aquatic", actual.files)] 
# file.remove(to.rm)
# actual.files <- actual.files[!actual.files %in% to.rm]

# select only files that have both forward and reverse reads
basenames <- sapply(strsplit(actual.files, "_R"), `[`, 1)
basenames.paired <- unique(basenames[basenames %in% basenames[duplicated(basenames)]])
paired <- actual.files[basenames %in% basenames.paired]
cmd.list <- list()
basenames.paired <- basenames.paired[c(86:87)]

# script to fix unpaired sequences using fastq_pair - can be added from https://github.com/linsalrob/fastq-pair
# fastq.pair.path <- "/location-of-fastq-pair-installation"
# for (i in 1:length(basenames.paired)){
#   cmd <- paste0("cd ", fastq.pair.path, " && ", "fastq_pair ", basenames.paired[[i]], "_R1.fastq ", basenames.paired[[i]], "_R2.fastq")
#   system(cmd)
# }
# problem.names <- paste0(paired, ".paired.fq")
# file.rename(from = problem.names, to = paired)

# clean up
cmd <- paste0("rm -f *single.fq") # removes any .single.fq
#cmd <- paste0("rm -f *paired.fq") # removes any .paired.fq
system(cmd)
