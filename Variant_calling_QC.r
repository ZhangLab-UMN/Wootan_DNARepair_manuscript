## filter_strelka2.R

# load in libraries
library(vcfR)
library(seqinr)
library(polymorphology)
library(tidyverse)

# Load in genome 
genome<- read.fasta("/data/varient_calling/10.fasta")
names(genome)<- gsub("At", "", names(genome))

# Load and process the GFF file
genes <- fread("/home/zhangumn/woota001/RNA_seq_alignments/genomes/Athaliana/50_annotation/10.gff")
#genes <- genes[V3 == "gene"]
genes$CHROM <- as.character(gsub("At", "", genes$V1))  # Assuming chromosome names like At1, At2, etc.

genes$ID <- 1:nrow(genes)
genes$START <- genes$V4
genes$STOP <- genes$V5
setkey(genes, CHROM, START, STOP)

bp<- c("A","C","G","T")

# List VCF files excluding .tbi files. These are the directory outputs of the strelka pipeline
files <- list.files("/data/all_comparisons_strelka2/", full.names = FALSE)
normals <- unique(sapply(files, function(x) paste(unlist(strsplit(x, split = "-"))[1], collapse="_")))
tumor <- unique(sapply(files, function(x) paste(unlist(strsplit(x, split = "-"))[2], collapse="_")))

results_wt<-rbindlist(lapply(normals, function(norm){
  cat(norm)
  mutations<-rbindlist(lapply(tumor, function(f){
    if(norm!=f){ #compares each tumor sample to the current norm sample if the norm and tumor aren't the same
      dat<-read.vcfR(paste0("data/all_comparisons_strelka2/",norm,"-",f,"/results/variants/somatic.snvs.vcf.gz"), verbose = F)
      calls<-cbind(data.table(dat@fix),data.table(dat@gt))
      calls$file<-f
      calls$norm<-norm
      calls$EVS<-as.numeric(gsub(".+SomaticEVS=","",calls$INFO)) #somatic varient score
      calls$POS<-as.numeric(calls$POS)
      calls$start<-calls$POS
      calls$stop<-calls$POS
      calls$CHROM<-as.character(calls$CHROM)
      calls$POSITION<-calls$POS
      return(calls)
      
    }
  }))
  
  mutations$unique<-paste(mutations$CHROM, mutations$POS, sep = "_") # creates a new column with the chromosome and position of variant
  mutations$ID<-1:nrow(mutations) #creates a unique identifer for each mutation
  tab<-data.table(table(mutations$unique))
  mutations$Mut_N<-tab$N[match(mutations$unique, tab$V1)] #counts the occurance of each unique mutation
  setkey(mutations, CHROM, start, stop)
  mutations$genic<-features_overlap_mutation(features = genes, mutations = mutations) #flags mutations overlapping with genes
  return(mutations)
}))

results_wt$norm_geno<-substr(results_wt$norm,1,2)

fwrite(results_wt, "/data/strelka2_results_all.csv")

setkey(results_wt, CHROM, start, stop)
results_wt$unique2 <- paste(results_wt$unique, results_wt$file, sep = "_")
tab_wt <- data.table(table(results_wt$unique2))
results_wt$results_N <- tab_wt$N[match(results_wt$unique2, tab_wt$V1)]
PASS_wt <- results_wt[Mut_N == 1 & CHROM %in% c("At1", "At2", "At3", "At4", "At5") & FILTER == "PASS"] # variants can't be in Atx, must pass strelka's qcs, and be a unique mutation
# Mut_N=1 meaning the unique mutation only occurs once in the dataset

PASS_wt$unique3 <- paste(PASS_wt$unique, PASS_wt$file, sep = "_") 
tab_wt <- data.table(table(PASS_wt$unique3))
PASS_wt$N <- tab_wt$N[match(PASS_wt$unique3, tab_wt$V1)] 
PASS_wt$geno <- substr(PASS_wt$file, 1, 2)

PASS_wt$tumor_ref_depth <- as.numeric(apply(PASS_wt, 1, function(x) unlist(strsplit(unlist(strsplit(x["TUMOR"], split = ":"))[4+which(bp == x["REF"])], split = ","))[1]))
PASS_wt$tumor_alt_depth <- as.numeric(apply(PASS_wt, 1, function(x) unlist(strsplit(unlist(strsplit(x["TUMOR"], split = ":"))[4+which(bp == x["ALT"])], split = ","))[1]))
PASS_wt$tumor_depth <- PASS_wt$tumor_ref_depth + PASS_wt$tumor_alt_depth
PASS_wt$depth_pct <- PASS_wt$tumor_alt_depth / PASS_wt$tumor_depth

PASS2_wt <- PASS_wt[N==11]
# the varient is called in 11 of the comparisons to norm. When N==1, 11 of the files being comared (the norm), also have that varient. When N==11, only one of the samples contains the variant

unique4 <- unique(PASS2_wt[, c("CHROM", "POS", "REF", "ALT", "unique", "unique3", "tumor_alt_depth", "tumor_depth", "genic", "file", "Mut_N", "POSITION", "depth_pct", "N", "results_N", "geno", "EVS"), with = FALSE])
unique4$dup <- duplicated(unique4$unique)

unique4$chi <- apply(unique4, 1, function(x){
  chisq.test(c(as.numeric(x["tumor_alt_depth"]), as.numeric(x["tumor_depth"]) - as.numeric(x["tumor_alt_depth"])))$p.value
})

unique4$type <- ifelse(unique4$depth_pct > 0.5 & unique4$chi < 0.05, "homozygous", "somatic")
unique4$type[unique4$chi > 0.001] <- "heterozygous"
PASS3 <- unique4[dup == FALSE  & tumor_depth>20 & tumor_alt_depth>4] # no duplicates, at least 21 reads per variant, at least 5 of those reads contain the variant

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Biostrings")
library(Biostrings)
genome <- readDNAStringSet("/data/10.fasta")
PASS3$context_SBS <- mapply(function(chrom, pos) {
  seq <- genome[[chrom]]  # Access the sequence for the chromosome
  start_pos <- max(1, pos - 3)  # Adjust to avoid negative positions
  end_pos <- min(length(seq), pos + 3)
  as.character(subseq(seq, start=start_pos, end=end_pos))
}, chrom = PASS3$CHROM, pos = PASS3$POS)

PASS3$context_SBS <- contexts(vars = PASS3, fasta = genome)
PASS3$mut <- paste(substr(PASS3$context, 2, 2), substr(PASS3$context, 6, 6), sep = ">")
PASS3$transversion <- !PASS3$mut %in% c("C>T", "T>C")
PASS3$CtoA <- PASS3$mut %in% c("C>A")



PASS3$ID <- 1:nrow(PASS3)
PASS3$start <- PASS3$POS
PASS3$stop <- PASS3$POS
setkey(PASS3, CHROM, start, stop)
PASS3$IG <- !PASS3$genic  # Intergenic regions
PASS3$loc <- ifelse(PASS3$genic, "genic", "intergenic")
fwrite(PASS3, "/data/WTMut_HS_mutations_all.csv")
