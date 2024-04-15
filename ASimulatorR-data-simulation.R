# Attempting to generate simulated RNA-seq data
setwd("C:/Users/jonan/Documents/1Work/Envisagenics/AS-Evaluation")

library(ASimulatoR)


GRCh38_gtf <- "../Data/GeneSet/Homo_sapiens.GRCh38.111.gtf"


# Example call to runASimulatoR

#get_file <- system.file('exdata', GRCh38_gtf, package = "ASimulatorR")

#exon_superset <- get_exon_supersets(get_file)

exon_superset <- get_exon_supersets(GRCh38_gtf)
save(exon_superset, file = "../Data/GeneSet/exon_superset.RData")
load("../Data/GeneSet/exon_superset.RData")

# Define your input_dir with the path to your FASTA data
input_dir <- "../Data/GeneSet/fasta/"


# define, how many groups and samples per group you analyze. 
# Here we create a small experiment with two groups with one 
# sample per group:
num_reps = c(1,1)

outdir <- "../Data/ASimulatoR_output"

# define the number of genes you want to work with. 
# If you want all exons, do not specify this parameter 
# or set it to NULL
# here we create splice variants from 9 exon supersets:
max_genes = 9



# in this example we use relative frequencies
# here we produce eight variants with one of each AS events as well as one variant containing every event
# if probs_as_freq was FALSE, a random number would be drawn for each event-superset combination and only if it was smaller than 1/9 the AS event would be created
probs_as_freq = T
event_freq = 
  setNames(rep(1/9, 9),
           c('es', 'mes', 'ir', 'a3', 'a5', 'afe', 'ale', 'mee', 'es,ir,mes,a3,a5,afe,ale,mee'))



# we use the previously created superset to simulate splice variants from, since it is saved in the same directory as the gtf
# if no superset is found, a new one will be created
simulate_alternative_splicing(input_dir = input_dir,
                              outdir = outdir, 
                              event_probs = event_freq,
                              probs_as_freq = probs_as_freq, 
                              max_genes = max_genes,
                              num_reps = num_reps,
                              verbose = TRUE)

install.packages("ggbio")
library(ggbio)

BiocManager::install("ggbio")

gtf = rtracklayer::import("../Data/ASimulatoR_output/splicing_variants.gtf")


# the gene id of the variant with all events
gene_id = gtf$gene_id[grep('es,ir,mes,a3,a5,afe,ale,mee', gtf$transcript_id, fixed = T)[1]]
exons = gtf[gtf$type == 'exon' & gtf$gene_id == gene_id]
suppressWarnings(ggbio::autoplot(split(exons, exons$transcript_id)))
exons$transcript_id

head(exons)



library(Gviz)


# Assuming 'genes' is a GRanges object containing gene regions


# Define a genome axis track
genomeAxis <- GenomeAxisTrack()

# Create an annotation track for exons
# annotationTrack <- AnnotationTrack(range = exons, name = "Exons", stack = TRUE, fill = "darkblue")
annotationTrack <- AnnotationTrack(range = exons, name = "Exons", stacking = "pack", fill = "purple")

# Plot tracks
plotTracks(list(genomeAxis, annotationTrack), from = min(start(exons)), to = max(end(exons)))




####################
geneModelTrack <- GeneRegionTrack(genes, name="Gene Model")














