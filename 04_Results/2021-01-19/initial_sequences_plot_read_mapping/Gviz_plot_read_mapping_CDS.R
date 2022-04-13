library(Gviz)
options(ucscChromosomeNames=FALSE)
library(GenomicRanges)



## Coords to plot - will plot one page per entry in list
## Bed format (seq_id <tab> start <tab> end <tab> name)
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=2) {
  stop("Supply list of seq coords to plot + plot name.n", call.=FALSE)
}
coords.list.file <- args[1]
plot.file <- args[2]

#coords.list.file <- "test_list.bed"
#plot.file <- paste(coords.list.file, ".read_mapping_plots.pdf", sep='')

Fasta.file <- "Paulinella_micropora_KR01_nuclear.augustus_190918v4.cds.fa"
Bam.file <- "All_combined.local_mapping.coordsorted.bam"

Annot.file.domains <- "Combined.CDS_coords.bed"
Annot.name.domains <- "Protein Features"
#Annot.file.predGenes <- "Paulinella_micropora_KR01_nuclear.augustus_190918v4.gff3"
#Annot.name.predGenes <- "Predicted Genes"


###################
###################
## Coords 2 Plot ##
###################
###################
coords.list <- read.table(coords.list.file, header = F, sep = '\t',
			col.names = c("seqid", "start", "end", "name"),
			stringsAsFactors=FALSE, comment.char = '#')



#####################
#####################
## Sequence tracks ##
#####################
#####################
track.list <- list()

## Genome track - Dont have to add an actual genome, just an object used for plotting.
Fasta.gTrack <- GenomeAxisTrack()
track.list <- append(track.list, Fasta.gTrack)

## Sequence track - Load fasta sequence so we can display it on plot
Fasta.sTrack <- SequenceTrack(Fasta.file, stream = FALSE)
track.list <- append(track.list, Fasta.sTrack)



############################
############################
## GFF3 Annotation tracks ##
############################
############################
#
# GFF3 format:
# .. mRNA .. ID=gene1.mRNA;Name=gene1.mRNA
# .. exon .. ID=gene1.mRNA.CDS;Parent=gene1.mRNA

## Annotation track - load predicted genes.
#Annot.aTrack.predGenes <- GeneRegionTrack(range = Annot.file.predGenes,
#                                     name=Annot.name.predGenes, transcriptAnnotation = 'transcript',
#                                     size=10,
#                                     showFeatureId=TRUE, showId=TRUE,
#                                     fontcolor.item="black", fontcolor.group="black")
#track.list <- append(track.list, Annot.aTrack.predGenes)



###########################
###########################
## BED Annotation tracks ##
###########################
###########################
#
# Bed format:
# seq_id <tab> start <tab> end <tab> feature_name

## Annotation track - load X.
Annot.aTrack.domains <- AnnotationTrack(range = Annot.file.domains,
                                 name=Annot.name.domains, groupAnnotation="group",
                                 cex=0.5, cex.group=0.5, size=5,
                                 showFeatureId=TRUE, showId=FALSE,
                                 fontcolor.item="black", fontcolor.group="black")
track.list <- append(track.list, Annot.aTrack.domains)



##################
##################
## Align tracks ##
##################
##################

## Alignment track - Bam file
Bam.alnTrack <- AlignmentsTrack(range=Bam.file, name="Reads", isPaired = TRUE, 
                                      col.mates="purple", col.gap="orange", size=0.5)
track.list <- append(track.list, Bam.alnTrack)



##########
##########
## PLOT ##
##########
##########
corrds_expantion_factor <- 0

## Setup 3x4 grid. Inner grid for plots and outer grids are margins
grip.template <- grid.layout(4, 3, widths = unit(c(0.5, 1, 0.5), c("cm", "null", "cm")), 
                             heights = unit(c(0.5, 1, 21, 0.5), c("cm", "cm", "null", "cm")))
#grid.show.layout(grip.template)

pdf(plot.file, width = 18, height = 12, onefile = TRUE)
for (idx in 1:length(coords.list$seqid)) {
  seq.seqid <- coords.list$seqid[idx]
  seq.start <- coords.list$start[idx]
  seq.end <- coords.list$end[idx]
  seq.name <- coords.list$name[idx]
  print (paste("[", Sys.time(), "] PLOT: ", seq.name, "   From: ", seq.start, "   To: ", seq.end, sep=''))
  
  #### Setup grid
  grid.newpage()
  pushViewport(viewport(layout = grip.template))
  
  ####
  #### START HEADER
  #### 
  pushViewport(viewport(layout.pos.col = 2, layout.pos.row = 2))
  grid.text(paste("Plot of gene ",seq.name,sep=''),
            gp = gpar(fontface = "bold", fontsize = 20), just = c("bottom"))
  popViewport()
  #### 
  #### END HEADER
  ####
  
  ####
  #### START TRACKS - plot fasta + bam + annotation tracks
  #### 
  pushViewport(viewport(layout.pos.col = 2, layout.pos.row = 3))
  plotTracks(track.list,
             chromosome=seq.seqid, 
             from=seq.start, extend.left=((seq.end - seq.start)*corrds_expantion_factor), 
             to=seq.end, extend.right=((seq.end - seq.start)*corrds_expantion_factor), add = TRUE)
  popViewport()
  #### 
  #### END TRACKS
  ####
}
dev.off()

print ("Script Finished!")



