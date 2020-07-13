#################################################
# Heat map plot on distance matrices
#
# based on script of Chiara Barbieri 
# changed by Alex Huebner 13.09.13
# Edited and reduced by Emil Karpinski - 2018-1211
#################################################

#Install packages
#install.packages("reshape")


#Load required libraries
library(reshape)
library(ggplot2)
library(ape)

#Set working directory

setwd("D:\\MastodonPhylogeography\\FinalAnalysis\\Substitution_Analyses\\HeatMap\\")

#Import the combined Blast+No Blast filtered Fasta
Fasta <- read.dna(file="3StDevFiltered_NoOutgroup_Ordered.muscle.fasta", format = "fasta")
#Fasta <- read.dna(file=file.choose(), format = "fasta")

#Check dimensions of input file. Just to make sure it imported correctly. Should be 35
dim(Fasta)

#Compute the pairwise distance matrix, model proportion of sites that differ, missing data deletion in a pairwise manner, and output as matrix
Distance_Matrix <- dist.dna(Fasta, model="F84", pairwise.deletion = T, as.matrix = T)

#Filtering Distance Matrices to just lower triangles. Also removes diagonal.
Distance_Matrix[upper.tri(Distance_Matrix, diag = T)] <- NA

##Not sure what these next two steps actually do, but without then you get a square not a triangle. 
# prepare for plotting with ggplot2
# convert from matrix to data.frame; necessary for ggplot
Dist_DataFrame <- melt(Distance_Matrix)

# reorder the levels according to the original order
Dist_DataFrame$X1 <- factor(Dist_DataFrame$X1, levels = rownames(Distance_Matrix))
Dist_DataFrame$X2 <- factor(Dist_DataFrame$X2, levels = rownames(Distance_Matrix))

#Generates the heatmap distance. Identical to the Alex's original command. 
heatmap <- ggplot(Dist_DataFrame, aes(X1,X2)) + geom_tile(aes(fill = value ), colour="white") + scale_fill_gradient(low= "blue", high = "orange", name="", na.value="white") + labs(x="", y="") + theme(axis.text.x = element_text(angle = 90, hjust = 1, size=7), axis.text.y = element_text(size=7))

#Plots the heatmap
heatmap


#Outputs the heatmap to pdf
pdf("AllMastodons_2Colour_F84.pdf")
heatmap
dev.off()
