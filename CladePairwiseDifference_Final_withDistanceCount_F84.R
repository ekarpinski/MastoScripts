##Pairwise distance measure for the 

#Load in the ape library
library(ape)

#Set the working directory to the local mastodon phylogeography folder
setwd("D:\\MastodonPhylogeography\\FinalAnalysis\\Substitution_Analyses\\")

#Import the combined Blast+No Blast filtered Fasta
Fasta <- read.dna(file="3StDevFiltered_NoOutgroup.muscle.fasta", format = "fasta")

#Check dimensions of input file. Just to make sure it imported correctly. Should be 35
dim(Fasta)

#Compute the pairwise distance matrix, model proportion of sites that differ, missing data deletion in a pairwise manner, and output as matrix
PW_DistMatrix <- dist.dna(Fasta, model="F84", pairwise.deletion = T, as.matrix = T)

#Compute the pairwise distance matrix, model N for raw counts, deletion in a pairwise manner, and output as matrix
Comb_PW_Sub_Matrix <- dist.dna(Fasta, model="N", pairwise.deletion = T, as.matrix = T)

#Compute the indel matrix, only counts indels, deletion in pairwise manner, output as matrix.
Comb_PW_Indel_Matrix <- dist.dna(Fasta, model="indel", pairwise.deletion = T, as.matrix = T)

#Create a combined distance matrix adding the Substituion matrix (Comb_PW_Sub_Matrix) and indel matrix (Comb_PW_Indel_Matrix) for a "total distance" (Note: ambigous sites not counted)
Comb_Final_Matrix <- Comb_PW_Sub_Matrix + Comb_PW_Indel_Matrix

#Creating clade vectors to filter the final matrices
#Pulls the rownames (identical to column names) from the distance matrix
All_OTUs <- rownames(Comb_Final_Matrix)

#Creates a list containing sample names, and clade designations based on number. 1 - Large Alaska Yukon; 2 - Great Lakes; 3 - Virginia; 4 - Alberta/Missouri; 5 - Small Alaska; 6 - RAM P94.16.1B; 9 - Everything Else
MasterClades <- list(sample=All_OTUs, clade = c(4,2,1,2,1,2,3,9,6,4,9,1,1,2,1,9,1,2,2,5,1,2,2,1,5,3,1,1,1,1,2,1,2,1,1))

#Creates Clade vectors from the above.
Lrg_AlaskaYukon <- MasterClades$sample[which(MasterClades$clade==1)]
Lrg_AlaskaYukon_AndRam <- MasterClades$sample[which(MasterClades$clade==1|MasterClades$clade==6)]
GreatLakes <- MasterClades$sample[which(MasterClades$clade==2)]
Virginia <- MasterClades$sample[which(MasterClades$clade==3)]
GLV<- MasterClades$sample[which(MasterClades$clade==2|MasterClades$clade==3)]
AlbertaMissouri<- MasterClades$sample[which(MasterClades$clade==4)]
SmallAlaska<- MasterClades$sample[which(MasterClades$clade==5)]


##
#Filtering the Pairwise Matrices above to just the clades listed above.
##

#Large Alaska Yukon Clade (no RAM P94.16.16B)
Lrg_AlaskaYukon_DistMatrix <- Comb_PW_Sub_Matrix[Lrg_AlaskaYukon,Lrg_AlaskaYukon]
Lrg_AlaskaYukon_IndelMatrix <- Comb_PW_Indel_Matrix[Lrg_AlaskaYukon,Lrg_AlaskaYukon]
Lrg_AlaskaYukon_FinalMatrix <- Comb_Final_Matrix[Lrg_AlaskaYukon,Lrg_AlaskaYukon]
Lrg_AlaskaYukon_PW_DistMatrix <- PW_DistMatrix[Lrg_AlaskaYukon,Lrg_AlaskaYukon]

#Large Alaska Yukon Clade (with RAM P94.16.16B)
Lrg_AlaskaYukon_AndRam_DistMatrix <- Comb_PW_Sub_Matrix[Lrg_AlaskaYukon_AndRam,Lrg_AlaskaYukon_AndRam]
Lrg_AlaskaYukon_AndRam_IndelMatrix <- Comb_PW_Indel_Matrix[Lrg_AlaskaYukon_AndRam,Lrg_AlaskaYukon_AndRam]
Lrg_AlaskaYukon_AndRam_FinalMatrix <- Comb_Final_Matrix[Lrg_AlaskaYukon_AndRam,Lrg_AlaskaYukon_AndRam]
Lrg_AlaskaYukon_AndRam_PW_DistMatrix <- PW_DistMatrix[Lrg_AlaskaYukon_AndRam,Lrg_AlaskaYukon_AndRam]

#Great Lakes Clade (No Virginia)
GreatLakes_DistMatrix <- Comb_PW_Sub_Matrix[GreatLakes,GreatLakes]
GreatLakes_IndelMatrix <- Comb_PW_Indel_Matrix[GreatLakes,GreatLakes]
GreatLakes_FinalMatrix <- Comb_Final_Matrix[GreatLakes,GreatLakes]
GreatLakes_PW_DistMatrix <- PW_DistMatrix[GreatLakes,GreatLakes]

#Virginia 
Virginia_DistMatrix <- Comb_PW_Sub_Matrix[Virginia,Virginia]
Virginia_IndelMatrix <- Comb_PW_Indel_Matrix[Virginia,Virginia]
Virginia_FinalMatrix <- Comb_Final_Matrix[Virginia,Virginia]
Virginia_PW_DistMatrix <- PW_DistMatrix[Virginia,Virginia]

#Great Lakes with Vrignia
GLV_DistMatrix <- Comb_PW_Sub_Matrix[GLV,GLV]
GLV_IndelMatrix <- Comb_PW_Indel_Matrix[GLV,GLV]
GLV_FinalMatrix <- Comb_Final_Matrix[GLV,GLV]
GLV_PW_DistMatrix <- PW_DistMatrix[GLV,GLV]

#AlbertaMissouri
AlbertaMissouri_DistMatrix <- Comb_PW_Sub_Matrix[AlbertaMissouri,AlbertaMissouri]
AlbertaMissouri_IndelMatrix <- Comb_PW_Indel_Matrix[AlbertaMissouri,AlbertaMissouri]
AlbertaMissouri_FinalMatrix <- Comb_Final_Matrix[AlbertaMissouri,AlbertaMissouri]
AlbertaMissouri_PW_DistMatrix <- PW_DistMatrix[AlbertaMissouri,AlbertaMissouri]

#SmallAlaska
SmallAlaska_DistMatrix <- Comb_PW_Sub_Matrix[SmallAlaska,SmallAlaska]
SmallAlaska_IndelMatrix <- Comb_PW_Indel_Matrix[SmallAlaska,SmallAlaska]
SmallAlaska_FinalMatrix <- Comb_Final_Matrix[SmallAlaska,SmallAlaska]
SmallAlaska_PW_DistMatrix <- PW_DistMatrix[SmallAlaska,SmallAlaska]


##
#Filtering Distance Matrices to just lower triangles
#

#Large Alaska Yukon Clade
Lrg_AlaskaYukon_DistMatrix[upper.tri(Lrg_AlaskaYukon_DistMatrix, diag = T)] <- NA
Lrg_AlaskaYukon_IndelMatrix[upper.tri(Lrg_AlaskaYukon_IndelMatrix, diag = T)] <- NA
Lrg_AlaskaYukon_FinalMatrix[upper.tri(Lrg_AlaskaYukon_FinalMatrix, diag = T)] <- NA
Lrg_AlaskaYukon_PW_DistMatrix[upper.tri(Lrg_AlaskaYukon_PW_DistMatrix, diag = T)] <- NA


#Large Alaska Yukon Clade (with RAM P94.16.1B)
Lrg_AlaskaYukon_AndRam_DistMatrix[upper.tri(Lrg_AlaskaYukon_AndRam_DistMatrix, diag = T)] <- NA
Lrg_AlaskaYukon_AndRam_IndelMatrix[upper.tri(Lrg_AlaskaYukon_AndRam_IndelMatrix, diag = T)] <- NA
Lrg_AlaskaYukon_AndRam_FinalMatrix[upper.tri(Lrg_AlaskaYukon_AndRam_FinalMatrix, diag = T)] <- NA
Lrg_AlaskaYukon_AndRam_PW_DistMatrix[upper.tri(Lrg_AlaskaYukon_AndRam_PW_DistMatrix, diag = T)] <- NA


#Great Lakes Clade (No Virginia)
GreatLakes_DistMatrix[upper.tri(GreatLakes_DistMatrix, diag = T)] <- NA
GreatLakes_IndelMatrix[upper.tri(GreatLakes_IndelMatrix, diag = T)] <- NA
GreatLakes_FinalMatrix[upper.tri(GreatLakes_FinalMatrix, diag = T)] <- NA
GreatLakes_PW_DistMatrix[upper.tri(GreatLakes_PW_DistMatrix, diag = T)] <- NA


#Virginia 
Virginia_DistMatrix[upper.tri(Virginia_DistMatrix, diag = T)] <- NA
Virginia_IndelMatrix[upper.tri(Virginia_IndelMatrix, diag = T)] <- NA
Virginia_FinalMatrix[upper.tri(Virginia_FinalMatrix, diag = T)] <- NA
Virginia_PW_DistMatrix[upper.tri(Virginia_PW_DistMatrix, diag = T)] <- NA

#Great Lakes with Vrignia
GLV_DistMatrix[upper.tri(GLV_DistMatrix, diag = T)] <- NA
GLV_IndelMatrix[upper.tri(GLV_IndelMatrix, diag = T)] <- NA
GLV_FinalMatrix[upper.tri(GLV_FinalMatrix, diag = T)] <- NA
GLV_PW_DistMatrix[upper.tri(GLV_PW_DistMatrix, diag = T)] <- NA

#AlbertaMissouri
AlbertaMissouri_DistMatrix[upper.tri(AlbertaMissouri_DistMatrix, diag = T)] <- NA
AlbertaMissouri_IndelMatrix[upper.tri(AlbertaMissouri_IndelMatrix, diag = T)] <- NA
AlbertaMissouri_FinalMatrix[upper.tri(AlbertaMissouri_FinalMatrix, diag = T)] <- NA
AlbertaMissouri_PW_DistMatrix[upper.tri(AlbertaMissouri_PW_DistMatrix, diag = T)] <- NA

#SmallAlaska
SmallAlaska_DistMatrix[upper.tri(SmallAlaska_DistMatrix, diag = T)] <- NA
SmallAlaska_IndelMatrix[upper.tri(SmallAlaska_IndelMatrix, diag = T)] <- NA
SmallAlaska_FinalMatrix[upper.tri(SmallAlaska_FinalMatrix, diag = T)] <- NA
SmallAlaska_PW_DistMatrix[upper.tri(SmallAlaska_PW_DistMatrix, diag = T)] <- NA


##
#Mean and Standard Deviation Calculations
##

#Large Alaska Yukon Clade
Lrg_AlaskaYukon_DistMatrix_Mean <- mean(Lrg_AlaskaYukon_DistMatrix, na.rm=TRUE)
Lrg_AlaskaYukon_DistMatrix_SD <- sd(Lrg_AlaskaYukon_DistMatrix, na.rm=TRUE)

Lrg_AlaskaYukon_IndelMatrix_Mean <- mean(Lrg_AlaskaYukon_IndelMatrix, na.rm=TRUE)
Lrg_AlaskaYukon_IndelMatrix_SD <- sd(Lrg_AlaskaYukon_IndelMatrix, na.rm=TRUE)

Lrg_AlaskaYukon_FinalMatrix_Mean <- mean(Lrg_AlaskaYukon_FinalMatrix, na.rm=TRUE)
Lrg_AlaskaYukon_FinalMatrix_SD <- sd(Lrg_AlaskaYukon_FinalMatrix, na.rm=TRUE)

Lrg_AlaskaYukon_PW_DistMatrix_Mean <- mean(Lrg_AlaskaYukon_PW_DistMatrix, na.rm=TRUE)
Lrg_AlaskaYukon_PW_DistMatrix_SD <- sd(Lrg_AlaskaYukon_PW_DistMatrix, na.rm=TRUE)

#Large Alaska Yukon Clade (with RAM P94.16.1B)
Lrg_AlaskaYukon_AndRam_DistMatrix_Mean <- mean(Lrg_AlaskaYukon_AndRam_DistMatrix, na.rm=TRUE)
Lrg_AlaskaYukon_AndRam_DistMatrix_SD <- sd(Lrg_AlaskaYukon_AndRam_DistMatrix, na.rm=TRUE)

Lrg_AlaskaYukon_AndRam_IndelMatrix_Mean <- mean(Lrg_AlaskaYukon_AndRam_IndelMatrix, na.rm=TRUE)
Lrg_AlaskaYukon_AndRam_IndelMatrix_SD <- sd(Lrg_AlaskaYukon_AndRam_IndelMatrix, na.rm=TRUE)

Lrg_AlaskaYukon_AndRam_FinalMatrix_Mean <- mean(Lrg_AlaskaYukon_AndRam_FinalMatrix, na.rm=TRUE)
Lrg_AlaskaYukon_AndRam_FinalMatrix_SD <- sd(Lrg_AlaskaYukon_AndRam_FinalMatrix, na.rm=TRUE)

Lrg_AlaskaYukon_AndRam_PW_DistMatrix_Mean <- mean(Lrg_AlaskaYukon_AndRam_PW_DistMatrix, na.rm=TRUE)
Lrg_AlaskaYukon_AndRam_PW_DistMatrix_SD <- sd(Lrg_AlaskaYukon_AndRam_PW_DistMatrix, na.rm=TRUE)

#Great Lakes Clade (No Virginia)
GreatLakes_DistMatrix_Mean <- mean(GreatLakes_DistMatrix, na.rm=TRUE)
GreatLakes_DistMatrix_SD <- sd(GreatLakes_DistMatrix, na.rm=TRUE)

GreatLakes_IndelMatrix_Mean <- mean(GreatLakes_IndelMatrix, na.rm=TRUE)
GreatLakes_IndelMatrix_SD <- sd(GreatLakes_IndelMatrix, na.rm=TRUE)

GreatLakes_FinalMatrix_Mean <- mean(GreatLakes_FinalMatrix, na.rm=TRUE)
GreatLakes_FinalMatrix_SD <- sd(GreatLakes_FinalMatrix, na.rm=TRUE)

GreatLakes_PW_DistMatrix_Mean <- mean(GreatLakes_PW_DistMatrix, na.rm=TRUE)
GreatLakes_PW_DistMatrix_SD <- sd(GreatLakes_PW_DistMatrix, na.rm=TRUE)

#Virginia 
Virginia_DistMatrix_Mean <- mean(Virginia_DistMatrix, na.rm=TRUE)
Virginia_DistMatrix_SD <- sd(Virginia_DistMatrix, na.rm=TRUE)

Virginia_IndelMatrix_Mean <- mean(Virginia_IndelMatrix, na.rm=TRUE)
Virginia_IndelMatrix_SD <- sd(Virginia_IndelMatrix, na.rm=TRUE)

Virginia_FinalMatrix_Mean <- mean(Virginia_FinalMatrix, na.rm=TRUE)
Virginia_FinalMatrix_SD <- sd(Virginia_FinalMatrix, na.rm=TRUE)

Virginia_PW_DistMatrix_Mean <- mean(Virginia_PW_DistMatrix, na.rm=TRUE)
Virginia_PW_DistMatrix_SD <- sd(Virginia_PW_DistMatrix, na.rm=TRUE)

#Great Lakes with Vrignia
GLV_DistMatrix_Mean <- mean(GLV_DistMatrix, na.rm=TRUE)
GLV_DistMatrix_SD <- sd(GLV_DistMatrix, na.rm=TRUE)

GLV_IndelMatrix_Mean <- mean(GLV_IndelMatrix, na.rm=TRUE)
GLV_IndelMatrix_SD <- sd(GLV_IndelMatrix, na.rm=TRUE)

GLV_FinalMatrix_Mean <- mean(GLV_FinalMatrix, na.rm=TRUE)
GLV_FinalMatrix_SD <- sd(GLV_FinalMatrix, na.rm=TRUE)

GLV_PW_DistMatrix_Mean <- mean(GLV_PW_DistMatrix, na.rm=TRUE)
GLV_PW_DistMatrix_SD <- sd(GLV_PW_DistMatrix, na.rm=TRUE)

#AlbertaMissouri
AlbertaMissouri_DistMatrix_Mean <- mean(AlbertaMissouri_DistMatrix, na.rm=TRUE)
AlbertaMissouri_DistMatrix_SD <- sd(AlbertaMissouri_DistMatrix, na.rm=TRUE)

AlbertaMissouri_IndelMatrix_Mean <- mean(AlbertaMissouri_IndelMatrix, na.rm=TRUE)
AlbertaMissouri_IndelMatrix_SD <- sd(AlbertaMissouri_IndelMatrix, na.rm=TRUE)

AlbertaMissouri_FinalMatrix_Mean <- mean(AlbertaMissouri_FinalMatrix, na.rm=TRUE)
AlbertaMissouri_FinalMatrix_SD <- sd(AlbertaMissouri_FinalMatrix, na.rm=TRUE)

AlbertaMissouri_PW_DistMatrix_Mean <- mean(AlbertaMissouri_PW_DistMatrix, na.rm=TRUE)
AlbertaMissouri_PW_DistMatrix_SD <- sd(AlbertaMissouri_PW_DistMatrix, na.rm=TRUE)

#SmallAlaska
SmallAlaska_DistMatrix_Mean <- mean(SmallAlaska_DistMatrix, na.rm=TRUE)
SmallAlaska_DistMatrix_SD <- sd(SmallAlaska_DistMatrix, na.rm=TRUE)

SmallAlaska_IndelMatrix_Mean <- mean(SmallAlaska_IndelMatrix, na.rm=TRUE)
SmallAlaska_IndelMatrix_SD <- sd(SmallAlaska_IndelMatrix, na.rm=TRUE)

SmallAlaska_FinalMatrix_Mean <- mean(SmallAlaska_FinalMatrix, na.rm=TRUE)
SmallAlaska_FinalMatrix_SD <- sd(SmallAlaska_FinalMatrix, na.rm=TRUE)

SmallAlaska_PW_DistMatrix_Mean <- mean(SmallAlaska_PW_DistMatrix, na.rm=TRUE)
SmallAlaska_PW_DistMatrix_SD <- sd(SmallAlaska_PW_DistMatrix, na.rm=TRUE)


##Creating the final data frames before writing to file.

FinalVales_DistMatrix <- data.frame("Clade" = c("Large Alaska Yukon", "Large Alaska Yukon and RAM", "Great Lakes (No Virginia)", "Virginia", "Great Lakes and Virginia", "Alberta/Missouri", "Small Alaska"), "Mean SNPs" = c(Lrg_AlaskaYukon_DistMatrix_Mean, Lrg_AlaskaYukon_AndRam_DistMatrix_Mean, GreatLakes_DistMatrix_Mean, Virginia_DistMatrix_Mean, GLV_DistMatrix_Mean, AlbertaMissouri_DistMatrix_Mean, SmallAlaska_DistMatrix_Mean), "St.Dev" = c(Lrg_AlaskaYukon_DistMatrix_SD, Lrg_AlaskaYukon_AndRam_DistMatrix_SD, GreatLakes_DistMatrix_SD, Virginia_DistMatrix_SD, GLV_DistMatrix_SD, AlbertaMissouri_DistMatrix_SD, SmallAlaska_DistMatrix_SD))

FinalVales_IndelMatrix <- data.frame("Clade" = c("Large Alaska Yukon", "Large Alaska Yukon and RAM", "Great Lakes (No Virginia)", "Virginia", "Great Lakes and Virginia", "Alberta/Missouri", "Small Alaska"), "Mean Indels" = c(Lrg_AlaskaYukon_IndelMatrix_Mean, Lrg_AlaskaYukon_AndRam_IndelMatrix_Mean, GreatLakes_IndelMatrix_Mean, Virginia_IndelMatrix_Mean, GLV_IndelMatrix_Mean, AlbertaMissouri_IndelMatrix_Mean, SmallAlaska_IndelMatrix_Mean), "St.Dev" = c(Lrg_AlaskaYukon_IndelMatrix_SD, Lrg_AlaskaYukon_AndRam_IndelMatrix_SD, GreatLakes_IndelMatrix_SD, Virginia_IndelMatrix_SD, GLV_IndelMatrix_SD, AlbertaMissouri_IndelMatrix_SD, SmallAlaska_IndelMatrix_SD))

FinalVales_FinalMatrix <- data.frame("Clade" = c("Large Alaska Yukon", "Large Alaska Yukon and RAM", "Great Lakes (No Virginia)", "Virginia", "Great Lakes and Virginia", "Alberta/Missouri", "Small Alaska"), "Mean Total Differences" = c(Lrg_AlaskaYukon_FinalMatrix_Mean, Lrg_AlaskaYukon_AndRam_FinalMatrix_Mean, GreatLakes_FinalMatrix_Mean, Virginia_FinalMatrix_Mean, GLV_FinalMatrix_Mean, AlbertaMissouri_FinalMatrix_Mean, SmallAlaska_FinalMatrix_Mean), "St.Dev" = c(Lrg_AlaskaYukon_FinalMatrix_SD, Lrg_AlaskaYukon_AndRam_FinalMatrix_SD, GreatLakes_FinalMatrix_SD, Virginia_FinalMatrix_SD, GLV_FinalMatrix_SD, AlbertaMissouri_FinalMatrix_SD, SmallAlaska_FinalMatrix_SD))

FinalVales_PW_DistMatrix <- data.frame("Clade" = c("Large Alaska Yukon", "Large Alaska Yukon and RAM", "Great Lakes (No Virginia)", "Virginia", "Great Lakes and Virginia", "Alberta/Missouri", "Small Alaska"), "Mean Pairwise Distance" = c(Lrg_AlaskaYukon_PW_DistMatrix_Mean, Lrg_AlaskaYukon_AndRam_PW_DistMatrix_Mean, GreatLakes_PW_DistMatrix_Mean, Virginia_PW_DistMatrix_Mean, GLV_PW_DistMatrix_Mean, AlbertaMissouri_PW_DistMatrix_Mean, SmallAlaska_PW_DistMatrix_Mean), "St.Dev" = c(Lrg_AlaskaYukon_PW_DistMatrix_SD, Lrg_AlaskaYukon_AndRam_PW_DistMatrix_SD, GreatLakes_PW_DistMatrix_SD, Virginia_PW_DistMatrix_SD, GLV_PW_DistMatrix_SD, AlbertaMissouri_PW_DistMatrix_SD, SmallAlaska_PW_DistMatrix_SD))


#Writing the output file.
write.table("Pairwise_Substitution_Matrix", file="Masto_Clade_Distances.csv", row.names = F, col.names = F)
write.table(FinalVales_DistMatrix, file="Masto_Clade_Distances_withPWDist_F84.csv", append = T, sep=",", row.names = T, col.names = T)

write.table("Indel_Matrix", file="Masto_Clade_Distances.csv",append = T, row.names = F, col.names = F)
write.table(FinalVales_IndelMatrix, file="Masto_Clade_Distances_withPWDist_F84.csv", append = T, sep=",", row.names = T, col.names = T)


write.table("Final_Matrix", file="Masto_Clade_Distances.csv",append = T, row.names = F, col.names = F)
write.table(FinalVales_FinalMatrix, file="Masto_Clade_Distances_withPWDist_F84.csv", append = T, sep=",", row.names = T, col.names = T)

write.table("PW_Dist_Matrix", file="Masto_Clade_Distances.csv",append = T, row.names = F, col.names = F)
write.table(FinalVales_PW_DistMatrix, file="Masto_Clade_Distances_withPWDist_F84.csv", append = T, sep=",", row.names = T, col.names = T)

