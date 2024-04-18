### Statistical analysis script in R v2023.06.1+524 

## Analysis of variance
# Data files are available in this Github project, and on Wilderlab NZ explore map (Accession codes: 603518, 603589, and 605453).
##########################
## Analysis of variance ##
##########################
# Analysis of variance for the 1.2 µm filter at the start of spring
SS1_12_start <- read.csv(file="SS1-12-start.csv")
res.aov <- aov(Count ~ Treatment, data = SS1_12_start)
summary(res.aov)
# Repeated for each set of eDNA sample replicates using corresponding data

## Non-metric multidimensional plots
# Data files are available in this Github project, and on Wilderlab NZ explore map (Accession codes: 603518, 603589, and 605453).
#######################################
## Non-metric multidimensional plots ##
#######################################
# Load in the R packages
library(vegan)
library(pander)
library(ggplot2)
library(BiodiversityR)
library("readxl")
# Read in the data
x1 <- read_excel("../WLJ603518startfilters.xlsx", sheet = 3)
x1 <- x1[, c(1,8:ncol(x1))]
x2 <- read_excel("../WLJ603518midfilters.xlsx", sheet = 3)
x2 <- x2[, c(1,8:ncol(x2))]
x3 <- read_excel("../WLJ603589endfilters.xlsx", sheet = 3)
x3 <- x3[, c(1,8:ncol(x3))]
# Merge start, middle, and end of spring data, and assign 'NA' to zero 
x <- merge(x1, x2, by = "Sequence")
x <- merge(x, x3, by= "Sequence", all = T)
x[is.na(x)] <- 0
# Read in the metadata
meta1 <- read_excel("../WLJ603518startfilters.xlsx", sheet = 1)
meta2 <- read_excel("../WLJ603518midfilters.xlsx", sheet = 1)
meta3 <- read_excel("../WLJ603589endfilters.xlsx", sheet = 1)
# Merge start, middle, and end of spring metadata
meta1 <- meta1[,c(1:2,11)]
colnames(meta1) <- c("UID", "metadata", "SpringTime")
meta2 <- meta2[,c(1:2,ncol(meta2))]
colnames(meta2) <- c("UID", "metadata", "SpringTime")
meta3 <- meta3[,c(1:2,ncol(meta3))]
colnames(meta3) <- c("UID", "metadata", "SpringTime")
metadata <- rbind(meta1, meta2)
metadata <- rbind(metadata, meta3)
metadata$SpringTime[metadata$SpringTime == 'Mid'] <- 'Middle'
# Working on merged data
community_matrix <- x[seq(2, ncol(x))]
UIDs <- colnames(community_matrix)
community_matrix <- t(as.matrix(community_matrix))
colnames(community_matrix) <- x$Sequence
rownames(community_matrix) <- UIDs
colsums <- apply(community_matrix, 2, function(v) sum(v > 0))
community_matrix <- community_matrix[, colsums >= 7L]
mymds <- vegan::metaMDS(community_matrix, weakties = FALSE, autotransform = F)
mysamples <- as.data.frame(scores(mymds, "sites"))
mysamples$UID <- rownames(mysamples)
mysamples$Sites <- metadata$metadata[match(mysamples$UID, metadata$UID)]
mysamples$SpringTime <- metadata$SpringTime[match(mysamples$UID, metadata$UID)]
myseqs <- as.data.frame(scores(mymds, "species"))
myseqs$species <- rownames(myseqs)
rownames(myseqs) <- NULL
# Reorder SpringTime for ggplot 
mysamples$SpringTime <- factor(mysamples$SpringTime, levels = c("Start","Middle","End"))
# Plot filter size with gradient sized legends
myplot <- ggplot(mysamples, aes(x=NMDS1, y=NMDS2, shape=SpringTime, color=SpringTime, size = factor(Sites, ordered = T))) +
  geom_point() +
  scale_color_manual(values = c("indianred1", "indianred", "indianred4")) + 
  labs(shape="Spring period", color="Spring period", size="Filter size (µm)") + 
  scale_size_manual(values= c(2.5, 5, 9)) +
  xlim(c(-1.04,1)) +
  ylim(c(-1.06,1)) + 
  guides(shape = guide_legend(order = 1, override.aes = list(size = 5)),
         col = guide_legend(order = 1, override.aes = list(size = 5)),
         size = guide_legend(order = 2)) +
  theme(
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.key = element_rect(size = 6),
    legend.key.size = unit(1, "cm"),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14), 
    axis.text.x = element_text(size = 14),  
    axis.text.y = element_text(size = 14)) +
  coord_equal() 
myplot
# Repeated for spatial and temporal analysis using corresponding data for the full dataset, at species and genus level, and key wetland species of interest

## PERMANOVA and PERMDISP
# Data files are available in this Github project, and on Wilderlab NZ explore map (Accession codes: 603518, 603589, and 605453).
############################
## PERMANOVA and PERMDISP ##
############################
# Load in the R packages
library(vegan)
# Formatting frequency table
pres_abs <- read.csv("WLJ603518startfilters_perm.csv")
pres_abs1 <- pres_abs[,1] 
pres_abs2 <- read.csv("WLJ603518startfilters_perm.csv")
pres_abs3 <- pres_abs2[,2:17]
pres_abs3[pres_abs3 > 0] <- 1
pres_abs_final <- cbind(pres_abs1, pres_abs3)
perm <- t(pres_abs_final) 
colnames(perm) <- perm[1,]
perm <- perm[-1, ]
write.csv(perm, file="finalmat_StartFilter.csv", row.names=T)
# Read in the data
metadata <- read.csv("metadata_startfilters_perm.csv")
StartFilters <- read.csv("finalmat_StartFilter.csv", header = TRUE, row.names = 1)
# Run PERMANOVA
adonis2(StartFilters ~ Filter, data = metadata, by = 'margin') 
# Test for dispersion effects
dis <- vegdist(StartFilters, method = 'jaccard')
groups <- factor(c(rep(1,5), rep(2,5), rep(3,6)),
                 labels = c('1.2', '5', 'Dacron'))
dis
groups
mod <- betadisper(dis, groups)
mod
boxplot(mod)
set.seed(25)
permutest(mod)
plot(mod)
set.seed(4)
permutest(mod, pairwise = TRUE)
# Repeated for the middle and end of spring filter experiments, as well as spatial and temporal analysis using corresponding data for the full dataset, at species and genus level, and key wetland species of interest

## Species accumulation curves
# Data files are available in this Github project, and on Wilderlab NZ explore map (Accession codes: 603518, 603589, and 605453).
#################################
## Species accumulation curves ##
#################################
# Load in the R packages
library(vegan)
# Formatting frequency table for 1.2 µm filter
StartOne <- read.csv("WLJ603518StartOne_SpAccu.csv")
StartOne1 <- StartOne[,1]
StartOne2 <- read.csv("WLJ603518StartOne_SpAccu.csv")
StartOne3 <- StartOne2[,2:6]
StartOne_final <- cbind(StartOne1, StartOne3)
StartOne_SpA <- t(StartOne_final) 
colnames(StartOne_SpA) <- StartOne_SpA[1,] 
StartOne_SpA <- StartOne_SpA[-1, ]
write.csv(StartOne_SpA, file="Finalmat_StartFilterOne_SpAccu.csv", row.names=T)
# Specaccum analysis 1.2 µm filter
StartOne_spec <- read.csv("Finalmat_StartFilterOne_SpAccu.csv")
sp1 <- specaccum(StartOne_spec)
sp2 <- specaccum(StartOne_spec, "random")
sp2
summary(sp2)
# Formatting frequency table for 5 µm filter
StartFive <- read.csv("WLJ603518StartFive_SpAccu.csv")
StartFive1 <- StartFive[,1] 
StartFive2 <- read.csv("WLJ603518StartFive_SpAccu.csv")
StartFive3 <- StartFive2[,2:6] 
StartFive_final <- cbind(StartFive1, StartFive3)
StartFive_SpA <- t(StartFive_final) 
colnames(StartFive_SpA) <- StartFive_SpA[1,] 
StartFive_SpA <- StartFive_SpA[-1, ]
write.csv(StartFive_SpA, file="Finalmat_StartFilterFive_SpAccu.csv", row.names=T)
# Specaccum analysis 5 µm filter
StartFive_spec <- read.csv("Finalmat_StartFilterFive_SpAccu.csv")
sp3 <- specaccum(StartFive_spec)
sp4 <- specaccum(StartFive_spec, "random")
sp4
summary(sp4)
# Formatting frequency table for dacron filter
StartDacron <- read.csv("WLJ603518StartDacron_SpAccu.csv")
StartDacron1 <- StartDacron[,1] 
StartDacron2 <- read.csv("WLJ603518StartDacron_SpAccu.csv")
StartDacron3 <- StartDacron2[,2:7] 
StartDacron_final <- cbind(StartDacron1, StartDacron3)
StartDacron_SpA <- t(StartDacron_final) 
colnames(StartDacron_SpA) <- StartDacron_SpA[1,] 
StartDacron_SpA <- StartDacron_SpA[-1, ] 
write.csv(StartDacron_SpA, file="Finalmat_StartFilterDacron_SpAccu.csv", row.names=T) 
# Specaccum analysis 5 µm filter for dacron filter
StartDacron_spec <- read.csv("Finalmat_StartFilterDacron_SpAccu.csv")
sp5 <- specaccum(StartDacron_spec)
sp6 <- specaccum(StartDacron_spec, "random")
sp6
summary(sp6)
# Formatting plot layout: three plots on same window
par(mfrow = c(1, 3), mar = c(5, 4, 1, 0.1))
# Define custom y-axis tick values and labels
y_ticks <- c(0, 250, 500, 750, 1000, 1250, 1500, 1750, 2000)
y_labels <- y_ticks
x_ticks <- c(1, 2, 3, 4, 5)
x_labels <- x_ticks
x2_ticks <- c(1, 2, 3, 4, 5, 6)
x2_labels <- x2_ticks
# 1.2 µm filter species accumulation curve
par(mgp=c(3, 0.6, 0), cex.lab = 1.3)
plot(sp1, ci.type="poly", col="darkmagenta", lwd=2, ci.lty=0, ci.col="lavender", xlim = c(1,5), ylim = c(0,2000), xlab = "", ylab = "Species richness", las = 1, xaxt = "n", yaxt = "n")  
axis(2, at = y_ticks, labels = y_labels, las = 1, cex.axis = 1.2)
axis(1, at = x_ticks, labels = x_labels, cex.axis = 1.2)
boxplot(sp2, col="lightgrey", add=TRUE, pch="+", boxwex = 0.25)
# 5 µm filter species accumulation curve
par(mgp=c(3, 0.6, 0), cex.lab = 1.3)
plot(sp3, ci.type="poly", col="darkviolet", lwd=2, ci.lty=0, ci.col="plum", xlim = c(1,5),  ylim = c(0,2000), xaxt = "n", yaxt = "n", xlab = "Replicates", ylab = "", las = 1)
axis(1, at = x_ticks, labels = x_labels, cex.axis = 1.2)
boxplot(sp4, col="lightgrey", add=TRUE, pch="+", boxwex = 0.25)
# Dacron filter species accumulation curve
par(mgp=c(3, 0.6, 0))
plot(sp5, ci.type="poly", col="#e75480", lwd=2, ci.lty=0, ci.col="pink", xlim = c(1,6),  ylim = c(0,2000), xaxt = "n", yaxt = "n", xlab = "", ylab = "", las = 1)
axis(1, at = x2_ticks, labels = x2_labels, cex.axis = 1.2)
boxplot(sp6, col="lightgrey", add=TRUE, pch="+", boxwex = 0.25)
# Repeated for the middle and end of spring filter experiments, as well as the number of sites sampled using the corresponding data
