##### mats prep and some descriptive stats ######

# read the OTU table
otumat <- read.table("input/OTU.txt", header = T, sep = "\t", row.names = 1)
# read the sequence hierarchical classification
hiera <- read.table("input/hiera_BLAST.txt", header = T, sep = "\t", row.names = 1, quote = "", comment.char = "", stringsAsFactors = F)
# read the experimental design
design <- read.table("input/design.txt", header = T, sep = "\t", row.names = 1)

# change the OTU names by using the lowest confident annotation
my_tax_tbl <- hiera
# select the taxa of interest
my_tax_tbl_2 <- my_tax_tbl[which(my_tax_tbl$Domain == "Archaea"),]
## Fix the annotation
# E.g. one issue with the fungal annotation is that they use unidentified also next to the ? and they add in the end of the unidientified  the lowest classified level along with unidentified XXX sp... I had to replace these annotations with ? in order to give the lowest possible classifications to the barplots 
my_tax_tbl_3 <- my_tax_tbl_2
for(myrow in 1:nrow(my_tax_tbl_3)){
  if(my_tax_tbl_3[myrow,"Genus"] == "unidentified"){
    my_tax_tbl_3[myrow,"Species"] <- "unidentified"
  }
}
## replace the "unidentified" with the unclassified leftmost (higher level annotated) taxa
my_tax_tbl_3b <- my_tax_tbl_3
for(myrow in 1:nrow(my_tax_tbl_3b)){
  for(mycol in 2:ncol(my_tax_tbl_3b)){
    if(my_tax_tbl_3b[myrow,mycol] == "unidentified" & length(grep("uncl\\.",my_tax_tbl_3b[myrow,mycol-1])) > 0){
      my_tax_tbl_3b[myrow,mycol] <- my_tax_tbl_3b[myrow,mycol-1]
    } else if (my_tax_tbl_3b[myrow,mycol] == "unidentified" & length(grep("uncl\\.",my_tax_tbl_3b[myrow,mycol-1])) == 0) {
      my_tax_tbl_3b[myrow,mycol] <- paste("uncl.",my_tax_tbl_3b[myrow,mycol-1])
    } else {
      my_tax_tbl_3b[myrow,mycol] <- my_tax_tbl_3b[myrow,mycol]
    }
  }
}
## replace the ? with the unclassified leftmost taxa
my_tax_tbl_3c <- my_tax_tbl_3b
for(myrow in 1:nrow(my_tax_tbl_3c)){
  for(mycol in 2:ncol(my_tax_tbl_3c)){
    if(my_tax_tbl_3c[myrow,mycol] == "?" & length(grep("uncl\\.",my_tax_tbl_3c[myrow,mycol-1])) > 0){
      my_tax_tbl_3c[myrow,mycol] <- my_tax_tbl_3c[myrow,mycol-1]
    } else if (my_tax_tbl_3c[myrow,mycol] == "?" & length(grep("uncl\\.",my_tax_tbl_3c[myrow,mycol-1])) == 0) {
      my_tax_tbl_3c[myrow,mycol] <- paste("uncl.",my_tax_tbl_3c[myrow,mycol-1])
    } else {
      my_tax_tbl_3c[myrow,mycol] <- my_tax_tbl_3c[myrow,mycol]
    }
  }
}
## replace the Incertae sedis with the taxon in the left plus IS (for easier visualization)
my_tax_tbl_4 <- my_tax_tbl_3c
for(myrow in 1:nrow(my_tax_tbl_4)){
  for(mycol in 2:ncol(my_tax_tbl_4)){
    if(my_tax_tbl_4[myrow,mycol] == "Incertae sedis"){
      my_tax_tbl_4[myrow,mycol] <- paste(my_tax_tbl_4[myrow,mycol-1], " IS", sep = "")
    } else {
      my_tax_tbl_4[myrow,mycol] <- my_tax_tbl_4[myrow,mycol]
    }
  }
}
my_tax_tbl_5 <- my_tax_tbl_4
# replace double IS with a single one coming up in some cases
for(myrow in 1:nrow(my_tax_tbl_5)){
  for(mycol in 2:ncol(my_tax_tbl_5)){
    my_tax_tbl_5[myrow,mycol] <- gsub("IS IS","IS", my_tax_tbl_5[myrow,mycol])
  }
}
# assign the table to a final name
my_tax_tbl_fin <- my_tax_tbl_5

# create a table containing both OTU abundances and their taxonomy and set the row.names
combomat <- merge(otumat, my_tax_tbl_5, by = "row.names", all.x = T)
row.names(combomat) <- combomat$Row.names
# select for the taxa of interest only
combomatfin <- combomat[grep("Archaea", combomat$Domain),]
# order the previous matrix and assign it to a new table
otumatfin <- combomatfin[order(as.numeric(gsub("OTU_","",row.names(combomatfin)))),grep(paste(row.names(design), collapse = "|", sep = ""), colnames(combomatfin))]
# check for empty OTUs (they should'n exist but good to validate)
if(length(which(rowSums(otumatfin) == 0)) > 0){
  otumatfin <- otumatfin[-which(rowSums(otumatfin) == 0),]
} else {
  otumatfin <- otumatfin
}

####### create an output directory #####
system("mkdir output")

####### rarefaction curves #########
library(vegan)
# obtain the minimum number of events throughout all samples in order to use as a subsample size for rarefying the community
raremax <- min(rowSums(t(otumatfin)[order(row.names(t(otumatfin))),]))
# obtain the rarefaction-based expected species richness
Srare <- rarefy(t(otumatfin)[order(row.names(t(otumatfin))),], raremax)
# get the rarefaction curve parameters
myrarecurve <- rarecurve(t(otumatfin)[order(row.names(t(otumatfin))),],step=50, sample = raremax, col = factor(design[order(row.names(t(otumatfin))),]$plant), cex = 0.6, xlim= c(0,150))
# set the factor for grouping the rarefaction curves
desgnfrrar <- factor(design[row.names(t(otumatfin)[order(row.names(t(otumatfin))),]),]$plant)
# initiate the graphics device for the rarefaction plot
cairo_pdf("output/rarefaction_curve.pdf", width = 7, height = 5)
# plot the first line
plot(x=as.numeric(gsub("N","",names(myrarecurve[[1]]))),y=myrarecurve[[1]], type = "lines", col = desgnfrrar[1], xlim=c(0,150), ylim=c(0,20), frame.plot = F, ylab = "OTUs", xlab = "Sequence #")
# add the rest rarefaction lines
for(i in 2:length(myrarecurve)) 
{
  # i-th element of `u1` squared into `i`-th position of `usq`
  lines(x=as.numeric(gsub("N","",names(myrarecurve[[i]]))),y=myrarecurve[[i]], col = desgnfrrar[i], add=T)
}
# add the legend/key
legend("topleft",legend = gsub("_"," ",levels(desgnfrrar)), lty = 1, lwd = 2, col = 1:length(levels(desgnfrrar)),cex=0.8, bty = "n")
# close the device for saving the plot
dev.off()







############### NMDS ###############

# assign the table to be used for the NMDS
mymicrmat <- t(otumatfin)
# run NMDS
mynmds <- metaMDS(mymicrmat)
# obtain the site scores
mynmdssit <- data.frame(scores(mynmds, display = "sites"))
# also the species scores (although we will not use it here)
mynmdsspe <- data.frame(scores(mynmds, display = "species"))
# obtain also the matrix the taxon portions for selecting the most abundant taxa to plot (not used here) 
mymicrmatre <- decostand(mymicrmat, "total")
# set the number of dominant taxa (to be plotted)
domOTUnums <- 15
# sort them by dataset-wise relative abundance
mymicrmatredomOTUs <- names(sort(colMeans(mymicrmatre), decreasing = T))[1:domOTUnums]
# calculate the lowest abundance of the selected OTUs
minradomOTUs <- round(100*colMeans(mymicrmatre)[mymicrmatredomOTUs[domOTUnums]],1)
# calculate the highest abundance of the selected OTUs
maxradomOTUs <- round(100*colMeans(mymicrmatre)[1],1)
# create a grouping factor according to the plant and season interaction
mynmdssit$fact <- gsub("_[a-z]$","",row.names(mynmdssit))
mynmdssit$fact <- as.factor(mynmdssit$fact)
# change the factor labels to assist with the selection of plotting symbols etc.
library(plyr)
meannms1 <- mapvalues(mynmdssit$fact, from=levels(mynmdssit$fact), to=1:length(levels(mynmdssit$fact)))
# calculate the plat/season interaction related centroids and create a dataframe
mynmdssitmeanssdpre <- aggregate(. ~ fact, mynmdssit, function(x) c(mean = mean(x), sd = sd(x)))
mynmdssitmeanssd <- data.frame(mynmdssitmeanssdpre$fact, mynmdssitmeanssdpre$NMDS1, mynmdssitmeanssdpre$NMDS2)
# convert the interaction factors to numbers useful for colouring
meannams2 <- mapvalues(mynmdssitmeanssd$mynmdssitmeanssdpre.fact, from = levels(mynmdssitmeanssd$mynmdssitmeanssdpre.fact), to=1:length(levels(mynmdssitmeanssd$mynmdssitmeanssdpre.fact)))
# create a colour list to plot the 16 different plant season combinations
myplotcols <- paste(rep(RColorBrewer::brewer.pal(8,"Dark2"),each=2), c("FF","B3"), sep = "")

# initiate the graphics device
cairo_pdf("output/NMDS.pdf", height = 8, width = 8)
# prep the plot (do not add the points yet in order to add the grouping ellipses first and add the points on top of the ellipses)
plot(mynmdssit[,1:2], bg = myplotcols[meannms1], frame = F, cex = 0, pch = 21, cex.main = 3)
# add some ellipses (also vegan has the ordiellipse function that could be used instead)
car::dataEllipse(as.matrix(mynmdssitmeanssd[,c(2,4)]), groups = factor(gsub("^.+_","",mynmdssitmeanssd$mynmdssitmeanssdpre.fact)), add = T, levels= 0.95, center.pch="", col = rep("black",2), lty = 2 , plot.points=F, group.labels = "")
# add the points
points(mynmdssit[,1:2], bg = myplotcols[meannms1], cex = 1.5, pch = 21)
# add the subtitle key for the plotted OTU gradients
par(adj = 0)
title(sub = paste("stress ", round(mynmds$stress,2), " (mean RA of presented OTUs ",minradomOTUs,"-",maxradomOTUs,"%)", sep = ""))
par(adj = .5)
# add the plant/season interaction centroids
text(mynmdssitmeanssd[,c(2,4)], labels = gsub("_"," ",gsub("Winter","W",gsub("Summer","S",mynmdssitmeanssd$mynmdssitmeanssdpre.fact))), pos = 3, col = myplotcols[meannams2])
# plot also the taxa... although too much to digest by a single pass, it provides some information on the gradients
arrows(0,0,mynmdsspe[mymicrmatredomOTUs,1] , mynmdsspe[mymicrmatredomOTUs,2], angle = 25, length = 0.15)
mynmdslabs <- gsub("___"," ",gsub("___.+___"," ",gsub("___[0-9].+$","",gsub("___uncultured [a-z, A-Z]+","",gsub("[___\\?]+$","",do.call(paste, c(combomatfin[mymicrmatredomOTUs,grep("Row.names|Domain|Phylum|Class|Order|Family|Genus", colnames(combomatfin))], sep="___")))))))
names(mynmdslabs) <- gsub(" .+","",mynmdslabs)
# plot the labels using the thigmophobe command which attempts to avoid overlaps as much as possible
plotrix::thigmophobe.labels(1.2*mynmdsspe[mymicrmatredomOTUs,1], 1.2*mynmdsspe[mymicrmatredomOTUs,2], labels = mynmdslabs[mymicrmatredomOTUs], cex = .7, font = 2, col = "grey60")
# add the season centroids
text(mean(mynmdssitmeanssd[which(factor(gsub("^.+_","",mynmdssitmeanssd$mynmdssitmeanssdpre.fact)) == "Summer"),c(2)]),mean(mynmdssitmeanssd[which(factor(gsub("^.+_","",mynmdssitmeanssd$mynmdssitmeanssdpre.fact)) == "Summer"),c(2)]), labels = "Summer", cex = 1.5, font = 2)
text(mean(mynmdssitmeanssd[which(factor(gsub("^.+_","",mynmdssitmeanssd$mynmdssitmeanssdpre.fact)) == "Winter"),c(2)]),mean(mynmdssitmeanssd[which(factor(gsub("^.+_","",mynmdssitmeanssd$mynmdssitmeanssdpre.fact)) == "Winter"),c(2)]), labels = "Winter", cex = 1.5, font = 2)
# shut the graphics device down
dev.off()




############ canonical analysis ############

### according to the plant aromatic character
# set the necessary tables
mycanmat <- t(otumatfin)
design2 <- design
# convert all character columns into factors
design2[sapply(design2,is.character)] <- lapply(design2[sapply(design2,is.character)],as.factor)

library(vegan)
#transform the with the Hellinger method values
mycanmattr <- decostand(mycanmat, method = "hellinger")
# Run DCA and calculate its first axis length
# and set the cutoff for the selection of the most appropriate method according to the gradients according to the 1st DCA axis 2.5-3 SDs suggested by Leps and Smilauer 2003 (the CANOCO manual) 
DCA <- decorana(mycanmattr)
DCA_1st_axis <- max(scores(DCA,display="sites",origin=FALSE)[,"DCA1"]) - min(scores(DCA,display="sites",origin=FALSE)[,"DCA1"])
myDCAaxiscutoff = 3
# set the conditional statement and run either method depending on the DCA outcome
if(DCA_1st_axis > myDCAaxiscutoff) {  
  #CCA
  myCCA <- cca(mycanmattr ~ Aromatic, design2)
  mysamp <- scores(myCCA, choices = c(1,2), display = "sites")

} else {
  #RDA
  myRDA <- rda(formula = mycanmattr ~ Aromatic, data = design2, scale = T) 
  mysamp <- scores(myRDA, choices = c(1,2), display = "sites")

}
# run the accompanying PERMANOVA according to plant aromatic character and save the output
mypermanova <- adonis2(mycanmattr ~ Aromatic, design2, method = "bray", permutations = 999)
capture.output(mypermanova, file = "output/permanova_aromatic.txt")

# prepare the canonical analysis plot 
cairo_pdf(height = 5, width = 6, file = "output/canonical_aromatic-season.pdf")
# add the PERMANOVA output on top of the plot
mymain <- paste("comm. ~ Aromatic character, shar. var. ",100*round(mypermanova$R2[1],3),"% (p-value ",round(mypermanova$`Pr(>F)`,3)[1],")", sep = "")
# plot the canonical analysis results as in NMDS but without the species scores
plot(mysamp, type = "n", xlim = c(-5,5), ylim = c(-5,7), bty="n", main = mymain, cex.main = 0.8)
points(mysamp, pch=c(21,24)[as.numeric(design2$season)], col="black" , bg = c("black","red")[as.numeric(design2$Aromatic)])
ordiellipse(mysamp, design2$Aromatic, col = c("black"), lty = c(1,2), lwd = 2, label = FALSE)
# prepare the colour/ellipse keys/legends
legend(x = -4, y = -0.7, cex=0.7, col = c("black"), lty = c(1,2), lwd = 2, bty = "n", legend = levels(design2$Aromatic), y.intersp = 1, text.font = 1)
legend(x = -4, y = -2, legend = levels(design2$Aromatic), cex=0.7, col=1:2, pch = 19, bty= "n", bg = "transparent", text.font = 1, y.intersp = 1)
legend(x = -4, y = -3.3, legend = levels(design2$season), cex=0.7, col="black", pch = c(21,24), bty = "n", text.font = 1, y.intersp = 1)
# close the graphics device
dev.off()





### same as above with colour grouping according to the plant canopy type
mycanmat <- t(otumatfin)
design2 <- design
# convert all character columns into factors
design2[sapply(design2,is.character)] <- lapply(design2[sapply(design2,is.character)],as.factor)
design2$type_canopy <- factor(design2$type_canopy, levels = c("Low_canopy","High_canopy"))
library(vegan)
#transform the with the Hellinger method values
mycanmattr <- decostand(mycanmat, method = "hellinger")
# calculate the DCA first axis length
# and set the cutoff for the selection of the most appropriate method according to the gradients according to the 1st DCA axis 2.5-3 SDs suggested by Leps and Smilauer 2003 (the CANOCO manual) 
DCA <- decorana(mycanmattr)
DCA_1st_axis <- max(scores(DCA,display="sites",origin=FALSE)[,"DCA1"]) - min(scores(DCA,display="sites",origin=FALSE)[,"DCA1"])
myDCAaxiscutoff = 2.5
if(DCA_1st_axis > myDCAaxiscutoff) {  
  #CCA
  myCCA <- cca(mycanmattr ~ type_canopy, design2)
  mysamp <- scores(myCCA, choices = c(1,2), display = "sites")
} else {
  #RDA
  myRDA <- rda(formula = mycanmattr ~ type_canopy, data = design2, scale = T) 
  mysamp <- scores(myRDA, choices = c(1,2), display = "sites")
}
mypermanova <- adonis2(mycanmattr ~ type_canopy, design2, method = "bray", permutations = 999)
capture.output(mypermanova, file = "output/permanova_canopy.txt")
# prepare the canonical analysis plot and the 
cairo_pdf(height = 5, width = 6, file = "output/canonical_canopy-season.pdf")
mymain <- paste("comm. ~ canopy type, shar. var. ",100*round(mypermanova$R2[1],3),"% (p-value ",round(mypermanova$`Pr(>F)`,3)[1],")", sep = "")
plot(mysamp, type = "n", xlim = c(-5,5), ylim = c(-5,7), bty="n", main = mymain, cex.main = 0.8)
points(mysamp, pch=c(21,24)[as.numeric(design2$season)], col="black" , bg = c("black","red")[as.numeric(design2$type_canopy)])
ordiellipse(mysamp, design2$type_canopy, col = c("black"), lty = c(1,2), lwd = 2, label = FALSE)
# prepare the colour/ellipse keys/legends
legend(x = -4, y = -0.7, cex=0.7, col = c("black"), lty = c(1,2), lwd = 2, bty = "n", legend = levels(design2$type_canopy), y.intersp = 1, text.font = 1)
legend(x = -4, y = -2, legend = levels(design2$type_canopy), cex=0.7, col=1:2, pch = 19, bty= "n", bg = "transparent", text.font = 1, y.intersp = 1)
legend(x = -4, y = -3.3, legend = levels(design2$season), cex=0.7, col="black", pch = c(21,24), bty = "n", text.font = 1, y.intersp = 1)

dev.off()




### Canonical analysis according to the season and plant
# load the table
mymicrmat <- t(otumatfin)
library(vegan)
# transform the table with the Hellinger transformation according to Legrndre and Gallagher (2001; doi 10.1007/s004420100716)
mymicrmathel <- decostand(mymicrmat, method = "hellinger")
# run the detrended correspondence analysis to select the canonical analysis type according to the DCA first axis length
mydecorana <- vegan::decorana(mymicrmathel)
mydecorana_1st_axis <- max(scores(mydecorana,display="sites",origin=FALSE)[,"DCA1"]) - min(scores(mydecorana,display="sites",origin=FALSE)[,"DCA1"])
myDCAaxiscutoff = 2.5

# Perform CCA or RDA depending on the 1st DCA axis length
if(mydecorana_1st_axis > myDCAaxiscutoff) {  
  eval(parse(text = paste("sol <- cca(mymicrmathel ~ season + plant ,data=design2)", sep = "")))

  mylresultslist <- list(DCA = mydecorana, message = paste("1st axis length more than",myDCAaxiscutoff, "SDs... I am using CCA, change myDCAaxiscutoff if otherwise desired"))
  
} else {
  eval(parse(text = paste("sol <- rda(mymicrmathel ~ season + plant ,data=design2)", sep = "")))

  mylresultslist <- list(DCA = mydecorana, message = paste("1st axis length less than",myDCAaxiscutoff, "SDs... I am using RDA, change myDCAaxiscutoff if otherwise desired"))
  
}
# save the canonical analysis summary
mysum <- summary(sol)
mylresultslist[["summary"]] <- mysum
# obtain the various scores for plotting
scrs <- vegan::scores(sol, display = c("sp", "wa", "lc", "bp", "cn"))
df_sites <- data.frame(scrs$sites, design[, which(colnames(design)%in%c("season", "plant"))])
df_specs <- data.frame(scrs$species)
colnames(df_sites)[1:2] <- c("x", "y")
# start preparing the plot
library(ggplot2)
p <- ggplot()
p <- p + geom_point(data = df_sites, aes(x, y, fill = plant), shape = 21, size=3)
p <- p + stat_ellipse(type = "norm", data = df_sites, aes(x, y, colour = season), linetype = 2)

# get the model variance and add it to the plot (PERMANOVA was used instead of the model below... but kept this part for saving the outputs)
mylresultslist[["AOV-type-output"]] <- mysolutaov <- anova(sol)
mylresultslist[["AOV-type-output-terms"]] <- mysolutaov_terms <- anova(sol, by = "terms", permutations = 999)
# if(mydecorana_1st_axis > myDCAaxiscutoff){
#   mod_exp_var <- round(100*mysolutaov$ChiSquare[1]/sum(mysolutaov$ChiSquare),1)
#   time_exp_var <- round(100*mysolutaov_terms$ChiSquare[1]/sum(mysolutaov_terms$ChiSquare[1:2]),1)
#   plnt_exp_var <- round(100*mysolutaov_terms$ChiSquare[2]/sum(mysolutaov_terms$ChiSquare[1:2]),1)
# } else {
#   mod_exp_var <- round(100*mysolutaov$Variance[1]/sum(mysolutaov$Variance),1)
#   time_exp_var <- round(100*mysolutaov_terms$Variance[1]/sum(mysolutaov_terms$Variance[1:2]),1)
#   plnt_exp_var <- round(100*mysolutaov_terms$Variance[2]/sum(mysolutaov_terms$Variance[1:2]),1)
# }
# mod_p_val <- round(mysolutaov$`Pr(>F)`[1],3)
# time_p_val <- round(mysolutaov_terms$`Pr(>F)`[1],3)
# plnt_p_val <- round(mysolutaov_terms$`Pr(>F)`[2],3)
# 
mypermanova <- adonis2(mymicrmathel ~ season + plant, design2, method = "bray", permutations = 999)

# add the PERMANOVA information
p <- p + labs(title = paste("comm. ~ season + plant, shared variace ",100*round(sum(mypermanova$R2[1:2]),1),"% \nvar. part. season/plant ",round(100*(mypermanova$R2[1]/sum(mypermanova$R2[1:2])),1),"/",round(100*(mypermanova$R2[2]/sum(mypermanova$R2[1:2])),1)," % and P-val. ",round(mypermanova$`Pr(>F)`[1],3),"/",round(mypermanova$`Pr(>F)`[2],3), sep = ""))

# add colours and labels
if(!is.null(myplotcols)) {
  p <- p + theme_bw() + xlab(paste(colnames(mysum$sites)[1]," (",round(100*mysum$cont$importance[2,1],1),"%)", sep = "")) + ylab(paste(colnames(mysum$sites)[2]," (",round(100*mysum$cont$importance[2,2],1),"%)", sep = "")) + scale_fill_manual(values=myplotcols)
} else {
  p <- p + theme_bw() + xlab(paste(colnames(mysum$sites)[1]," (",round(100*mysum$cont$importance[2,1],1),"%)", sep = "")) + ylab(paste(colnames(mysum$sites)[2]," (",round(100*mysum$cont$importance[2,2],1),"%)", sep = ""))
}
# save the stats
capture.output(mylresultslist, file = paste("output/canonical_plant_season_stats_output.txt", sep = "")) 

capture.output(mypermanova, file = paste("output/permanova_plant_season_stats_output.txt", sep = "")) 

cairo_pdf("output/canonical_plant_season.pdf", height = 5, width = 7)
print(p)
dev.off()






######## Alpha diversity analysis #####

my_otu_tbl <- otumatfin
## prepare the alpha diversity indices individually in a format suitable for long to wide table conversion
# calculate the per sample Shannon index and prep the long table version
Shannon <- data.frame(diversity(t(my_otu_tbl)))
Shannon$idx <- "Shannon"
Shannon$rnames <- row.names(Shannon)
Shannon$samplenms
colnames(Shannon)[1] <- "value"
# calculate the per sample Simpson index and prep the long table version
Simpson <- data.frame(diversity(t(my_otu_tbl), "simpson"))
Simpson$idx <- "Simpson"
Simpson$rnames <- row.names(Simpson)
Simpson$samplenms
colnames(Simpson)[1] <- "value"
# calculate the per sample Fisher's a index and prep the long table version
Fishers_a <- data.frame(fisher.alpha(t(my_otu_tbl)))
Fishers_a$idx <- "Fishers_a"
Fishers_a$rnames <- row.names(Fishers_a)
Fishers_a$samplenms
colnames(Fishers_a)[1] <- "value"
# calculate the per sample richness and prep the long table version
Richness <- data.frame(specnumber(t(my_otu_tbl))) ## rowSums(BCI > 0) does the same...
Richness$idx <- "Richness"
Richness$rnames <- row.names(Richness)
Richness$samplenms
colnames(Richness)[1] <- "value"
# calculate the per sample Pielou's eveness and prep the long table version
Pielous_eveness <- data.frame(diversity(t(my_otu_tbl))/log(specnumber(t(my_otu_tbl))))
Pielous_eveness$idx <- "Pielous_eveness"
Pielous_eveness$rnames <- row.names(Pielous_eveness)
Pielous_eveness$samplenms
colnames(Pielous_eveness)[1] <- "value"

# row-bind the tables, convert the final table form long to wide, and add the necessary design factors for using perwise or ANOVA and alike methods
library(reshape2)
my_setup <- merge(dcast(rbind(Shannon,Simpson,Fishers_a,Richness,Pielous_eveness), rnames ~ idx), design, by.x = "rnames", by.y = "row.names", all = T)
my_setup$season <- factor(my_setup$season,levels = c("Summer","Winter"))
my_setup$plant <- factor(my_setup$plant)

# run the per variable tests and save them in a list object
library(agricolae) 
library(car) 
# create an empty list that will host the results
mystatsout <- list()
# set the to be tested variables
myalphaindices <- c("Shannon","Simpson","Fishers_a","Richness","Pielous_eveness")
# loop through the various tests
for(myalphaindex in myalphaindices){
  # ... and run and save in the list the t.tests comparing seasons per each plant (non-parametric ANOVA equivalent) through a loop
  for(myplant in levels(my_setup$plant)){
    mystatsout[[paste(myalphaindex,"Student's t tests")]][[paste("t.test",myplant)]] <- t.test(my_setup[,myalphaindex][my_setup$plant == myplant & my_setup$season == "Summer"],my_setup$Shannon[my_setup$plant == myplant & my_setup$season == "Winter"])
  }
  # then loop through summer and winter and test with ANOVA the plant differences
  for(myseason in c("Summer","Winter")){
    # ... and run and save in the list the Shapiro test for normality of the used values
    mystatsout[[paste(myalphaindex,"Shapiro",myseason)]] <- shapiro.test(my_setup[my_setup$season == myseason,myalphaindex])
    # ... and run and save in the list the Levene test for homogenity of variances
    mystatsout[[paste(myalphaindex,"Levene",myseason)]] <- leveneTest(as.formula(paste(myalphaindex,"~ plant")), data=my_setup[my_setup$season == myseason,])
    # ... and run and save in the list the ANOVA
    model <- aov(as.formula(paste(myalphaindex,"~ plant")), data = my_setup[my_setup$season == myseason,])
    mystatsout[[paste(myalphaindex,"ANOVA",myseason)]] <- HSD.test(model,"plant", group=TRUE,console=TRUE, main=paste(myalphaindex,"\ plant"))
    # ... and run and save in the list the Kruskal (non-parametric ANOVA equivalent)
    mystatsout[[paste(myalphaindex,"Kruskal",myseason)]] <- kruskal(my_setup[my_setup$season == myseason,myalphaindex], my_setup[my_setup$season == myseason,]$plant, alpha = 0.05, p.adj="fdr", group=TRUE)
  }
}
# save the output list used for statistical significance grouping in Fig. S4 (the associated bar-plot was prepared in excel with the significance letters being added manually)
capture.output(mystatsout, file = "output/alpha_div_stats.txt")




########### Barplots ##################
# the phyloseq, together with other packages listed below were used for producing the barplots
library(phyloseq) # used to produce easily subsetable OTU matrix containing objects
library(RColorBrewer) # necessary for colours
library(PerformanceAnalytics) # necessary for the rainbowXequal command
library(plyr) # necessary for easily renaming factors with the mapvalues command
library(ggplot2) # necessary for plots


### prepare the dataset in the phyloseq format
# set the OTU table
my_otu_tbl <- t(otumat)
# set the design table
my_samp_dat <- design
# set some design table variables used in this plotting section as factors
my_samp_dat$plant <- factor(my_samp_dat$plant)
my_samp_dat$season <- factor(my_samp_dat$season, levels = c("Summer","Winter"))
# set the OTU classification table
my_tax_tbl <- hiera
# select the taxa of interest
my_tax_tbl_2 <- my_tax_tbl[which(my_tax_tbl$Domain == "Archaea"),]
## Fix the annotation
# E.g. one issue with the fungal annotation is that they use unidentified also next to the ? and they add in the end of the unidientified  the lowest classified level along with unidentified XXX sp... I had to replace these annotations with ? in order to give the lowest possible classifications to the barplots 
my_tax_tbl_3 <- my_tax_tbl_2
for(myrow in 1:nrow(my_tax_tbl_3)){
  if(my_tax_tbl_3[myrow,"Genus"] == "unidentified"){
    my_tax_tbl_3[myrow,"Species"] <- "unidentified"
  }
}
## replace the "unidentified" with the unclassified leftmost (higher level annotated) taxa
my_tax_tbl_3b <- my_tax_tbl_3
for(myrow in 1:nrow(my_tax_tbl_3b)){
  for(mycol in 2:ncol(my_tax_tbl_3b)){
    if(my_tax_tbl_3b[myrow,mycol] == "unidentified" & length(grep("uncl\\.",my_tax_tbl_3b[myrow,mycol-1])) > 0){
      my_tax_tbl_3b[myrow,mycol] <- my_tax_tbl_3b[myrow,mycol-1]
    } else if (my_tax_tbl_3b[myrow,mycol] == "unidentified" & length(grep("uncl\\.",my_tax_tbl_3b[myrow,mycol-1])) == 0) {
      my_tax_tbl_3b[myrow,mycol] <- paste("uncl.",my_tax_tbl_3b[myrow,mycol-1])
    } else {
      my_tax_tbl_3b[myrow,mycol] <- my_tax_tbl_3b[myrow,mycol]
    }
  }
}
## replace the ? with the unclassified leftmost taxa
my_tax_tbl_3c <- my_tax_tbl_3b
for(myrow in 1:nrow(my_tax_tbl_3c)){
  for(mycol in 2:ncol(my_tax_tbl_3c)){
    if(my_tax_tbl_3c[myrow,mycol] == "?" & length(grep("uncl\\.",my_tax_tbl_3c[myrow,mycol-1])) > 0){
      my_tax_tbl_3c[myrow,mycol] <- my_tax_tbl_3c[myrow,mycol-1]
    } else if (my_tax_tbl_3c[myrow,mycol] == "?" & length(grep("uncl\\.",my_tax_tbl_3c[myrow,mycol-1])) == 0) {
      my_tax_tbl_3c[myrow,mycol] <- paste("uncl.",my_tax_tbl_3c[myrow,mycol-1])
    } else {
      my_tax_tbl_3c[myrow,mycol] <- my_tax_tbl_3c[myrow,mycol]
    }
  }
}
## replace the Incertae sedis with the taxon in the left plus IS (for easier visualization)
my_tax_tbl_4 <- my_tax_tbl_3c
for(myrow in 1:nrow(my_tax_tbl_4)){
  for(mycol in 2:ncol(my_tax_tbl_4)){
    if(my_tax_tbl_4[myrow,mycol] == "Incertae sedis"){
      my_tax_tbl_4[myrow,mycol] <- paste(my_tax_tbl_4[myrow,mycol-1], " IS", sep = "")
    } else {
      my_tax_tbl_4[myrow,mycol] <- my_tax_tbl_4[myrow,mycol]
    }
  }
}
my_tax_tbl_5 <- my_tax_tbl_4
# replace double IS with a single one coming up in some cases
for(myrow in 1:nrow(my_tax_tbl_5)){
  for(mycol in 2:ncol(my_tax_tbl_5)){
    my_tax_tbl_5[myrow,mycol] <- gsub("IS IS","IS", my_tax_tbl_5[myrow,mycol])
  }
}
# assign the table to a final name
my_tax_tbl_fin <- my_tax_tbl_5

# match OTUs between my_otu_tbl and my_tax_tbl
my_otu_tbl_1 <- my_otu_tbl[,which(row.names(t(my_otu_tbl))%in%row.names(my_tax_tbl_fin))]
my_otu_tbl_2 <- my_otu_tbl_1 # save to new table for security

# check the OTU table for samples with 0 sequences due to the removal of OTUs and remove these samples
if(length(which(rowSums(my_otu_tbl_2) == 0)) > 0){
  my_otu_tbl_3 <- my_otu_tbl_2[-which(rowSums(my_otu_tbl_2) == 0),]
} else {
  my_otu_tbl_3 <- my_otu_tbl_2
}

# assign as matrix in order to be able to later on import the file as phyloseq object
my_otu_tbl_fin <- data.matrix(my_otu_tbl_3)
# Make sure that samples match between the OTU matrix and the sample data
my_samp_dat_fin <- my_samp_dat[which(row.names(my_samp_dat)%in%row.names(my_otu_tbl_fin)),]

# create a phyloseq object (keep in mind that the tax_table command accepts matrices and not data.frames)
my_phydata <- phyloseq(otu_table(my_otu_tbl_fin, taxa_are_rows = FALSE), tax_table(as.matrix(my_tax_tbl_fin)), sample_data(my_samp_dat_fin))



### prepare a barplot for the relative abundance data
# convert to relative abundances
my_phydata_ra <- transform_sample_counts(my_phydata, function(x) {x*100/sum(x)})

# get the taxonomy levels from the column names of the phyloseq table
my_tax_lvls <- colnames(my_phydata_ra@tax_table@.Data)

# create the directory where to save the generated files
system("mkdir output/barplots")

# iterate through the taxa and produce the taxon-level barplots
for(my_tax_lev in my_tax_lvls[3:length(my_tax_lvls)]){
  
  # modify the relative abundance object by merging the samples according to the "planttreat" parameter and dividing by the number of replicates (since the mean function seems not to work in the case of merging for some reason... check if this changes in the future... in case this happens add the fun parameter and omit the line after the following)
  my_phydata_rafr_plot_bar <- merge_samples(t(my_phydata_ra), "planttreat") 
  otu_table(my_phydata_rafr_plot_bar) <- otu_table(data.frame(my_phydata_rafr_plot_bar@otu_table, check.names = F)/(rowSums(data.frame(my_phydata_rafr_plot_bar@otu_table, check.names = F))/100), taxa_are_rows = F)

  ## amend the plot_bar function by setting the legend order according to the relative abundances of events, add the number of plotted features, and remove the black slices feature from plot_bar as indicated in the hashed text of the plot_bar function version below
  plot_bar2 <- function (physeq, x = "Sample", y = "Abundance", fill = NULL, toptaxa = NULL, 
                         title = NULL, facet_grid = NULL) 
  {
    mdf = psmelt(physeq)
    # set the legend order according to the dataset-wise relative abundance of taxa by inserting the two lines below
    myorddf <- aggregate(mdf[,"Abundance"], by = list(mdf[,fill]), sum)
    mdf[,fill] <- factor(mdf[,fill], levels = myorddf$Group.1[order(myorddf$x, decreasing = T)])
    # set the toptaxa number to plot (e.g. according to the colours) by labeling the rest as "Other" in case the selected taxon number is lower than those existing in the dataset
    if(length(levels(mdf[,fill]))-toptaxa > 0){
      mdf[,fill] <- plyr::mapvalues(mdf[,fill], from = levels(mdf[,fill]), to = c(levels(mdf[,fill])[1:toptaxa], rep("Others", length(levels(mdf[,fill]))-toptaxa))) 
    }
    
    p = ggplot(mdf, aes_string(x = x, y = y, fill = fill))
    # remove the ", color = "black"" parameter from the line below
    # p = p + geom_bar(stat = "identity", position = "stack", color = "black")
    p = p + geom_bar(stat = "identity", position = "stack")
    p = p + theme(axis.text.x = element_text(angle = -90, hjust = 0))
    if (!is.null(facet_grid)) {
      p <- p + facet_grid(facet_grid)
    }
    if (!is.null(title)) {
      p <- p + ggtitle(title)
    }
    return(p)
  }
  
  
  
  # use up to 11 colours for labeling the taxa
  mycols <- c(rainbow10equal,"grey70")
  # run the plotting function
  p <- plot_bar2(my_phydata_rafr_plot_bar, fill = "Class", toptaxa = 11)
  # replace the underscores with spaces in the sample names
  p$data$Sample <- gsub("_"," ",p$data$Sample)
  # set the black and white theme of the plot background
  theme_set(theme_bw())
  # save the plot using the cairo graphics device (Greek-letter friendly)
  cairo_pdf(height = 4, width = 8, file = paste("output/barplots/",my_tax_lev,"_means.pdf", sep = ""))
  # print the plot to the device and amend the colours to those of choice, add the appropriate axis labels, rotate the sample labels
  print(p + scale_fill_manual(values=mycols) + ylab("Relative abundance (%)") + xlab("") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
  # close the graphics device
  dev.off()
  
  ## prepare also the facet wrap version of the above plot
  # first add the season info in the sample_data table which for some reason was converted into numeric values (possibly during the merging)
  mysmpldata <- data.frame(sample_data(my_phydata_rafr_plot_bar))
  mysmpldata$season <- gsub(".+_","",row.names(mysmpldata))
  sample_data(my_phydata_rafr_plot_bar) <- sample_data(mysmpldata)
  # the following facet_wrap differentiates the two plots
  p <- plot_bar2(my_phydata_rafr_plot_bar, fill = "Class", toptaxa = 11) + facet_wrap(~season, scales="free_x") + theme(strip.background =element_rect(fill ="white"))
  p$data$Sample <- gsub(" Summer| Winter","",gsub("_"," ",p$data$Sample))
  theme_set(theme_bw())
  cairo_pdf(height = 4, width = 8, file = paste("output/barplots/",my_tax_lev,"_facet.pdf", sep = ""))
  print(p + scale_fill_manual(values=mycols) + ylab("Relative abundance (%)") + xlab("") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
  dev.off()
}






########### Heatmaps ##################


library(vegan)

### prepare the dataset in the phyloseq format

my_otu_tbl <- t(otumat)
my_samp_dat <- design
my_samp_dat <- data.frame(my_samp_dat,plant= interaction(my_samp_dat$plant,my_samp_dat$season))
my_samp_dat$plant <- factor(my_samp_dat$plant, levels = c("L._stoechas","C._incanus","A._unedo","P._latifolia","P._lentiscus","M._communis","Q._coccifera","M._officinalis"))
my_samp_dat$treat <- factor(my_samp_dat$season, levels = c("Summer","Winter"))
my_tax_tbl <- hiera
# select the taxa of interest
my_tax_tbl_2 <- my_tax_tbl[which(my_tax_tbl$Domain == "Archaea"),]
## Fix the annotation
# E.g. one issue with the fungal annotation is that they use unidentified also next to the ? and they add in the end of the unidientified  the lowest classified level along with unidentified XXX sp... I had to replace these annotations with ? in order to give the lowest possible classifications to the barplots 
my_tax_tbl_3 <- my_tax_tbl_2
for(myrow in 1:nrow(my_tax_tbl_3)){
  if(my_tax_tbl_3[myrow,"Genus"] == "unidentified"){
    my_tax_tbl_3[myrow,"Species"] <- "unidentified"
  }
}
## replace the "unidentified" with the unclassified leftmost (higher level annotated) taxa
my_tax_tbl_3b <- my_tax_tbl_3
for(myrow in 1:nrow(my_tax_tbl_3b)){
  for(mycol in 2:ncol(my_tax_tbl_3b)){
    if(my_tax_tbl_3b[myrow,mycol] == "unidentified" & length(grep("uncl\\.",my_tax_tbl_3b[myrow,mycol-1])) > 0){
      my_tax_tbl_3b[myrow,mycol] <- my_tax_tbl_3b[myrow,mycol-1]
    } else if (my_tax_tbl_3b[myrow,mycol] == "unidentified" & length(grep("uncl\\.",my_tax_tbl_3b[myrow,mycol-1])) == 0) {
      my_tax_tbl_3b[myrow,mycol] <- paste("uncl.",my_tax_tbl_3b[myrow,mycol-1])
    } else {
      my_tax_tbl_3b[myrow,mycol] <- my_tax_tbl_3b[myrow,mycol]
    }
  }
}
## replace the ? with the unclassified leftmost taxa
my_tax_tbl_3c <- my_tax_tbl_3b
for(myrow in 1:nrow(my_tax_tbl_3c)){
  for(mycol in 2:ncol(my_tax_tbl_3c)){
    if(my_tax_tbl_3c[myrow,mycol] == "?" & length(grep("uncl\\.",my_tax_tbl_3c[myrow,mycol-1])) > 0){
      my_tax_tbl_3c[myrow,mycol] <- my_tax_tbl_3c[myrow,mycol-1]
    } else if (my_tax_tbl_3c[myrow,mycol] == "?" & length(grep("uncl\\.",my_tax_tbl_3c[myrow,mycol-1])) == 0) {
      my_tax_tbl_3c[myrow,mycol] <- paste("uncl.",my_tax_tbl_3c[myrow,mycol-1])
    } else {
      my_tax_tbl_3c[myrow,mycol] <- my_tax_tbl_3c[myrow,mycol]
    }
  }
}
## replace the Incertae sedis with the taxon in the left plus IS (for easier visualization)
my_tax_tbl_4 <- my_tax_tbl_3c
for(myrow in 1:nrow(my_tax_tbl_4)){
  for(mycol in 2:ncol(my_tax_tbl_4)){
    if(my_tax_tbl_4[myrow,mycol] == "Incertae sedis"){
      my_tax_tbl_4[myrow,mycol] <- paste(my_tax_tbl_4[myrow,mycol-1], " IS", sep = "")
    } else {
      my_tax_tbl_4[myrow,mycol] <- my_tax_tbl_4[myrow,mycol]
    }
  }
}
my_tax_tbl_5 <- my_tax_tbl_4
# replace double IS with a single one coming up in some cases
for(myrow in 1:nrow(my_tax_tbl_5)){
  for(mycol in 2:ncol(my_tax_tbl_5)){
    my_tax_tbl_5[myrow,mycol] <- gsub("IS IS","IS", my_tax_tbl_5[myrow,mycol])
  }
}
# assign the table to a final name
my_tax_tbl_fin <- my_tax_tbl_5

### prepare the relative abundance tables 
my_otu_tbl_1 <- t(my_otu_tbl[,which(colnames(my_otu_tbl)%in%row.names(my_tax_tbl_fin))])
my_otu_tbl_2_ra <- decostand(my_otu_tbl_1, MARGIN = 2, method = "total")
# set a percent of lowest max relative participation of OTU for considering in plotting and remove OTUs with lower such values
for(myreqrelabun in c(0.02,0.05)){
  my_otu_tbl_2_ra_min_SETmax <- my_otu_tbl_2_ra[which(apply(my_otu_tbl_2_ra,1,max) >= myreqrelabun),]
  # merge the OTU table with the taxonomy table
  my_otu_tbl_2_ra_min_SETmax_and_tax <- merge(my_tax_tbl_fin,my_otu_tbl_2_ra_min_SETmax, by = "row.names",all.y = T)
  row.names(my_otu_tbl_2_ra_min_SETmax_and_tax) <- paste(my_otu_tbl_2_ra_min_SETmax_and_tax$Row.names,my_otu_tbl_2_ra_min_SETmax_and_tax$Genus)
  # select the numeric columns for aggregation etc
  mynums <- unlist(lapply(my_otu_tbl_2_ra_min_SETmax_and_tax, is.numeric))
  mat_pre <- my_otu_tbl_2_ra_min_SETmax_and_tax[,mynums]
  
  
  ### prepare the heatmaps
  # prep the colours and breaks of the colour key
  library(gplots)
  mybreaks <- log10(0.1+c(seq(0,1,length.out = 40),seq(3,10,length.out = 20),seq(11,25,length.out = 5),seq(26,max(mat),length.out = 5)))
  mycolors <- colorRampPalette(colors = c("steelblue4","grey60","grey80","white","orange","orangered"))(n=69)
  
  ## per plant/season mean calculation and plotting of values
  mat_pre_aggr_t <- aggregate(t(mat_pre), by = list(my_samp_dat[row.names(t(mat_pre)),]$planttreat), FUN = mean)
  row.names(mat_pre_aggr_t) <- mat_pre_aggr_t$Group.1
  mat_pre_aggr <- mat_pre_aggr_t[,-which(colnames(mat_pre_aggr_t)%in%"Group.1")]
  mat <- 100*t(mat_pre_aggr)
  
  # plot and save the plot
  library(pheatmap)
  cairo_pdf(height = 14, width = 12, file= paste("output/Pheatmap_per_plant_season_aggr_",myreqrelabun,".pdf", sep = ""))
  pheatmap(as.matrix(log10(mat+0.1)), trace = "none", Rowv = F, Colv = F, breaks = mybreaks , col = mycolors, clustering_method = "ward.D", cluster_rows = F, cluster_cols = T, cellwidth = 12, cellheight = 12, cutree_cols = 2, angle_col=90, display_numbers = matrix(ifelse(mat > 5, "*", ""), nrow(mat)), legend_breaks = c(-1,0,0.30103,0.69897,1,1.8), legend_labels = c("0","1","2","5","10","63"))
  dev.off()
  
  ## per plant mean calculation and plotting of values
  mat_pre_aggr_t <- aggregate(t(mat_pre), by = list(my_samp_dat[row.names(t(mat_pre)),]$plant), FUN = mean)
  row.names(mat_pre_aggr_t) <- mat_pre_aggr_t$Group.1
  mat_pre_aggr <- mat_pre_aggr_t[,-which(colnames(mat_pre_aggr_t)%in%"Group.1")]
  mat <- 100*t(mat_pre_aggr)
  
  # plot and save the plot
  library(pheatmap)
  cairo_pdf(height = 14, width = 12, file= paste("output/Pheatmap_per_plant_aggr_",myreqrelabun,".pdf", sep = ""))
  pheatmap(as.matrix(log10(mat+0.1)), trace = "none", Rowv = F, Colv = F, breaks = mybreaks , col = mycolors, clustering_method = "ward.D", cluster_rows = F, cluster_cols = T, cellwidth = 12, cellheight = 12, cutree_cols = 2, angle_col=90, display_numbers = matrix(ifelse(mat > 5, "*", ""), nrow(mat)), legend_breaks = c(-1,0,0.30103,0.69897,1,1.8), legend_labels = c("0","1","2","5","10","63"))
  dev.off()
  
  ## per plant aromatic character mean calculation and plotting of values
  mat_pre_aggr_t <- aggregate(t(mat_pre), by = list(my_samp_dat[row.names(t(mat_pre)),]$Aromatic), FUN = mean)
  row.names(mat_pre_aggr_t) <- mat_pre_aggr_t$Group.1
  mat_pre_aggr <- mat_pre_aggr_t[,-which(colnames(mat_pre_aggr_t)%in%"Group.1")]
  mat <- 100*t(mat_pre_aggr)
  
  # plot and save the plot
  library(pheatmap)
  cairo_pdf(height = 14, width = 12, file= paste("output/Pheatmap_per_plant_arom_char_aggr_",myreqrelabun,".pdf",sep = ""))
  pheatmap(as.matrix(log10(mat+0.1)), trace = "none", Rowv = F, Colv = F, breaks = mybreaks , col = mycolors, clustering_method = "ward.D", cluster_rows = F, cluster_cols = T, cellwidth = 12, cellheight = 12, cutree_cols = 2, angle_col=90, display_numbers = matrix(ifelse(mat > 5, "*", ""), nrow(mat)), legend_breaks = c(-1,0,0.30103,0.69897,1,1.8), legend_labels = c("0","1","2","5","10","63"))
  dev.off()
  
  ## per plant canopy type mean calculation and plotting of values
  mat_pre_aggr_t <- aggregate(t(mat_pre), by = list(my_samp_dat[row.names(t(mat_pre)),]$type_canopy), FUN = mean)
  row.names(mat_pre_aggr_t) <- mat_pre_aggr_t$Group.1
  mat_pre_aggr <- mat_pre_aggr_t[,-which(colnames(mat_pre_aggr_t)%in%"Group.1")]
  mat <- 100*t(mat_pre_aggr)
  
  # plot and save the plot
  library(pheatmap)
  cairo_pdf(height = 14, width = 12, file= paste("output/Pheatmap_per_plant_canopy_type_aggr_",myreqrelabun,".pdf", sep = ""))
  pheatmap(as.matrix(log10(mat+0.1)), trace = "none", Rowv = F, Colv = F, breaks = mybreaks , col = mycolors, clustering_method = "ward.D", cluster_rows = F, cluster_cols = T, cellwidth = 12, cellheight = 12, cutree_cols = 2, angle_col=90, display_numbers = matrix(ifelse(mat > 5, "*", ""), nrow(mat)), legend_breaks = c(-1,0,0.30103,0.69897,1,1.8), legend_labels = c("0","1","2","5","10","63"))
  dev.off()
}
