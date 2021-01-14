library(agricolae)
library(car)
my_setup <- read.table("qpcr_g_leaf.txt", header = T, row.names = 1, check.names = F, sep  ="\t", quote = "", comment.char = "")


# create an empty list that will host the results
mystatsout <- list()
# set the to be tested dependent variables
testvars <- c("T.Fungi","T.Bacteria","Crenarchaea","Cladosporium","Alternaria")
# set the independent variables
indpvars <- c("Type_plant_season","Aromatic_season","Plant","Season","planttreat","Type_plant","Aromatic")
# loop through the various tests
for (indpvar in indpvars){
  for(testvar in testvars){
    # ... and run and save in the list the Shapiro test for normality of the used values
    mystatsout[[paste(indpvar,testvar,"Shapiro")]] <- shapiro.test(my_setup[,testvar])
    # ... and run and save in the list the Levene test for homogenity of variances
    mystatsout[[paste(indpvar,testvar,"Levene")]] <- leveneTest(as.formula(paste(testvar,"~",indpvar)), data=my_setup)
    # ... and run and save in the list the ANOVA
    model <- aov(as.formula(paste(testvar,"~",indpvar)), data = my_setup)
    mystatsout[[paste(indpvar,testvar,"ANOVA")]] <- HSD.test(model,indpvar, group=TRUE,console=TRUE, main=paste(testvar,"\\",indpvar))
    # ... and run and save in the list the Kruskal (non-parametric ANOVA equivalent)
    mystatsout[[paste(indpvar,testvar,"Kruskal")]] <- kruskal(my_setup[,testvar], my_setup[,indpvar], alpha = 0.05, p.adj="fdr", group=TRUE)

  }
  
}
# save the outputs into a list (the associated bar-plot was prepared in excel with the significance letters being added manually)
capture.output(mystatsout, file = "qpcr_stats_out.txt")

#### figure 1
#### select according to season and test per plant

#### pairise comparisons according to season within each plant


#### figure S2
#### test within each season the aromatic character 

#### test seasons within each aromatic character group