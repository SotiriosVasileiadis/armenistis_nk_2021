library(agricolae)
library(car)
my_setup <- read.table("qpcr_g_leaf.txt", header = T, row.names = 1, check.names = F, sep  ="\t", quote = "", comment.char = "")
my_setup$Plant <- factor(my_setup$Plant)
my_setup$Aromatic <- factor(my_setup$Aromatic)
my_setup$Type_plant <- factor(my_setup$Type_plant)
my_setup$Season <- factor(my_setup$Season)


#### Figure 1 ----
# create an empty list that will host the results
mystatsout <- list()
# set the to be tested dependent variables
testvars <- c("T.Fungi","T.Bacteria","Crenarchaea","Cladosporium","Alternaria")

# loop through the various tests
for(testvar in testvars){
  # ... and run and save in the list the t.tests comparing seasons per each plant through a loop
  for(myplant in levels(my_setup$Plant)){
    mystatsout[[paste(testvar,"Student's t tests")]][[paste("t.test",myplant,"for Sumer vs Winter group comparison")]] <- t.test(my_setup[my_setup$Plant == myplant & my_setup$Season == "Summer",testvar],my_setup$Shannon[my_setup$Plant == myplant & my_setup$Season == "Winter"])
  }
  # then loop through summer and winter and test with ANOVA the plant differences
  for(myseason in c("Summer","Winter")){
    # ... and run and save in the list the Shapiro test for normality of the used values
    mystatsout[[paste(testvar,"Shapiro",myseason,"samples for between plant species comparison")]] <- shapiro.test(my_setup[my_setup$Season == myseason,testvar])
    # ... and run and save in the list the Levene test for homogenity of variances
    mystatsout[[paste(testvar,"Levene",myseason,"samples for between plant species comparison")]] <- leveneTest(as.formula(paste(testvar,"~ Plant")), data=my_setup[my_setup$Season == myseason,])
    # ... and run and save in the list the ANOVA
    model <- aov(as.formula(paste(testvar,"~ Plant")), data = my_setup[my_setup$Season == myseason,])
    mystatsout[[paste(testvar,"ANOVA",myseason,"samples for between plant species comparison")]] <- HSD.test(model,"Plant", group=TRUE,console=TRUE, main=paste(testvar,"\ Plant"))
    # ... and run and save in the list the Kruskal (non-parametric ANOVA equivalent)
    mystatsout[[paste(testvar,"Kruskal",myseason,"samples for between plant species comparison")]] <- kruskal(my_setup[my_setup$Season == myseason,testvar], my_setup[my_setup$Season == myseason,]$Plant, alpha = 0.05, p.adj="fdr", group=TRUE)
  }
}
# save the output list used for statistical significance grouping in Fig. 1 (the associated bar-plot was prepared in excel with the significance letters being added manually)
capture.output(mystatsout, file = "Figure_1_stats.txt")





#### Figure S2 ----
# create an empty list that will host the results
mystatsout <- list()
# set the to be tested dependent variables
testvars <- c("T.Fungi","T.Bacteria","Crenarchaea","Cladosporium","Alternaria")

# loop through the various tests
for(testvar in testvars){
  # ... and run and save in the list the t.tests comparing seasons per each aromatic character through a loop
  for(aromaticchar in levels(my_setup$Aromatic)){
    mystatsout[[paste(testvar,"Student's t tests")]][[paste("t.test for",aromaticchar,"comparing the Summer with the Winter group")]] <- t.test(my_setup[my_setup$Aromatic == aromaticchar & my_setup$Season == "Summer",testvar],my_setup$Shannon[my_setup$Aromatic == aromaticchar & my_setup$Season == "Winter"])
  }
  # ... and run and save in the list the t.tests comparing aromatic characters per each season through a loop
  for(myseason in levels(my_setup$Season)){
    mystatsout[[paste(testvar,"Student's t tests")]][[paste("t.test for",myseason,"comparing the Aromatic with the non-aromatic group")]] <- t.test(my_setup[my_setup$Season == myseason & my_setup$Aromatic == "Aromatic",testvar],my_setup$Shannon[my_setup$Season == myseason & my_setup$Aromatic == "No_aromatic"])
  }
  # ... run the same tests for the canopy type instead of aromatic character
  for(typeplant in levels(my_setup$Type_plant)){
    mystatsout[[paste(testvar,"Student's t tests")]][[paste("t.test for",typeplant,"comparing the Summer with the Winter group")]] <- t.test(my_setup[my_setup$Type_plant == typeplant & my_setup$Season == "Summer",testvar],my_setup$Shannon[my_setup$Type_plant == typeplant & my_setup$Season == "Winter"])
  }
  for(myseason in levels(my_setup$Season)){
    mystatsout[[paste(testvar,"Student's t tests")]][[paste("t.test for",myseason,"comparing the High with the low canopy group")]] <- t.test(my_setup[my_setup$Season == myseason & my_setup$Type_plant == "Higher_canopy",testvar],my_setup$Shannon[my_setup$Season == myseason & my_setup$Type_plant == "Lower_canopy"])
  }
}
# save the output list used for statistical significance grouping in Fig. S2 (the associated bar-plot was prepared in excel with the significance letters being added manually)
capture.output(mystatsout, file = "Figure_S2_stats.txt")
