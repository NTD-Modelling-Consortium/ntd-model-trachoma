library(readr)
library("RColorBrewer")
source("analyze_ihme_ipm_data.R")

# set the path for where the simulations are stored
path_header <- '/Users/user/Documents/'


# full set of IUs
trachomaIUs <- read.csv("trachomaIUs.csv")
# select subset of IUs to test plot with
which_IUs <- c("BDI06385","BDI06386","BDI06388","BDI06392",
               "BEN03235","BEN03236","BEN03237","BEN03238","BEN03239",
               "COD14320","COD14321","COD14324","COD14330")
trachomaIUs <- trachomaIUs[which(trachomaIUs$IU_ID2 %in% which_IUs),]  

# do scenario 3's
count_1 = 0
count_2a = 0
count_2b = 0
count_2c = 0




for(i in 1:nrow(trachomaIUs)){

  
  ihme_name = paste0(path_header,'endgame-ihme-ipm-outputs-trachoma-export-20220908a/ihme-',trachomaIUs$IU_ID2[i],"-trachoma-scenario_1-200_simulations.csv")
  if(file.exists(ihme_name)){
    ihme1 = read.csv(ihme_name)

    if(count_1 == 0){
      
      prop_inf_1 = matrix(0, nrow(trachomaIUs),  length(get_data_over_years_trachoma(ihme1, 2020:2040)/1000))
      prop_inf_2a = prop_inf_1
      prop_inf_2b = prop_inf_1
      prop_inf_2c = prop_inf_1
      count_1 = 1
                          
    }
    prop_inf_1[i,] = get_data_over_years_trachoma(ihme1, 2020:2040)/1000
    
  }
  
  
  ihme_name = paste0(path_header,'endgame-ihme-ipm-outputs-trachoma-export-20220908a/ihme-',trachomaIUs$IU_ID2[i],"-trachoma-scenario_2a-200_simulations.csv")
  if(file.exists(ihme_name)){
    ihme1 = read.csv(ihme_name)
    prop_inf_2a[i,] = get_data_over_years_trachoma(ihme1, 2020:2040)/1000
    
  }
  
  ihme_name = paste0(path_header,'endgame-ihme-ipm-outputs-trachoma-export-20220908a/ihme-',trachomaIUs$IU_ID2[i],"-trachoma-scenario_2b-200_simulations.csv")
  if(file.exists(ihme_name)){
    ihme1 = read.csv(ihme_name)
    prop_inf_2b[i,] = get_data_over_years_trachoma(ihme1, 2020:2040)/1000
    
  }
  
  ihme_name = paste0(path_header,'endgame-ihme-ipm-outputs-trachoma-export-20220908a/ihme-',trachomaIUs$IU_ID2[i],"-trachoma-scenario_2c-200_simulations.csv")
  if(file.exists(ihme_name)){
    ihme1 = read.csv(ihme_name)
    prop_inf_2c[i,] = get_data_over_years_trachoma(ihme1, 2020:2040)/1000
    
  }

  
  print(paste("Done", i,"of", nrow(trachomaIUs)))
  }
  
  
############ code runs to here no errors ############ 
# infer max age from files
age_groups = 0:max(ihme1$age_end)

# group_pop <- 300000
# trachomaIUs_scen3 <- data.frame(pop = rep(1000, nrow(trachomaIUs)))
group_pop = sum(trachomaIUs_scen3$pop)

total_infs1 = trachomaIUs_scen3$pop*prop_inf_1
total_infs2a = trachomaIUs_scen3$pop*prop_inf_2a
total_infs2b = trachomaIUs_scen3$pop*prop_inf_2b
total_infs2c = trachomaIUs_scen3$pop*prop_inf_2c

total_infs1 = colSums(total_infs1/group_pop)
total_infs2a = colSums(total_infs2a/group_pop)
total_infs2b = colSums(total_infs2b/group_pop)
total_infs2c = colSums(total_infs2c/group_pop)


cols = c("#0098FF", "#1b9e77", "#d95f02","#7570b3", '#d4d133')
png("trachoma_Group_trajectory.png", height = 8, width = 12, res = 300, units = "in")
plot(2020:2040, total_infs1, type = 'l', lwd = 4, 
     xlab = 'year', ylab = 'prevalence all infections',
     col = cols[1],
     bty = 'n', ylim = c(0, max(total_infs1)),
     
     cex = 1.7, cex.axis = 1.7, cex.lab = 1.7,cex.main = 1.7,
     main = 'trachoma Group 1')


lines(2020:2040,total_infs2a, type = 'l', lwd = 4, col = cols[2])

lines(2020:2040,total_infs2b, type = 'l', lwd = 4, col = cols[3])

lines(2020:2040,total_infs2c, type = 'l', lwd = 4, 
      col = cols[4],
      bty = 'n')


legend('topright',  legend = c("Scenario 1", "Scenario 2a", "Scenario 2b", "Scenario 2c"),
       col = cols, lwd = c(4,4,4,4),
       bty = 'n', cex = 1.7)
dev.off()


###############################################################################


group1_pop = sum(trachomaIUs_scen3$pop[which(trachomaIUs_scen3$CoverageGroup == 1)])
group2_pop = sum(trachomaIUs_scen3$pop[which(trachomaIUs_scen3$CoverageGroup == 2)])

total_infs0_1 = trachomaIUs_scen3$pop*prop_inf_0_1
total_infs1_1 = trachomaIUs_scen3$pop*prop_inf_1_1
total_infs2_1 = trachomaIUs_scen3$pop*prop_inf_2_1
total_infs3a_1 = trachomaIUs_scen3$pop*prop_inf_3a_1
total_infs3b_1 = trachomaIUs_scen3$pop*prop_inf_3b_1

total_infs0_1 = colSums(total_infs0_1/group1_pop)
total_infs1_1 = colSums(total_infs1_1/group1_pop)
total_infs2_1 = colSums(total_infs2_1/group1_pop)
total_infs3a_1 = colSums(total_infs3a_1/group1_pop)
total_infs3b_1 = colSums(total_infs3b_1/group1_pop)

total_infs0_2 = trachomaIUs_scen3$pop*prop_inf_0_2
total_infs1_2 = trachomaIUs_scen3$pop*prop_inf_1_2
total_infs2_2 = trachomaIUs_scen3$pop*prop_inf_2_2
total_infs3a_2 = trachomaIUs_scen3$pop*prop_inf_3a_2
total_infs3b_2 = trachomaIUs_scen3$pop*prop_inf_3b_2

total_infs0_2 = colSums(total_infs0_2/group2_pop)
total_infs1_2 = colSums(total_infs1_2/group2_pop)
total_infs2_2 = colSums(total_infs2_2/group2_pop)
total_infs3a_2 = colSums(total_infs3a_2/group2_pop)
total_infs3b_2 = colSums(total_infs3b_2/group2_pop)


cols = c("#0098FF", "#1b9e77", "#d95f02","#7570b3", '#d4d133')
png("trachoma_Group1_trajectory.png", height = 8, width = 12, res = 300, units = "in")
plot(2020:2040, colSums(prop_inf_1), type = 'l', lwd = 4, 
     xlab = 'year', ylab = 'prevalence all infections',
     col = cols[1],
     bty = 'n', ylim = c(0, max(colSums(prop_inf_1))),
     
     cex = 1.7, cex.axis = 1.7, cex.lab = 1.7,cex.main = 1.7,
     main = 'trachoma Group 1')


lines(2020:2040,colSums(prop_inf_2a), type = 'l', lwd = 4, col = cols[2])

lines(2020:2040,colSums(prop_inf_2b), type = 'l', lwd = 4, col = cols[3])


lines(2020:2040,colSums(prop_inf_2c), type = 'l', lwd = 4, 
      col = cols[4],
      bty = 'n')


legend('topright',  legend = c("Scenario 1", "Scenario 2a", "Scenario 2b", "Scenario 2c"),
       col = cols, lwd = c(4,4,4,4),
       bty = 'n', cex = 1.7)
dev.off()

png("trachoma_Group2_trajectory.png", height = 8, width = 12, res = 300, units = "in")
plot(2020:2040, total_infs1_2, type = 'l', lwd = 4, 
     xlab = 'year', ylab = 'prevalence all infections',
     col = cols[1],
     bty = 'n', ylim = c(0, max(total_infs1_2)),
     
     cex = 1.7, cex.axis = 1.7, cex.lab = 1.7,cex.main = 1.7,
     main = 'trachoma Group 2')

lines(2020:2040,total_infs2_2, type = 'l', lwd = 4, col = cols[2])


lines(2020:2040,total_infs3a_2, type = 'l', lwd = 4, 
      col = cols[3],
      bty = 'n')
lines(2020:2040,total_infs3b_2, type = 'l', lwd = 4, 
      col = cols[4],
      bty = 'n')

legend('topright',  legend = c("Scenario 1", "Scenario 2", "Scenario 3a","Scenario 3b"),
       col = cols, lwd = c(4,4,4,4),
       bty = 'n', cex = 1.7)
dev.off()

total_infs1_1 = trachomaIUs_scen3$pop*prop_inf_1_1
total_infs2_1 = trachomaIUs_scen3$pop*prop_inf_2_1
total_infs3a_1 = trachomaIUs_scen3$pop*prop_inf_3a_1
total_infs3b_1 = trachomaIUs_scen3$pop*prop_inf_3b_1
total_infs1_2 = trachomaIUs_scen3$pop*prop_inf_1_2
total_infs2_2 = trachomaIUs_scen3$pop*prop_inf_2_2
total_infs3a_2 = trachomaIUs_scen3$pop*prop_inf_3a_2
total_infs3b_2 = trachomaIUs_scen3$pop*prop_inf_3b_2


all_1 = colSums((total_infs1_1 + total_infs1_2)/sum(trachomaIUs_scen3$pop))
all_2 = colSums((total_infs2_1 + total_infs2_2)/sum(trachomaIUs_scen3$pop))
all_3a = colSums((total_infs3a_1 + total_infs3a_2)/sum(trachomaIUs_scen3$pop))
all_3b = colSums((total_infs3b_1 + total_infs3b_2)/sum(trachomaIUs_scen3$pop))




png("trachoma_All_trajectory.png", height = 8, width = 12, res = 300, units = "in")
plot(2020:2040, all_1, type = 'l', lwd = 4, 
     xlab = 'year', ylab = 'prevalence all infections',
     col = cols[1],
     bty = 'n', ylim = c(0, max(all_1)),
     
     cex = 1.7, cex.axis = 1.7, cex.lab = 1.7,cex.main = 1.7,
     main = 'trachoma')

lines(2020:2040,all_2, type = 'l', lwd = 4, col = cols[2])


lines(2020:2040,all_3a, type = 'l', lwd = 4, 
      col = cols[3],
      bty = 'n')
lines(2020:2040,all_3b, type = 'l', lwd = 4, 
      col = cols[4],
      bty = 'n')

legend('topright',  legend = c("Scenario 1", "Scenario 2", "Scenario 3a","Scenario 3b"),
       col = cols, lwd = c(4,4,4,4),
       bty = 'n', cex = 1.7)
dev.off()

prop_infected_1_1_tri = totalNumInfecteds_1_1/(3000*count1_1)
prop_infected_2_1_tri = totalNumInfecteds_2_1/(3000*count2_1)
prop_infected_3a1_tri = totalNumInfecteds_3a1/(3000*count3a_1)
prop_infected_3b1_tri = totalNumInfecteds_3b1/(3000*count3b_1)
cols = c("#F8766D", "#00BA38", "#619CFF","#ca6df8")
#prop_infeced = totalNumInfecteds/3000
# png("trichuriasis_trajectory.png", units = "in",
#     width = 12, height = 8, res = 300)
plot(2020:2040,rowMeans(prop_infected_1_1_tri), type = 'l', lwd = 4, 
     xlab = 'year', ylab = 'prevalence all infections',
     col = cols[1],
     bty = 'n', 
     ylim = c(0, max(rowMeans(prop_infected_1_1_tri),
                     rowMeans(prop_infected_2_1_tri),
                     rowMeans(prop_infected_3a1_tri))),
     cex = 1.7, cex.axis = 1.7, cex.lab = 1.7)
# 
# lines(2020:2040,rowMeans(prop_infected_2_tri), type = 'l', lwd = 4, col = cols[2])


lines(2020:2040,rowMeans(prop_infected_2_1_tri), type = 'l', lwd = 4, 
      col = cols[2],
      bty = 'n')
lines(2020:2040,rowMeans(prop_infected_3a1_tri), type = 'l', lwd = 4, 
      col = cols[3],
      bty = 'n')

lines(2020:2040,rowMeans(prop_infected_3b1_tri), type = 'l', lwd = 4, 
      col = cols[4],
      bty = 'n')
legend('topright', title = "", legend = c("Scenario 1_1", "Scenario 2_1", "Scenario 3a_1","Scenario 3b_1"),
       col = cols, lwd = c(4,4,4,4),
       bty = 'n', cex = 1.7)
# dev.off()






prop_infected_1_2_tri = totalNumInfecteds_1_2/(3000*count1_2)
prop_infected_2_2_tri = totalNumInfecteds_2_2/(3000*count2_2)
prop_infected_3a2_tri = totalNumInfecteds_3a2/(3000*count3a_2)
prop_infected_3b2_tri = totalNumInfecteds_3b2/(3000*count3b_2)
cols = c("#F8766D", "#00BA38", "#619CFF","#ca6df8")
#prop_infeced = totalNumInfecteds/3000
# png("trichuriasis_trajectory.png", units = "in",
#     width = 12, height = 8, res = 300)
plot(2020:2040,rowMeans(prop_infected_1_2_tri), type = 'l', lwd = 4, 
     xlab = 'year', ylab = 'prevalence all infections',
     col = cols[1],
     bty = 'n', 
     ylim = c(0, max(rowMeans(prop_infected_1_2_tri),
                     rowMeans(prop_infected_2_2_tri),
                     rowMeans(prop_infected_3a2_tri))),
     cex = 1.7, cex.axis = 1.7, cex.lab = 1.7)
# 
# lines(2020:2040,rowMeans(prop_infected_2_tri), type = 'l', lwd = 4, col = cols[2])


lines(2020:2040,rowMeans(prop_infected_2_2_tri), type = 'l', lwd = 4, 
      col = cols[2],
      bty = 'n')
lines(2020:2040,rowMeans(prop_infected_3a2_tri), type = 'l', lwd = 4, 
      col = cols[3],
      bty = 'n')

lines(2020:2040,rowMeans(prop_infected_3b2_tri), type = 'l', lwd = 4, 
      col = cols[4],
      bty = 'n')
legend('topright', title = "", legend = c("Scenario 1_2", "Scenario 2_2", "Scenario 3a_2","Scenario 3b_2"),
       col = cols, lwd = c(4,4,4,4),
       bty = 'n', cex = 1.7)


#######################################################################################################################
#######################################################################################################################
#######################################################################################################################


plot(2020:2040, (rowMeans(t1)[1:23])/rowMeans(allNums1[1:23,]), type = 'l', ylim = c(0, max((rowMeans(t1)[1:23])/rowMeans(allNums1[1:23,]))),bty = 'n',
     lwd = 3, col = cols[1])
lines(2020:2040, (rowMeans(t3b1)[1:23])/rowMeans(allNums3b1[1:23,]), type = 'l', ylim = c(0, max(rowMeans(totalNumInfecteds_1)/3000)), bty = 'n',
      lwd = 3, col = cols[2])

lines(2020:2040, (rowMeans(t3b2)[1:23])/rowMeans(allNums3b2[1:23,]), type = 'l', ylim = c(0, max(rowMeans(totalNumInfecteds_1)/3000)), bty = 'n',
      lwd = 3, col = cols[3])
# 
# lines(2020:2040, (rowMeans(t3b)[1:23])/rowMeans(allNums3b[1:23,]), type = 'l', ylim = c(0, max(rowMeans(totalNumInfecteds_1)/3000)), bty = 'n',
#       lwd = 3, col = cols[4])

legend('topright', title = "", legend = c("Scenario 1_1", "Scenario 2_1", "Scenario 3a_1","Scenario 3b_1"),
       col = cols, lwd = c(4,4,4,4),
       bty = 'n', cex = 1.7)
