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

# # add dummy column so code can be tested
# trachomaIUs$num_ius <-1

cumPop = 0
totalIUS = 0
count = 0

for (i in 1:nrow(trachomaIUs)){
  

  ipm_name = paste0(path_header,'endgame-ihme-ipm-outputs-trachoma-export-20220908a/ipm-',trachomaIUs$IU_ID2[i],"-trachoma-scenario_1-200_simulations.csv")
  if(file.exists(ipm_name)){
    ipm1 = read.csv(ipm_name)
    
    sPass = which(ipm1$measure == "surveyPass")
    k = which(colnames(ipm1) == 'draw_0')
    ipm1[sPass,k:ncol(ipm1)]
    if(i == 1){
      num_mda_finished_1 = rowMeans(ipm1[sPass,k:ncol(ipm1)])
    }else{
      num_mda_finished_1 = num_mda_finished_1 + 
        rowMeans(ipm1[sPass,k:ncol(ipm1)])
    }
    
    ipm_name = paste0(path_header,'endgame-ihme-ipm-outputs-trachoma-export-20220908a/ipm-',trachomaIUs$IU_ID2[i],"-trachoma-scenario_2a-200_simulations.csv")
    ipm2a = read.csv(ipm_name)
    
    sPass = which(ipm2a$measure == "surveyPass")
    
    k = which(colnames(ipm2a) == 'draw_0')
    ipm2a[sPass,k:ncol(ipm2a)]
    if(i == 1){
      num_mda_finished_2a = rowMeans(ipm2a[sPass,k:ncol(ipm2a)]) * trachomaIUs$num_ius[i]
    }else{
      num_mda_finished_2a = num_mda_finished_2a + 
        rowMeans(ipm2a[sPass,k:ncol(ipm2a)]) * trachomaIUs$num_ius[i]
    }
    
    ipm_name = paste0(path_header,'endgame-ihme-ipm-outputs-trachoma-export-20220908a/ipm-',trachomaIUs$IU_ID2[i],"-trachoma-scenario_2b-200_simulations.csv")
    
    ipm2b = read.csv(ipm_name)
    
    sPass = which(ipm2b$measure == "surveyPass")
    k = which(colnames(ipm2b) == 'draw_0')
    ipm2b[sPass,k:ncol(ipm2b)]
    if(i == 1){
      num_mda_finished_2b = rowMeans(ipm2b[sPass,k:ncol(ipm2b)]) * trachomaIUs$num_ius[i]
    }else{
      num_mda_finished_2b = num_mda_finished_2b + 
        rowMeans(ipm2b[sPass,k:ncol(ipm2b)]) * trachomaIUs$num_ius[i]
    }
    
    
    ipm_name = paste0(path_header,'endgame-ihme-ipm-outputs-trachoma-export-20220908a/ipm-',trachomaIUs$IU_ID2[i],"-trachoma-scenario_2c-200_simulations.csv")
    
    ipm2c = read.csv(ipm_name)
    
    sPass = which(ipm2c$measure == "surveyPass")
    k = which(colnames(ipm2c) == 'draw_0')
    ipm2c[sPass,k:ncol(ipm2c)]
    if(i == 1){
      num_mda_finished_2c = rowMeans(ipm2c[sPass,k:ncol(ipm2c)]) * trachomaIUs$num_ius[i]
    }else{
      num_mda_finished_2c = num_mda_finished_2c + 
        rowMeans(ipm2c[sPass,k:ncol(ipm2c)]) * trachomaIUs$num_ius[i]
    }
    count = count + 1
    totalIUS = totalIUS + trachomaIUs$num_ius[i]
  }
  print(paste("Done", i,"of", nrow(trachomaIUs)))
}




cols = c("#0098FF", "#1b9e77", "#d95f02","#7570b3")
prop_finished_1 = num_mda_finished_1/totalIUS
prop_finished_2a = num_mda_finished_2a/totalIUS
prop_finished_2b = num_mda_finished_2b/totalIUS
prop_finished_2c = num_mda_finished_2c/totalIUS
png("trachoma_mda_stopping.png", units = "in",
    width = 12, height = 8, res = 300)
plot(2020:2040, prop_finished_1 , type = 'l', lwd = 4, 
     xlab = 'year', ylab = 'proportion IUs stopped MDA',
     col = cols[1],
     bty = 'n', ylim = c(0, 1),
     cex = 1.7, cex.axis = 1.7, cex.lab = 1.7,
     cex.main = 1.7, main = "trachoma")

lines(2020:2040, prop_finished_2a, type = 'l', lwd = 4, 
      xlab = 'year', ylab = 'proportion IUs stop MDA',
      col = cols[2],
      bty = 'n')

lines(2020:2040, prop_finished_2b , type = 'l', lwd = 4, 
      xlab = 'year', ylab = 'proportion IUs stop MDA',
      col = cols[3],
      bty = 'n')

lines(2020:2040, prop_finished_2c , type = 'l', lwd = 4, 
      xlab = 'year', ylab = 'proportion IUs stop MDA',
      col = cols[4],
      bty = 'n')

legend('bottomright', title = "Scenario", legend = c("Scenario 1", "Scenario 2a", "Scenario 2b", "Scenario 2c"),
       col = cols, lwd = c(3,3,3, 3),
       bty = 'n',cex = 1.7)
dev.off()
