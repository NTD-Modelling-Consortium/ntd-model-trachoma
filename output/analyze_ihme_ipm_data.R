
get_prevs <- function(d, year, age_group, cat){
  return(d[1 + age_group*69 + (year-2018) + (cat-1)*23,])
}


get_numbers <- function(d, year, age_group){
  return(d[5521 + age_group + (year-2018)*80,])
}

get_data_over_years <- function(d, years,  age_group, cat, numAges){
  #age_groups = 0:79
  ks1 = 1 + age_group* length(years)* 4  + (years-min(years)) + (cat-1)*length(years) 
  numberStartPoint=  which(d$measure == "number")[1]
  ks2 = numberStartPoint  + age_group + (years-min(years))*numAges
  all_years_prev = d[ks1, ]
  all_years_number = d[ks2, ]
  
  return(list(all_years_prev, all_years_number))
}



get_number_infecteds <- function(ihme1, years){
  k = which(colnames(ihme1) == 'draw_0')
  age_groups = unique(ihme1$age_start)
  for(i in 1 : length(age_groups)){
    for(cat in 1:3){
      # output ihme and ipm data for given year, age group and prevalence intensity
      a = get_data_over_years(d = ihme1, years,  age_group = age_groups[i], cat = cat,  numAges = length(age_groups))
      
      # construct array to hold all the data
      if((i == 1 & cat == 1)){
        b = a[[1]][, k:ncol(ihme1)] * a[[2]][, k:ncol(ihme1)]
        num_infecteds = array(0, dim = c(dim(b)[1], dim(b)[2], length(age_groups), 3),
                              dimnames = list(rownames(b),
                                              colnames(b)))
        num_infecteds[, , 1, cat] = array(unlist(b))

      }else{
        num_infecteds[, , i, cat]= array(unlist(a[[1]][, k:ncol(ihme1)] * a[[2]][, k:ncol(ihme1)]))

      }
    }
    
    
  }
  return(num_infecteds)
}






get_data_over_years_trachoma<- function(d, years, measure = "TruePrevalence"){
  k = which(colnames(ihme1) == 'draw_0')
  for(i in 1:length(years)){
    year =  years[i]
    ks2 = (181+ (year-min(years))*240) :  (240 + (year-min(years))*240)
    
    if(measure == "TruePrevalence"){
      ks1 = (1 + (year-min(years))*240) :  (60 + (year-min(years))*240)
    }
    else if(measure == "ObservedTF"){
      ks1 <- (61 + (year-min(years))*240) :  (120 + (year-min(years))*240)
    }else if(measure == "heavyInfections"){
      ks1 <- (121 + (year-min(years))*240) :  (180 + (year-min(years))*240)
    }
    all_years_prev = d[ks1, k:ncol(d)]
    all_years_number = d[ks2, k:ncol(d)]
    if(i == 1){
      mean_infs = mean(colSums(all_years_prev * all_years_number))
    }else{
      mean_infs = c(mean_infs, mean(colSums(all_years_prev * all_years_number)))
    }
  }
  return(mean_infs)
}



get_total_number_infected <- function(num_infecteds, age_groups){
  total_number_infected_age_group = array(0, dim = c(dim(num_infecteds)[1], dim(num_infecteds)[2], length(age_groups)))
  for(i in 1 :length(age_groups)){
    total_number_infected_age_group[,,i] =  num_infecteds[, , i, 1] + num_infecteds[, , i, 2] +
      num_infecteds[, , i, 3]
  }
  
  total_number_infected = total_number_infected_age_group[,,1]
  #all_people = numbers_of_people[,,1]
  for(i in 2:length(age_groups)){
    total_number_infected = total_number_infected + total_number_infected_age_group[,,i]
    # all_people = all_people + numbers_of_people[,,i]
  }
  browser(expr = any(is.na(total_number_infected)) == TRUE)
  return(total_number_infected)
}


get_total_number_infected_age_range <- function(num_infecteds, age_range){
  total_number_infected_age_group = array(0, dim = c(dim(num_infecteds)[1], dim(num_infecteds)[2], length(age_range)))
  for(j in 1 :length(age_range)){
    i = age_range[j]
    total_number_infected_age_group[,,i] =  num_infecteds[, , i, 1] + num_infecteds[, , i, 2] +
      num_infecteds[, , i, 3]
  }
  
  total_number_infected = total_number_infected_age_group[,,1]
  #all_people = numbers_of_people[,,1]
  for(i in 2:length(age_range)){
    total_number_infected = total_number_infected + total_number_infected_age_group[,,i]
    # all_people = all_people + numbers_of_people[,,i]
  }
  return(total_number_infected)
}

get_surveyPass_data<- function(d){
  x = which(d$measure == 'surveyPass')
  survey  = d[x, ]
  
  years = unique(survey$year_id)
  
  c1 = which(colnames(nums) == "draw_0")
  
  return(survey[, c1:ncol(survey)])
}


get_number_people_in_age_range <- function(d, minAge, maxAge){
  x = which(d$measure == 'number')
  nums  = d[x, ]
  y = which(nums$age_start >= minAge & nums$age_end <= maxAge)
  nums = nums[y, ]
  years = unique(nums$year_id)
  allNums = matrix(0, length(years), 200)
  c1 = which(colnames(nums) == "draw_0")
  for(i in 1:length(years)){
    j = which(nums$year_id == years[i])
    allNums[i, ] = colSums(nums[j, c1:ncol(nums)])
  }
  return(allNums)
}




