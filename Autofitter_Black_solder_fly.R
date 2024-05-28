#===============================
#..........Load data............
#===============================
rm(list = ls())
library("ggplot2")
library("readxl")
library("scales")
library("Rmisc")

data <- read.csv("000_Data_for_R.csv")

V_bluecap <- 78.38 / 10^6 #m^3
T <- 25 #deg C
s_to_day <- 60 * 60 * 24 #s/day
ppm_to_Pa <- 101325 / 1000000 #Pa/1000000 (ppm = atm/10^6)
Pa_to_mol <- 1 / ((273 + T) * 8.314) # mol/m^3
M_c <- 12.01 #g/mol


Respiration.C <- data$Respiration * s_to_day * ppm_to_Pa * Pa_to_mol * V_bluecap * M_c / 10 #g C/day/lv
Weight.C <- data$Weight_Indv * 0.49 #mg C lv
q <- Respiration.C*1000/Weight.C # mg C/(day * mg C-lv)

D <- cbind(data, Weight.C, Respiration.C, q)
D$Culture <- as.factor(D$Culture)

m.Weight <- summarySE(data = D, measurevar = "Weight_Indv", groupvars = c("Day", "Density", "Substrate"))


#=========================================================
#.........Finding growth parameters with fitter4..........     
#=========================================================

fitter4 = function(D){
  
  curve = seq(from = 0, to = 12, by = 0.05)
  names = unique(substr(D$Sample, 1, 9)) 
  Y = rep(0, length(names))
  m = rep(0, length(names))
  NGE_avg = rep(0, length(names))
  
  for (i in 1:length(names)){
    for (j in 1:nrow(D)){
      if (substr(D[[j, 1]], 1, 9) == names[i]){
        day_input = D[j, 5]
        weight_input = D[j, 9]
        q_input = D[j, 13]
        if (day_input == 1){
          day_list = day_input
          weight_list = weight_input
          q_list = q_input
        }else {
          day_list = rbind(day_list, day_input)
          weight_list = rbind(weight_list, weight_input)
          q_list = rbind(q_list, q_input)
        }
      
      }
    }
    
    input = data.frame(day_list, weight_list)
    
    
    fit = lm(weight_list ~ poly(day_list, 4, raw = TRUE))
    
    
    X_modeller<- function(curve, fit){
      x_list = rep(0, length(curve))
      for (i in 1:length(curve)) {
        x_list[i] <- coef(fit)[[1]] + coef(fit)[[2]] * curve[i] + coef(fit)[[3]]  * curve[i]^2 + coef(fit)[[4]] * curve[i]^3 + coef(fit)[[5]] * curve[i]^4
      }
      return(x_list)
    }
    
    X_mod_list <- X_modeller(curve, fit)

    
    u_modeller<- function(day_list, fit){
      u_list = rep(0, length(day_list))
      for (i in 1:length(day_list)) {
        u_list[i] <- (coef(fit)[[2]] + coef(fit)[[3]] * 2 * day_list[i] + coef(fit)[[4]] * 3 * day_list[i]^2 + coef(fit)[[5]] * 4 * day_list[i]^3)/(coef(fit)[[1]] + coef(fit)[[2]] * day_list[i] + coef(fit)[[3]]  * day_list[i]^2 + coef(fit)[[4]] * day_list[i]^3 + coef(fit)[[5]] * day_list[i]^4)
        }
      return(u_list)
    }
    
    
    u_mod_list <- u_modeller(day_list, fit)
    day_mod_list <- day_list
    
    png(paste(c("w_",names[i], "weigh.png"), collapse = "_"), width = 800, height = 800)
    plot(day_list, weight_list, pch = 19, main = paste(c(names[i], "Modelled vs measured weight"), collapse = "_"), xlim = c(1, 12), ylim = c(0, 80), ylab = "weight [mg/lv]", xlab = "day")
    lines(curve, X_mod_list, pch = 19, col = "red")
    dev.off()
    
    
    u_mod_list[c(1, length(u_mod_list))] <- NA
    q_list[c(1, length(q_list))] <- NA
    q_list[u_mod_list < 0] <- NA
    day_mod_list[c(1, length(u_mod_list))] <- NA
    day_mod_list[u_mod_list < 0] <- NA
    u_mod_list[u_mod_list < 0] <- NA
    u_mod_list <- na.omit(u_mod_list)
    u_mod_list <- u_mod_list
    q_list  <-  na.omit(q_list) 
    day_mod_list <- na.omit(day_mod_list)

    
    resp_fit = lm(q_list ~ u_mod_list)
    
    Y[i] = resp_fit$coefficients[[2]]
    m[i] = resp_fit$coefficients[[1]]
   
    png(paste(c("x_",names[i], "q_vs_u.png"), collapse = "_"), width = 800, height = 800)
    plot(u_mod_list, q_list, pch = 19, main = paste(c(names[i], "respiration_vs_growth"), collapse = "_"), xlim = c(0, 0.4), ylim = c(0, 0.1), abline(m[i], Y[i]), xlab = "Modelled specific growth rate [1/day]", ylab = "Specific respiration rate [1/day]")
    dev.off()

    
    NGE_calc <- function(q_list, u_mod_list){
      NGE_list = rep(0, length(q_list))
      for (i in 1:length(q_list)) {
        NGE_list[i] <- u_mod_list[i]/(u_mod_list[i]+q_list[i])
          }
      return(NGE_list)
    }
    
    NGE_list <- NGE_calc(q_list, u_mod_list)
    NGE_avg[i] <- mean(NGE_list)
    
    
    rm(day_list)
    rm(weight_list)
    rm(q_list)
    rm(NGE_list)
    }
  return(data.frame(names, Y, m, NGE_avg))
}

fitter4(D)

param <- fitter4(D)
param$m[param$m < 0] <- NA
param$Y[param$Y < 0] <- NA

term_data <- read.csv("Phenotypes_Termination.csv")
term_data <- term_data[-1]


Culture_data <- cbind(term_data[,-4], param[,-1])

#-----------------------------------------------
#.........Save data for next analysis...........
#-----------------------------------------------

write.csv(Culture_data, file = "culture_data.csv", row.names = FALSE)

