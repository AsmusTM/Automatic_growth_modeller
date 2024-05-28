#===============================
#..........Load data............
#===============================
rm(list = ls())
library("ggplot2")
library("readxl")
library("scales")
library("minpack.lm")
library("cowplot")
library("Rmisc")


data <- read_excel("00 Datatable for R.xlsx")


V_bluecap <- 78.38 / 10^6 #m^3
T <- 25 #deg C
s_to_day <- 60 * 60 * 24 #s/day
ppm_to_Pa <- 101325 / 1000000 #Pa/1000000 (ppm = atm/10^6)
Pa_to_mol <- 1 / ((273 + T) * 8.314) # mol/m^3
M_c <- 12.01 #g/mol


Respiration.C <- data$Respiration * s_to_day * ppm_to_Pa * Pa_to_mol * V_bluecap * M_c / 10 #g C/day/lv
Weight.C <- data$Dry.weight * 0.49 #mg C
q <- Respiration.C*1000/Weight.C # mg C/(day * mg lv)

D <- cbind(data, Weight.C, Respiration.C, q)
D$Culture <- as.factor(D$Culture)
m.Weight <- summarySE(data = D, measurevar = "Dry.weight", groupvars = c("Day", "Density", "Medium"))


#=========================================================
#.........Finding growth parameters with fitter3..........     
#=========================================================



fitter3 = function(D){
  f = function(Day, X_max, u_max, X_0 = 0.5388) {
    weight_list ~ X_max/(1+(X_max-X_0)/X_0*(exp(-u_max*(day_list-3))))
  }
  names = unique(substr(D$Sample, 1, 6)) 
  X_0_fit = rep(0, length(names))
  X_max_fit = rep(0, length(names))
  u_max_fit = rep(0, length(names))
  Y = rep(0, length(names))
  m = rep(0, length(names))
  NGE_max = rep(0, length(names))
  NGE_avg = rep(0, length(names))
  
  
  for (i in 1:length(names)){
    for (j in 1:nrow(D)){
      if (substr(D[[j, 1]], 1, 6) == names[i]){
        day_input = D[j, 5]
        weight_input = D[j, 10]
        q_input = D[j, 13]
        if (day_input == 3){
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
    
    fit = nlsLM(f(X_0 = weight_list[1]), data = input, 
                start=list(X_max=8, u_max= 2))
    
    X_0_fit[i] = weight_list[1]
    X_max_fit[i] = coef(fit)[[1]]
    u_max_fit[i] = coef(fit)[[2]]
    
    X_mod_list = rep(0, length(day_list))
    u_mod_list = rep(0, length(day_list))
    NGE_list = rep(0, length(day_list))
    
    for (k in 1:length(day_list)){
      X_mod_list[k] = X_max_fit[i]/(1+((X_max_fit[i] - X_0_fit[i])/X_0_fit[i])*exp(-u_max_fit[i]*(day_list[k]-3)))
      u_mod_list[k] = u_max_fit[i] * (1 - X_mod_list[k]/X_max_fit[i])
      NGE_list[k] = u_mod_list[k]/(u_mod_list[k]+q_list[k])
    }
    
    png(paste(c("w_",names[i], "weigh.png"), collapse = "_"), width = 800, height = 800)
    plot(day_list, weight_list, pch = 19, main = paste(c(names[i], "Modelled vs measured weight"), collapse = "_"), xlim = c(3, 8), ylim = c(0, 10), ylab = "weight [mg/lv]", xlab = "day")
    lines(day_list, X_mod_list, pch = 19, col = "red")
    dev.off()
    
    
    resp_fit = lm(q_list ~ u_mod_list)
    Y[i] = resp_fit$coefficients[[2]]
    m[i] = resp_fit$coefficients[[1]]
    NGE_avg[i] = mean(NGE_list)
    
    
    png(paste(c("x_",names[i], "q_vs_u.png"), collapse = "_"), width = 800, height = 800)
    plot(u_mod_list, q_list, pch = 19, main = paste(c(names[i], "respiration_vs_growth"), collapse = "_"), xlim = c(0, 2.5), ylim = c(0, 0.6), abline(m[i], Y[i]), xlab = "Modelled specific growth rate [1/day]", ylab = "Specific respiration rate [1/day]")
    dev.off()
    
    
    rm(day_list)
    rm(weight_list)
    rm(q_list)
    rm(NGE_list)
    
  }
  
  return(data.frame(names, X_0_fit, X_max_fit, u_max_fit, Y, m, NGE_avg))
}


param <- fitter3(D)

term_data <- read.csv("Phenotypes_Termination.csv")
term_data <- term_data[-1]


Culture_data <- cbind(term_data[,-4], param[,-1])

#-----------------------------------------------
#.........Save data for next analysis...........
#-----------------------------------------------


write.csv(Culture_data, file = "culture_data.csv", row.names = FALSE)

