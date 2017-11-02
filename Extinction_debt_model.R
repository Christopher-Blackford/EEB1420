##############################################################################################################################
##############################################################################################################################
###Extinction_debt_model.R
#Code by: Christopher Blackford and Fielding Montgomery
#Copyright © 2017

###READ.ME:
#This file models metapopulation dynamics of extinction debt

###################TABLE OF CONTENTS
###Initializing
###[1] Setting up the model
###[2] Running model


#Loading in libraries
require(deSolve)
require(ggplot2)
#
##
###
####
#####
#Clear workspace
rm(list=ls())
full.run.time <- proc.time() # time for one run-through


########################################################################
########################################################################
#Initializing showing how you want r to vary over time, and what your intial population size will be

#How do you want the r to vary over time? (r > 0 = population increasing, r < 0 = population decreasing)
r_loop <- rnorm(20, mean = 0.02, sd = 0.1)
r_loop[1] <- 0.02

#Setting up dataframe to capture intial and final population sizes to feed into runs when you loop
population_loop <- matrix(0, nrow = length(r_loop), ncol = 4)

#These are the initial population sizes for the different patches
population_loop[1,1] <- 20 
population_loop[1,2] <- 5
population_loop[1,3] <- 0
population_loop[1,4] <- 0

#This a a vector that I'm using to give names to the different loop outputs
Model_output_names <- NULL
for (i in 1:length(r_loop)){
Model_output_names[i] <- append(paste0("loop_", i), Model_output_names)
}

#How long should each loop run for?
time_of_loop <- 20

########################################################################
########################################################################
#[1] Setting up the model

for (i in 1:length(r_loop)){

##### Model description  #####
SISmodel=function(t,y,parameters){ 
  ## Variables
  S1_A=y[1]; S1_B=y[2]; S1_C=y[3]; S1_D=y[4]
  
  ## Parameters
  r = parameters[1];
  K = parameters[2];
  m_AB = parameters[3];
  m_AD = parameters[4];
  m_BA = parameters[5];
  m_BC = parameters[6];
  m_CB = parameters[7];  
  m_CD = parameters[8];  
  m_DA = parameters[9];  
  m_DC = parameters[10]; 
  
  ## Ordinary differential equations
  dS1_Adt <- r*S1_A*(1 - S1_A/K) + S1_B*m_BA + S1_D*m_DA - S1_A*(m_AB + m_AD)
  
  dS1_Bdt <- r*S1_B*(1 - S1_B/K) + S1_A*m_AB + S1_C*m_CB - S1_B*(m_BA + m_BC)
  
  dS1_Cdt <- r*S1_C*(1 - S1_C/K) + S1_B*m_BC + S1_D*m_DC - S1_C*(m_CB + m_DC)
  
  dS1_Ddt <- r*S1_D*(1 - S1_D/K) + S1_A*m_AD + S1_C*m_CD - S1_D*(m_DA + m_DC)
  
  
  return(list(c(dS1_Adt,
                dS1_Bdt,
                dS1_Cdt,
                dS1_Ddt))); 
}  


########################################################################
########################################################################
#[2] Parameter values
r <- r_loop[i] #r value will change with each loop
K <- 5000
m_AB <- 0.01 #Setting all migration rates equal
m_AD = m_AB
m_BA = m_AB
m_BC = m_AB
m_CB = m_AB
m_CD = m_AB
m_DA = m_AB
m_DC = m_AB



########################################################################
########################################################################
#[3] Running the model

## Initial state
variables0=c(S1_A0=population_loop[i,1],
             S1_B0=population_loop[i,2],
             S1_C0=population_loop[i,3],
             S1_D0=population_loop[i,4]) #Initial population size will be updated for each loop

## Times at which estimates of the variables are returned
timevec=seq(0,time_of_loop,1) #Defines how long to run the model before looping

parameters=c(r,
             K,
             m_AB,
             m_AD,
             m_BA,
             m_BC,
             m_CB,
             m_CD,
             m_DA,
             m_DC)

#Model run
output=lsoda(y = variables0,    # intial values  
             times = timevec,   # time vector
             func = SISmodel,   # model
             parms = parameters # constant parameters
)


########################################################################
########################################################################
#[4] Getting output into dataframe
colnames(output)=c("time","S1_A", "S1_B", "S1_C", "S1_D")
Metapop_model=as.data.frame(output)
Metapop_model$loop_number <- i

assign(Model_output_names[i], Metapop_model)


if (i < length(r_loop)){
population_loop[i+1,1] <- Metapop_model$S1_A[nrow(Metapop_model)]
population_loop[i+1,2] <- Metapop_model$S1_B[nrow(Metapop_model)]
population_loop[i+1,3] <- Metapop_model$S1_C[nrow(Metapop_model)]
population_loop[i+1,4] <- Metapop_model$S1_D[nrow(Metapop_model)]
}

else print("Done")

}



Model_output <- lapply(ls(pattern=paste0("loop_")), function(x) get(x))
Model_output <- do.call(rbind, Model_output)
Model_output <- Model_output[with(Model_output, order(loop_number, time)),]
Model_output$time_loop <- Model_output$time
Model_output$time <- 1:nrow(Model_output)


Linear_plot <- ggplot(Model_output, aes(time)) +
  geom_line(aes(y = S1_A, colour = "Species 1, patch A")) + 
  geom_line(aes(y = S1_B, colour = "Species 1, patch B")) +
  geom_line(aes(y = S1_C, colour = "Species 1, patch C")) +
  geom_line(aes(y = S1_D, colour = "Species 1, patch D")) +
  labs(title = "", x = "Time", y = "Population size")

Linear_plot <- Linear_plot + theme(
  plot.title = element_text(size = 16), 
  axis.text = element_text(size = 16),
  axis.title = element_text(size = 16),
  axis.line = element_line("black"),
  legend.title = element_blank(),
  panel.background = element_blank()
)

Linear_plot
