##############################################################################################################################
##############################################################################################################################
###Extinction_debt_model.R
#Code by: Christopher Blackford and Fielding Montgomery
#Copyright © 2017

###READ.ME:
#This file models metapopulation dynamics of extinction debt

###################TABLE OF CONTENTS
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
#[1] Setting up the model

##### Model description  #####
SISmodel=function(t,y,parameters) 
{ 
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
r <- 1.05
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
variables0=c(S1_A0=20, S1_B0=5, S1_C0=0, S1_D0=0)

## Times at which estimates of the variables are returned
timevec=seq(0,20,1)

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
#[4] Plotting
colnames(output)=c("time","S1_A", "S1_B", "S1_C", "S1_D")
Metapop_model=as.data.frame(output)


Linear_plot <- ggplot(Metapop_model, aes(time)) +
  geom_line(aes(y = S1_A, colour = "Species 1, patch A")) + 
  geom_line(aes(y = S1_B, colour = "Species 1, patch B")) +
  geom_line(aes(y = S1_C, colour = "Species 1, patch C")) +
  geom_line(aes(y = S1_D, colour = "Species 1, patch D")) +
  labs(title = "", x = "Time", y = "Population size")

Linear_plot + theme(
  plot.title = element_text(size = 16), 
  axis.text = element_text(size = 16),
  axis.title = element_text(size = 16),
  axis.line = element_line("black"),
  legend.title = element_blank(),
  panel.background = element_blank()
)



