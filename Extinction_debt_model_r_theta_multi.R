##############################################################################################################################
##############################################################################################################################
###Extinction_debt_model.R
#Code by: Christopher Blackford and Fielding Montgomery
#Copyright © 2017

###READ.ME:
#This file models metapopulation dynamics of extinction debt
#It looks at how long it takes a metapopulation to go to extinction after disturbance depending on:
# Species growth rates (r)
# Metapopulation connectivity
#
#

###################TABLE OF CONTENTS
###[1] Defining model parameters
###[2] Looping across simulations
###[3] Looping model parameters
###[4] Building model structure
###[5] Running the model
###[6] Getting model output into dataframe/plotting single simulation results
###[7] Getting time to extinction for a run
###[8] Plotting time to extinction across model simulations

###################Loading in libraries
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


pnorm(0,mean=0.2, sd=0.35) #calculate percent of bad years (where r< 0) given mean and sd

########################################################################
########################################################################
#[1] Defining Model Parameters
number_of_simulations <- 5000  #How many simulations to do
years_each_run <- 500 #How long should each simulation run for

r_mean <- c(0.2, 0.15, 0.10) #r grows at "x" percent per time step

r_sd <- 0.35

#Defining K
K_all <- c(5000, 3750, 2500)

#Defining m (migration)
m_all <- 0.05

#How long should each loop run for?
time_of_loop <- 1

#Quasi-extinction threshold
extinction_threshold <- 500 #10% of carrying capacity


for (Multiple_r_scenarios in 1:length(r_mean)){
  for (Multiple_K_scenarios in 1:length(K_all)){
    

########################################################################
########################################################################
#[2] Looping across simulations

time_quasi_extinct <- NULL #Time to quasi-extinction vector
for (simulation_dummy in 1:number_of_simulations){

########################################################################
########################################################################
#[3] Looping model parameters
  
#For different r in each patch
r_list <- c("rA_loop", "rB_loop", "rC_loop", "rD_loop")
for (r_listing in 1:length(r_list)){
r_loop <- rnorm(years_each_run, mean = r_mean[Multiple_r_scenarios], sd = r_sd) 
assign(r_list[r_listing], r_loop)}


#Setting up dataframe to capture intial and final population sizes to feed into runs when you loop
population_loop <- matrix(0, nrow = length(r_loop), ncol = 4)
#These are the initial population sizes for the different patches
population_loop[1,1] <- K_all[Multiple_K_scenarios]/2 
population_loop[1,2] <- K_all[Multiple_K_scenarios]/2
population_loop[1,3] <- K_all[Multiple_K_scenarios]/2
population_loop[1,4] <- K_all[Multiple_K_scenarios]/2
  
#This a a vector that I'm using to give names to the different loop outputs
Model_output_names <- NULL
for (i in 1:length(r_loop)){Model_output_names[i] <- append(paste0("loop_", i), Model_output_names)}

rm(r_list, r_listing) #removing usless variables
#######################################################################
########################################################################
#[4] Building model structure
for (i in 1:length(r_loop)){

##### Model description  #####
SISmodel=function(t,y,parameters){ 
  ## Variables
  S1_A=y[1]; S1_B=y[2]; S1_C=y[3]; S1_D=y[4]
  
  ## Parameters
  rA = parameters[1];
  rB = parameters[2];
  rC = parameters[3];
  rD = parameters[4];
  K = parameters[5];
  m_AB = parameters[6];
  m_AD = parameters[7];
  m_BA = parameters[8];
  m_BC = parameters[9];
  m_CB = parameters[10];  
  m_CD = parameters[11];  
  m_DA = parameters[12];  
  m_DC = parameters[13]; 
  theta_A = parameters[14];
  theta_B = parameters[15];
  theta_C = parameters[16];
  theta_D = parameters[17];
  
  ## Ordinary differential equations
  dS1_Adt <- rA*S1_A*(1 - S1_A/K)^theta_A + S1_B*m_BA + S1_D*m_DA - S1_A*(m_AB + m_AD)
  
  dS1_Bdt <- rB*S1_B*(1 - S1_B/K)^theta_B + S1_A*m_AB + S1_C*m_CB - S1_B*(m_BA + m_BC)
  
  dS1_Cdt <- rC*S1_C*(1 - S1_C/K)^theta_C + S1_B*m_BC + S1_D*m_DC - S1_C*(m_CB + m_DC)
  
  dS1_Ddt <- rD*S1_D*(1 - S1_D/K)^theta_D + S1_A*m_AD + S1_C*m_CD - S1_D*(m_DA + m_DC)
  
  
  return(list(c(dS1_Adt,
                dS1_Bdt,
                dS1_Cdt,
                dS1_Ddt))); 
}  

###[4b] Parameter values
rA <- rA_loop[i] #r value will change with each loop
rB <- rB_loop[i]
rC <- rC_loop[i]
rD <- rD_loop[i]
K <- K_all[Multiple_K_scenarios]
m_AB <- m_all #Setting all migration rates equal
m_AD = m_AB
m_BA = m_AB
m_BC = m_AB
m_CB = m_AB
m_CD = m_AB
m_DA = m_AB
m_DC = m_AB
if (rA_loop[i] >= 0){theta_A = 1}else (theta_A = 0)
if (rB_loop[i] >= 0){theta_B = 1}else (theta_B = 0)
if (rC_loop[i] >= 0){theta_C = 1}else (theta_C = 0)
if (rD_loop[i] >= 0){theta_D = 1}else (theta_D = 0)

########################################################################
########################################################################
#[5] Running the model

#Initial state
#Need to do this so that it takes the initial starting size the first time, and then the output from the last fun the next time
variables0=c(S1_A0=population_loop[i,1],
                     S1_B0=population_loop[i,2],
                     S1_C0=population_loop[i,3],
                     S1_D0=population_loop[i,4]) #Initial population size will be updated for each loop

#Times at which estimates of the variables are returned
timevec=seq(0,time_of_loop,1) #Defines how long to run the model before looping

parameters=c(rA,
             rB,
             rC,
             rD,
             K,
             m_AB,
             m_AD,
             m_BA,
             m_BC,
             m_CB,
             m_CD,
             m_DA,
             m_DC,
             theta_A,
             theta_B,
             theta_C,
             theta_D)

#Model run
output=lsoda(y = variables0,    # intial values  
             times = timevec,   # time vector
             func = SISmodel,   # model
             parms = parameters # constant parameters
)


########################################################################
########################################################################
#[6] Getting model output into dataframe/plotting single simulation results
colnames(output)=c("time","S1_A", "S1_B", "S1_C", "S1_D")
Metapop_model=as.data.frame(output)
Metapop_model$loop_number <- i
#Need to remove the way the output is structure since the first row is a repetition of previous abundances
if (i > 1){Metapop_model = Metapop_model[-1,]}
assign(Model_output_names[i], Metapop_model) #adding on time period to previous time step


#Start next simulation and end of previous simulation
if (i < length(r_loop)){
population_loop[i+1,1] <- Metapop_model$S1_A[nrow(Metapop_model)]
population_loop[i+1,2] <- Metapop_model$S1_B[nrow(Metapop_model)]
population_loop[i+1,3] <- Metapop_model$S1_C[nrow(Metapop_model)]
population_loop[i+1,4] <- Metapop_model$S1_D[nrow(Metapop_model)]
} #else do nothing

} #close r_loop loop

#Synthesizing across multiple r values in a single simulation
Model_output <- lapply(ls(pattern=paste0("loop_")), function(x) get(x))
Model_output <- do.call(rbind, Model_output)
Model_output <- Model_output[with(Model_output, order(loop_number, time)),]
Model_output$time_loop <- Model_output$time
Model_output$time <- 1:nrow(Model_output)

rm(list=ls(pattern="loop_")) #Removing needless loop files
########################################################################
########################################################################
#[7] Getting time to extinction for a run

for (i in 1:nrow(Model_output)){
  if (Model_output[i,2] < extinction_threshold){
    time_quasi_extinct <- append(time_quasi_extinct, Model_output[i,1])
    break
  } #don't need an else statement
}

###Progress bar
if (simulation_dummy == ceiling(number_of_simulations*0.25)){
  print("25%")}
else if (simulation_dummy == ceiling(number_of_simulations*0.50)){
  print("50%")}
else if (simulation_dummy == ceiling(number_of_simulations*0.75)){
  print("75%")}
else if (simulation_dummy == ceiling(number_of_simulations)){
  print("Done")}

} #closing simulation loop


########################################################################
########################################################################
#[8] Plotting time to extinction across model simulations
base_title <- paste0("r_mean=", r_mean[Multiple_r_scenarios], " r_sd=", r_sd, " K=", K_all[Multiple_K_scenarios], " m=", m_all, " Num_sims=", number_of_simulations)
base_title

Time_to_extinction_df <- as.data.frame(time_quasi_extinct)
write.csv(Time_to_extinction_df, file = paste0("./output/time_to_extinction_df/TTE_", base_title,".csv"))

###Histogram
plot_title <- paste0("Hist ", base_title)

Time_to_extinction_plot <- ggplot(Time_to_extinction_df, aes(time_quasi_extinct)) +
  geom_histogram(colour = "Black", bins = nrow(Model_output))+
  labs(title = plot_title, x = "Time to extinction", y = "Frequency")

Time_to_extinction_plot <- Time_to_extinction_plot + theme(
  plot.title = element_text(size = 16), 
  axis.text = element_text(size = 14),
  axis.title = element_text(size = 14),
  axis.line = element_line("black"),
  legend.title = element_blank(),
  panel.background = element_blank()
)
Time_to_extinction_plot

ggsave(paste0("./output/figures/r_theta/", plot_title, ".png"), width = 10, height = 6)

###Frequency histogram
plot_title <- paste0("Freq ", base_title)

Time_to_extinction_freq <- ggplot(Time_to_extinction_df, aes(time_quasi_extinct)) +
  geom_freqpoly(colour = "Black", bins = nrow(Model_output)/10)+
  labs(title = plot_title, x = "Time to extinction", y = "Frequency")

Time_to_extinction_freq <- Time_to_extinction_freq + theme(
  plot.title = element_text(size = 16), 
  axis.text = element_text(size = 14),
  axis.title = element_text(size = 14),
  axis.line = element_line("black"),
  legend.title = element_blank(),
  panel.background = element_blank()
)
Time_to_extinction_freq

ggsave(paste0("./output/figures/r_theta/", plot_title, ".png"), width = 10, height = 6)


  } #closing multiple K scenarios
} #closing multiple r scenarios


proc.time() - full.run.time
#####
####
###
##
#END

