###Plotting

#
##
###
####
#####
#Clear workspace
rm(list=ls())

low_df <- read.csv("./output/time_to_extinction_df/m_0.005/TTE_r_mean=0.2 r_sd=0.35 K=5000 m=0.005 Num_sims=5000.csv")
low_df <- low_df[,-1]

med_df <- read.csv("./output/time_to_extinction_df/m_0.025/TTE_r_mean=0.2 r_sd=0.35 K=5000 m=0.025 Num_sims=5000.csv")
med_df <- med_df[,-1]

high_df <- read.csv("./output/time_to_extinction_df/m_0.05/TTE_r_mean=0.2 r_sd=0.35 K=5000 m=0.05 Num_sims=5000.csv")
high_df <- high_df[,-1]


mydiff <- function(data, diff){
  c(diff(data, lag = diff), rep(NA, diff))
}

x2 <- data.frame(v1=low_df,  v2 = mydiff(med_df, (length(med_df)-length(med_df))))


                 
                 
                 
x <- data.frame(v1=low_df,v2=med_df,v3=high_df)
library(ggplot2);library(reshape2)
data<- melt(x)
ggplot(data,aes(x=value, fill=variable)) + geom_density(alpha=0.25)
ggplot(data,aes(x=value, fill=variable)) + geom_histogram(alpha=0.25)
ggplot(data,aes(x=variable, y=value, fill=variable)) + geom_boxplot()
