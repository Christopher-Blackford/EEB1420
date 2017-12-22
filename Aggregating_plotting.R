###Aggregating plotting

file1 <- read.csv("./output/time_to_extinction_df/TTE_r_mean=0.2 r_sd=0.35 K=5000 m=0.025 Num_sims=25000.csv")
file2 <- read.csv("./output/time_to_extinction_df/TTE_r_mean=0.2 r_sd=0.35 K=5000 m=0.025 Num_sims=75000.csv")
file3 <- read.csv("./output/time_to_extinction_df/TTE_r_mean=0.2 r_sd=0.35 K=5000 m=0.025 Num_sims=150000.csv")

merged <- rbind(file1, file2)
merged <- rbind(merged, file3)

#merged <- as.data.frame(merged[,2])


###Histogram
plot_title <- "Merged plot, r_mean=0.2 r_sd=0.35 K=5000 m=0.025 Num_sims=250000"

Time_to_extinction_plot <- ggplot(merged, aes(time_quasi_extinct)) +
  geom_histogram(colour = "Black", bins = nrow(merged))+
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
