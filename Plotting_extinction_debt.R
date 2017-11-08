########################################################################
########################################################################
#[5] Plotting results
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
