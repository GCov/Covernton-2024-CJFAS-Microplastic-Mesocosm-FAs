library(ggplot2)

theme1 <-
  theme_bw() +
  theme(
    panel.spacing = unit(1, "lines"), 
    text = element_text(size = 7),
    axis.text = element_text(size = 7),
    strip.background = element_blank(),
    strip.text = element_text(size = 8),
    legend.text = element_text(size = 8),
    panel.grid = element_blank()
  )

perch_diet <- read.csv("perch_diet.csv", header = TRUE)[1:20,]

perch_diet$dose <- as.numeric(perch_diet$dose)

ggplot(data = perch_diet,
       aes(x = dose,
           y = total.animals)) +
  geom_point() +
  labs(x = expression(paste("Dose (particles"~L^-1*")")),
       y = "Number of  animals in stomach") +
  scale_x_continuous(trans = "log1p",
                     breaks = c(0, 1, 10, 100, 1000, 10000, 30000)) +
  theme1
