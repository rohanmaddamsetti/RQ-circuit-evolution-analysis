## plot-slim-fitness.R
## by Rohan Maddamsetti
## Plot mean fitness over time from SLiM simulation output.

library(tidyverse)

## Read the data
DR.data <- read.csv("SLiM-output-DR-50Kgen.csv") %>%
    select(generation, mean_fitness) %>%
    mutate(Treatment="RQ-")

RQ.DR.data <- read.csv("SLiM-output-RQ-DR-50Kgen.csv") %>%
    select(generation, mean_intrinsic_fitness) %>%
    ## CONFUSING: we want mean intrinsic fitness, will
    ## rename to mean fitness for comparisons.
    rename(mean_fitness = mean_intrinsic_fitness) %>%
    mutate(Treatment="RQ+")

fitness.data <- full_join(DR.data, RQ.DR.data)

fitness.plot <- fitness.data %>%
    ggplot(aes(x = generation, y = mean_fitness, color=Treatment)) +
    geom_line(linewidth=1) +
    labs(
        x = "Generation",
        y = "Mean Fitness",
        title = "Red Queen (RQ) gene circuit speeds adaptation"
    ) +
    theme_classic(base_size = 10) +
    theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position="bottom"
    )

## Save the plot
ggsave("RQ-fitness-comparison.pdf", fitness.plot, width=5, height = 3)
