## plot-slim-substutions.R
## by Rohan Maddamsetti
## Examine substitutions over time from SLiM simulation output.

library(tidyverse)

## Read the data
DR.sub.data <- read.csv("SLiM-substitutions-DR-50Kgen.csv") %>%
    mutate(Treatment="RQ-")

RQ.DR.sub.data <- read.csv("SLiM-substitutions-RQ-DR-50Kgen.csv") %>%
    mutate(Treatment="RQ+")

sub.data <- full_join(DR.sub.data, RQ.DR.sub.data)

 # --- Plot 1: Number of mutations per mutation_type ---
  p1 <- sub.data %>%
      count(Treatment, mutation_type) %>%
      ggplot(aes(x = mutation_type, y = n, fill = Treatment)) +
      geom_col(position = "dodge") +
      labs(x = "Mutation Type", y = "Number of Mutations",
           title = "Mutation Count by Type") +
      theme_minimal()

  # --- Plot 2: Number of fixation events (unique fixation_gen values) ---
  p2 <- sub.data %>%
      distinct(Treatment, fixation_gen) %>%
      count(Treatment) %>%
      ggplot(aes(x = Treatment, y = n, fill = Treatment)) +
      geom_col() +
      labs(x = "Treatment", y = "Number of Fixation Events",
           title = "Fixation Events (Unique Fixation Generations)") +
      theme_minimal()

  # --- Plot 3: Number of mutations per fixation event ---
p3 <- sub.data %>%                                                                                                                                                                 
    count(Treatment, fixation_gen) %>%                                                                                                                                             
    ggplot(aes(x = factor(fixation_gen), y = n, fill = Treatment)) +
    geom_col() +                                                                                                                                                                   
    facet_wrap(~Treatment, scales = "free_x") +
    labs(x = "Fixation Generation", y = "Mutations per Fixation",                                                                                                                  
         title = "Mutations per Fixation Event") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6))

p1

p2

p3


p3_alt <- sub.data %>%
      count(Treatment, fixation_gen) %>%
      ggplot(aes(x = Treatment, y = n, fill = Treatment)) +
      geom_boxplot() +
      labs(x = "Treatment", y = "Mutations per Fixation Event",
           title = "Distribution of Mutations per Fixation") +
      theme_minimal()

p3_alt

 # Table 1: Mutation count by type per treatment                                                                                                                                    
  t1 <- sub.data %>%                                                                                                                                                                 
      count(Treatment, mutation_type) %>%
      pivot_wider(names_from = Treatment, values_from = n, values_fill = 0)                                                                                                          
                  
  # Table 2: Number of unique fixation events per treatment                                                                                                                          
  t2 <- sub.data %>%
      distinct(Treatment, fixation_gen) %>%
      count(Treatment, name = "n_fixation_events")

  # Table 3: Mutations per fixation event per treatment
  t3 <- sub.data %>%
      count(Treatment, fixation_gen, name = "mutations_per_fixation")

t1

t2

t3

  - t1 pivots to wide format so each treatment is a column, making mutation type counts easy to compare side by side.
  - t2 gives one row per treatment with the count of unique fixation generations.
  - t3 gives one row per fixation event per treatment — you could further summarize it with:

  # Summary statistics of mutations per fixation
  t3_summary <- t3 %>%
      group_by(Treatment) %>%
      summarise(
          mean   = mean(mutations_per_fixation),
          median = median(mutations_per_fixation),
          min    = min(mutations_per_fixation),
          max    = max(mutations_per_fixation),
          sd     = sd(mutations_per_fixation)
      )
