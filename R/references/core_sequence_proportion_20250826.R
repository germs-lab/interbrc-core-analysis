library(tidyverse)

# ==== load file ====
raw <- read_csv(
  "/Users/jaejinlee/Files/Data/2023SABR_amplicon/analysis/seq-percentage_count_data.csv"
)

if (is.na(names(raw)[1]) || names(raw)[1] == "") {
  names(raw)[1] <- "Clustering"
}

# long format
long_df <- raw %>%
  pivot_longer(
    cols = -1,
    names_to = "Occurrence",
    values_to = "Proportion"
  ) %>%
  mutate(
    Occurrence = gsub("%", "", Occurrence),
    Occurrence = as.numeric(Occurrence),
    Clustering = as.character(.[[1]])
  )

# ==== axis, legend order ====
occ_levels <- c(50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100)
clust_levels <- c(85, 87, 90, 93, 95, 97, 100)

plot_df <- long_df %>%
  mutate(
    Occurrence = factor(
      Occurrence,
      levels = occ_levels,
      ordered = TRUE,
      labels = paste0(occ_levels, "%")
    ),
    Clustering = factor(
      Clustering,
      levels = paste0(clust_levels, "%"),
      ordered = TRUE
    )
  )

# color designation
custom_colors <- c(
  "85%" = "#FF68A1",
  "87%" = "#E68613",
  "90%" = "#0CB702",
  "93%" = "#00A9FF",
  "95%" = "#C77CFF",
  "97%" = "#999999",
  "100%" = "#ABA300"
)

p <- ggplot(
  plot_df,
  aes(x = Occurrence, y = Proportion, color = Clustering, group = Clustering)
) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2) +
  labs(
    title = "Core OTU Sequence Proportions across Clustering and Occurrence Thresholds",
    x = "Occurrence Threshold (%)",
    y = "Proportion of Sequences (%)",
    color = "Clustering\nThreshold"
  ) +
  scale_color_manual(values = custom_colors) +
  scale_y_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 10)
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 20, hjust = 0.5),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 14),
    legend.position = "right"
  )

print(p)
save(p, file = "/Users/jaejinlee/Files/Data/2023SABR_amplicon/analysis/p.rda")
