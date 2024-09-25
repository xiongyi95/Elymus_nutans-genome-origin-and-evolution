# Load necessary libraries
library(dplyr)
library(ggtern)
library(fields) # for the rdist function

# Read in the data
# Assume triads_data.csv contains columns: Gene_ID, St_tpm, Y_tpm, H_tpm, Stress (UV-B, Cold, Drought, Control)
triads_data <- read.csv("triads_data.csv")

# Function to calculate Total TPM and normalize gene expression values
normalize_expression <- function(df) {
  df <- df %>%
    mutate(Total_tpm = St_tpm + Y_tpm + H_tpm) %>%
    filter(Total_tpm > 0.5) %>%
    mutate(
      St_norm = St_tpm / Total_tpm,
      Y_norm = Y_tpm / Total_tpm,
      H_norm = H_tpm / Total_tpm
    )
  return(df)
}

# Normalize the expression data for each stress condition
normalized_data <- normalize_expression(triads_data)

# Define ideal categories for expression patterns
ideal_categories <- data.frame(
  Category = c("St-dominant", "St-suppressed", "Y-dominant", "Y-suppressed", 
               "H-dominant", "H-suppressed", "Balanced"),
  St_ideal = c(1, 0, 0, 0, 0.5, 0.5, 1/3),
  Y_ideal = c(0, 0.5, 1, 0, 0.5, 0, 1/3),
  H_ideal = c(0, 0.5, 0, 0.5, 1, 0, 1/3)
)

# Function to assign categories based on Euclidean distance
assign_category <- function(normalized_df, ideal_df) {
  normalized_df$Category <- apply(normalized_df, 1, function(row) {
    dist <- sapply(1:nrow(ideal_df), function(i) {
      rdist(c(row["St_norm"], row["Y_norm"], row["H_norm"]), 
            c(ideal_df$St_ideal[i], ideal_df$Y_ideal[i], ideal_df$H_ideal[i]))
    })
    ideal_df$Category[which.min(dist)]
  })
  return(normalized_df)
}

# Assign categories to each triad
categorized_data <- assign_category(normalized_data, ideal_categories)

# Visualize relative expressions using ggtern
ggtern(data = categorized_data, aes(x = St_norm, y = Y_norm, z = H_norm, color = Category)) +
  geom_point() +
  theme_bw() +
  labs(title = "Triad Expression Patterns", x = "St", y = "Y", z = "H") +
  theme_showarrows() +
  scale_color_manual(values = c("red", "blue", "green", "orange", "purple", "brown", "black"))

# Save the plot
ggsave("triad_expression_patterns.png")

# Save the categorized data for further analysis
write.csv(categorized_data, "categorized_triads.csv", row.names = FALSE)

print("Analysis completed and results saved.")
