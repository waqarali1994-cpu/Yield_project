
library(dplyr)

# Define which columns you want to keep
columns_to_keep <- c("genotype", "year", "location", "nitrogenTreatment",
                     "kernelRowNumber", "hundredKernelMass", "kernelMassPerEar",
                     "earWeight", "earWidth", "earLength", "kernelsPerRow")

# Apply filtering and selection
phenotypic_data_filtered <- phenotypic_data %>%
  filter(nitrogenTreatment == "Medium") %>%
  select(all_of(columns_to_keep))


write.csv(phenotypic_data_filtered,"filtered_phenotype_data.csv")



# Load dataset
snp_trait_raw <- fread("newyieldgwas_z.csv")  # Replace with your actual filename

# Combine CHROM and POS into a new SNP column
snp_trait_clean <- snp_trait_raw %>%
  mutate(SNP = paste0(CHROM, "_", POS)) %>%
  select(SNP, Trait)  # Reorder: SNP first, then Trait

# View result
print(snp_trait_clean)

# Save to CSV
fwrite(snp_trait_clean, "combinedsnp.csv", row.names = FALSE)



setwd("~/Documents/PhD work/Yield_project/SNP effect/SNP_effect_v2/")

library(data.table)
library(dplyr)

# ---- Step 1: Load your raw SNP–trait mapping file ----
snp_trait_raw <- fread("combinedsnp.csv")  # Replace with your actual filename

# ---- Step 2: Combine CHROM and POS to create SNP ID ----
snp_trait_clean <- snp_trait_raw %>%
  mutate(SNP = paste0(CHROM, "_", POS)) %>%
  select(SNP, Trait)

# ---- Step 3: Define a mapping from readable trait names to exact column names ----
trait_name_map <- c(
  "Kernel Weight" = "hundredKernelMass",
  "Plot Weight" = "kernelMassPerEar",
  "Kernel Row Number" = "kernelRowNumber",
  "Kernels Per Row" = "kernelsPerRow",
  "Ear Width" = "earWidth",
  "Ear Length" = "earLength",
  "Ear Weight" = "earWeight"
)

# ---- Step 4: Standardize trait names ----
snp_trait_clean <- snp_trait_clean %>%
  mutate(Trait = trait_name_map[Trait]) %>%
  filter(!is.na(Trait))  # Remove any unmatched trait names

# ---- Step 5: Save the cleaned SNP–Trait map ----
fwrite(snp_trait_clean, "combinedsnp.csv", row.names = FALSE)



# ---- Load Required Libraries ----
library(data.table)
library(dplyr)
library(tidyr)

# ---- Step 1: Load SNP-Trait Mapping ----
snp_trait_map <- fread("combinedsnp.csv") %>%
  filter(!is.na(Trait) & Trait != "")  # Remove empty or NA trait names

# ---- Step 2: Load Genotype and Phenotype Data ----
allele_data <- fread("transposed.csv")           # SNPs as columns, 'genotype' in first column
phenotypic_data <- fread("filtered_phenotype_data.csv")

# ---- Step 3: Merge and Create Environment ----
merged_data <- inner_join(phenotypic_data, allele_data, by = "genotype") %>%
  filter(nitrogenTreatment == "Medium") %>%
  mutate(Environment = paste(location, year))

# ---- Step 4: No global homozygous filtering ----
filtered_data <- merged_data  # Keep full dataset, filter per SNP inside loop

write.csv(filtered_data, "filtered_data_withsnp.csv")
# ---- Step 5: Function to Compute Per-SNP Effect Sizes ----
compute_all_effects <- function(data, snp_trait_map) {
  results <- list()
  
  for (i in seq_len(nrow(snp_trait_map))) {
    snp <- snp_trait_map$SNP[i]
    trait <- snp_trait_map$Trait[i]
    
    # Skip invalid mappings
    if (is.na(snp) || is.na(trait) || trait == "" || !(trait %in% colnames(data))) next
    
    message("Processing: ", snp, " → ", trait)
    
    # 🔸 Filter to homozygotes for this SNP only
    snp_data <- data %>%
      filter(.data[[snp]] %in% c("0|0", "1|1"))
    
    # Compute mean trait per genotype per environment
    effect_df <- snp_data %>%
      group_by(Env = Environment, Genotype = .data[[snp]]) %>%
      summarise(mean_trait = mean(.data[[trait]], na.rm = TRUE), .groups = "drop") %>%
      pivot_wider(names_from = Genotype, values_from = mean_trait)
    
    # Assign values safely
    ref_val <- if ("0|0" %in% colnames(effect_df)) effect_df$`0|0` else NA_real_
    alt_val <- if ("1|1" %in% colnames(effect_df)) effect_df$`1|1` else NA_real_
    
    # Compute percent effect size
    effect_df <- effect_df %>%
      mutate(
        SNP = snp,
        Trait = trait,
        REF = ref_val,
        ALT = alt_val,
        Effect_Size = ((ALT / REF) - 1) * 100
      ) %>%
      select(SNP, Trait, Env, REF, ALT, Effect_Size)
    
    results[[paste(snp, trait, sep = "_")]] <- effect_df
  }
  
  bind_rows(results)
}

# ---- Step 6: Run and Save Final Output ----
all_effects <- compute_all_effects(filtered_data, snp_trait_map)
fwrite(all_effects, "Final_Effect_Summary_SpecificSNPs.csv")













# ---- Load Required Libraries ----
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)

# ---- Step 1: Preprocess Phenotype Data ----
phenotypic_data <- fread("your_original_phenotype_file.csv")  # Update if needed

columns_to_keep <- c("genotype", "year", "location", "nitrogenTreatment",
                     "kernelRowNumber", "hundredKernelMass", "kernelMassPerEar",
                     "earWeight", "earWidth", "earLength", "kernelsPerRow")

phenotypic_data_filtered <- phenotypic_data %>%
  filter(nitrogenTreatment == "Medium") %>%
  select(all_of(columns_to_keep))

fwrite(phenotypic_data_filtered, "filtered_phenotype_data.csv")


# ---- Step 2: Prepare SNP-Trait Mapping ----
snp_trait_raw <- fread("newyieldgwas_z.csv")

trait_name_map <- c(
  "Kernel Weight" = "hundredKernelMass",
  "Plot Weight" = "kernelMassPerEar",
  "Kernel Row Number" = "kernelRowNumber",
  "Kernels Per Row" = "kernelsPerRow",
  "Ear Width" = "earWidth",
  "Ear Length" = "earLength",
  "Ear Weight" = "earWeight"
)

snp_trait_clean <- snp_trait_raw %>%
  mutate(SNP = paste0(CHROM, "_", POS)) %>%
  mutate(Trait = trait_name_map[Trait]) %>%
  filter(!is.na(Trait)) %>%
  select(SNP, Trait)

fwrite(snp_trait_clean, "combinedsnp.csv")




# ---- Step 3: Load and Merge Datasets ----
allele_data <- fread("transposed.csv")
phenotypic_data <- fread("filtered_phenotype_data.csv")
snp_trait_map <- fread("combinedsnp.csv")

merged_data <- inner_join(phenotypic_data, allele_data, by = "genotype") %>%
  filter(nitrogenTreatment == "Medium") %>%
  mutate(Environment = paste(location, year))


# ---- Step 4: Compute Effect Sizes (Per-SNP, Homozygous Only) ----
compute_all_effects <- function(data, snp_trait_map) {
  results <- list()
  
  for (i in seq_len(nrow(snp_trait_map))) {
    snp <- snp_trait_map$SNP[i]
    trait <- snp_trait_map$Trait[i]
    
    if (is.na(snp) || is.na(trait) || trait == "" || !(trait %in% colnames(data))) next
    message("Processing: ", snp, " → ", trait)
    
    snp_data <- data %>% filter(.data[[snp]] %in% c("0|0", "1|1"))
    
    effect_df <- snp_data %>%
      group_by(Env = Environment, Genotype = .data[[snp]]) %>%
      summarise(mean_trait = mean(.data[[trait]], na.rm = TRUE), .groups = "drop") %>%
      pivot_wider(names_from = Genotype, values_from = mean_trait)
    
    ref_val <- if ("0|0" %in% colnames(effect_df)) effect_df$`0|0` else NA_real_
    alt_val <- if ("1|1" %in% colnames(effect_df)) effect_df$`1|1` else NA_real_
    
    effect_df <- effect_df %>%
      mutate(
        SNP = snp,
        Trait = trait,
        REF = ref_val,
        ALT = alt_val,
        Effect_Size = ((ALT / REF) - 1) * 100
      ) %>%
      select(SNP, Trait, Env, REF, ALT, Effect_Size)
    
    results[[paste(snp, trait, sep = "_")]] <- effect_df
  }
  
  bind_rows(results)
}

all_effects <- compute_all_effects(merged_data, snp_trait_map)
fwrite(all_effects, "Final_Effect_Summary_SpecificSNPs.csv")





#install.packages('R.utils')
library(R.utils)







library(data.table)
library(dplyr)
library(R.utils)

# Load effect data
effect_data <- fread("Final_Effect_Summary_SpecificSNPs.csv")

# Get unique Trait–SNP pairs
trait_snp_list <- unique(effect_data[, .(Trait, SNP)])

# Function to read p-values for a trait
read_trait_pvalues <- function(trait) {
  file <- paste0("mlm_yield_results_", trait, ".y.csv.gz")
  if (!file.exists(file)) return(NULL)
  
  gwas_data <- fread(file)
  gwas_data <- gwas_data %>%
    mutate(SNP = paste0(CHR, "_", POS)) %>%
    select(SNP, P.value = `P`) %>%
    mutate(Trait = trait)
  
  return(gwas_data)
}

# Run for all 7 traits
all_pvals <- bind_rows(lapply(unique(trait_snp_list$Trait), read_trait_pvalues))

# Merge with effect data
effect_annotated <- effect_data %>%
  left_join(all_pvals, by = c("SNP", "Trait")) %>%
  mutate(
    Significance = case_when(
      is.na(P.value) ~ "",
      P.value < 0.001 ~ "***",
      P.value < 0.01  ~ "**",
      P.value < 0.05  ~ "*",
      TRUE ~ ""
    )
  )

fwrite(effect_annotated, "Effect_Summary_With_Pvalues.csv")







read_trait_pvalues <- function(trait) {
  file <- paste0("mlm_yield_results_", trait, ".y.csv.gz")
  if (!file.exists(file)) return(NULL)
  
  gwas_data <- fread(file)
  
  p_col <- paste0(trait, ".y.MLM")  # e.g., earLength.y.MLM
  
  if (!p_col %in% colnames(gwas_data)) {
    warning("Column ", p_col, " not found in ", file)
    return(NULL)
  }
  
  gwas_data <- gwas_data %>%
    select(SNP, P.value = all_of(p_col)) %>%
    mutate(Trait = trait)
  
  return(gwas_data)
}



all_pvals <- bind_rows(lapply(unique(trait_snp_list$Trait), read_trait_pvalues))

effect_annotated <- effect_data %>%
  left_join(all_pvals, by = c("SNP", "Trait")) %>%
  mutate(
    Significance = case_when(
      is.na(P.value) ~ "",
      P.value < 0.001 ~ "***",
      P.value < 0.01  ~ "**",
      P.value < 0.05  ~ "*",
      TRUE ~ ""
    )
  )

fwrite(effect_annotated, "Effect_Summary_With_Pvalues.csv")














library(R.utils)

# List all gzipped CSV files in the current directory
gz_files <- list.files(pattern = "\\.csv\\.gz$")

# Unzip each file ONLY if the corresponding .csv does not already exist
lapply(gz_files, function(f) {
  csv_file <- sub("\\.gz$", "", f)
  if (!file.exists(csv_file)) {
    message("Unzipping: ", f)
    gunzip(f, overwrite = FALSE, remove = FALSE)  # keep .gz file, do not overwrite .csv if it exists
  } else {
    message("Already unzipped: ", csv_file)
  }
})





# Set working directory
setwd("~/Documents/PhD work/Yield_project/SNP effect/SNP_effect_v2/cleaned_data/")

library(data.table)
library(dplyr)
library(stringr)


snp_data <- fread('2023_Missouri_Valley_Medium_cleaned_hundredKernelMass.y_mlm.csv')

# Load full effect summary
effect_data <- fread("Final_Effect_Summary_SpecificSNPs.csv")

#Set the environment you want to process
env <- "Missouri Valley 2023"

#  Clean SNP format once
effect_data$SNP <- paste0("chr", str_replace(effect_data$SNP, "^chr", ""))

# Subset all rows for this environment
env_data <- effect_data[Env == env]

# Get list of traits in this environment
traits <- unique(env_data$Trait)

env_split <- str_split(env, " ", simplify = TRUE)
year <- env_split[, ncol(env_split)]  # Year is the last word
loc <- paste(env_split[, -ncol(env_split)], collapse = "_")  # Everything else is location

# Loop through each trait
for (trait in traits) {
  message("🔄 Processing: ", env, " — ", trait)
  
  # Filter data for this trait
  df <- env_data[Trait == trait]
  
  # Make filename
  file <- paste0(year, "_", loc, "_Medium_cleaned_", trait, ".y_mlm.csv")
  
  # Skip if file missing
  if (!file.exists(file)) {
    message("File not found: ", file)
    next
  }
  
  # Read GWAS file
  p_col <- "value.MLM"
  gwas_data <- fread(file, select = c("SNP", p_col))
  gwas_data$SNP <- as.character(gwas_data$SNP)
  
  # Merge p-values
  df <- df %>%
    left_join(gwas_data, by = "SNP") %>%
    rename(P_value = value.MLM) %>%
    mutate(
      Significance = case_when(
        is.na(P_value) ~ "",
        P_value < 0.001 ~ "***",
        P_value < 0.01  ~ "**",
        P_value < 0.05  ~ "*",
        TRUE ~ ""
      ),
      minus_log10_P = ifelse(is.na(P_value), NA, -log10(P_value))
    )
  
  # Save per-trait result
  out_name <- paste0("Effect_Summary_", str_replace_all(env, " ", "_"), "_", trait, ".csv")
  fwrite(df, out_name)
  message("Done: ", out_name, "\n")
}


# Combine and remove duplicates
all_data <- rbindlist(lapply(files, fread), fill = TRUE)

# Define duplicates based on SNP, Trait, Env
all_data_unique <- all_data[!duplicated(all_data[, .(SNP, Trait, Env)]), ]

# Save cleaned version
fwrite(all_data_unique, "Effect_Summary_ALL_Envs_Traits.csv")











library(data.table)
library(dplyr)
library(ggplot2)

# ---- Load data ----
data <- fread("Effect_Summary_ALL_Envs_Traits.csv")

# ---- Prepare data ----
data_long <- data %>%
  mutate(
    Env_Name = case_when(
      grepl("Ames.*2022", Env) ~ "Ames 2022",
      grepl("Ames.*2023", Env) ~ "Ames 2023",
      grepl("Lincoln.*2022", Env) ~ "Lincoln 2022",
      grepl("Lincoln.*2023", Env) ~ "Lincoln 2023",
      grepl("Missouri.*2022", Env) ~ "Missouri Valley 2022",
      grepl("Missouri.*2023", Env) ~ "Missouri Valley 2023",
      grepl("Crawfordsville.*2022", Env) ~ "Crawfordsville 2022",
      grepl("Crawfordsville.*2023", Env) ~ "Crawfordsville 2023",
      TRUE ~ Env
    ),
    Significance = ifelse(Significance == "", NA, Significance)
  )

env_colors <- c(
  "Ames 2022"             = "#4477AA",  # Blue
  "Ames 2023"             = "#66CCEE",  # Sky Blue
  "Missouri Valley 2022" = "#228833",  # Green
  "Missouri Valley 2023" = "#CCBB44",  # Olive/Yellow
  "Crawfordsville 2022"  = "#EE6677",  # Coral
  "Crawfordsville 2023"  = "#AA3377",  # Plum
  "Lincoln 2022"         = "#BBBBBB",  # Light Gray
  "Lincoln 2023"         = "#000000"   # Black
)

# ---- Loop over each trait ----
for (trait_name in unique(data_long$Trait)) {
  trait_data <- data_long %>% filter(Trait == trait_name)
  
  p <- ggplot(trait_data, aes(x = Effect_Size, y = factor(SNP, levels = unique(SNP)))) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
    
    # Solid circles, filled by environment, no border
    geom_point(
      aes(color = Env_Name),
      shape = 19,
      size = 4
    ) +
    
    # Stars above significant points
    geom_text(
      data = trait_data %>% filter(!is.na(Significance)),
      aes(label = Significance),
      vjust = -1.2,
      fontface = "bold",
      size = 3.5,
      color = "black",
      show.legend = FALSE
    ) +
    
    # Trait name inside plot at top right
    annotate(
      "text",
      x = max(trait_data$Effect_Size, na.rm = TRUE),
      y = Inf,
      label = trait_name,
      hjust = 1.1,
      vjust = 1.5,
      fontface = "bold",
      size = 5
    ) +
    
    scale_color_manual(values = env_colors, name = "Environment") +
    
    labs(
      x = "Effect Size (% Change)",
      y = "SNP"
    ) +
    
    theme_classic(base_size = 14) +
    theme(
      panel.grid.major.y = element_blank(),
      axis.line.y = element_line(color = "grey70"),
      legend.position = "bottom",
      legend.title = element_blank(),  # 👈 Remove "Environment" label
      legend.text = element_text(size = 11)
    ) +
    guides(
      color = guide_legend(override.aes = list(size = 5, shape = 19))
    )
  
  # Save the plot
  ggsave(
    filename = paste0("effect_size_", trait_name, "_clean_final.png"),
    plot = p,
    width = 12, height = 7, dpi = 300
  )
}












library(data.table)
library(dplyr)
library(stringr)
library(ggplot2)

# ---- Load data ----

getwd()
setwd('~/Documents/PhD work/Yield_project/SNP effect/SNP_effect_v2/cleaned_data/')
data <- fread("Effect_Summary_ALL_Envs_Traits.csv")

# ---- Prepare data ----
data_long <- data %>%
  mutate(
    Env_Name = case_when(
      grepl("Ames.*2022", Env) ~ "Ames 2022",
      grepl("Ames.*2023", Env) ~ "Ames 2023",
      grepl("Lincoln.*2022", Env) ~ "Lincoln 2022",
      grepl("Lincoln.*2023", Env) ~ "Lincoln 2023",
      grepl("Missouri.*2022", Env) ~ "Missouri Valley 2022",
      grepl("Missouri.*2023", Env) ~ "Missouri Valley 2023",
      grepl("Crawfordsville.*2022", Env) ~ "Crawfordsville 2022",
      grepl("Crawfordsville.*2023", Env) ~ "Crawfordsville 2023",
      TRUE ~ Env
    ),
    Significance = ifelse(Significance == "", NA, Significance)
  )

# ---- Set colors for environments ----
env_colors <- c(
  "Ames 2022"             = "#4477AA",
  "Ames 2023"             = "#66CCEE",
  "Missouri Valley 2022"  = "#228833",
  "Missouri Valley 2023"  = "#CCBB44",
  "Crawfordsville 2022"   = "#EE6677",
  "Crawfordsville 2023"   = "#AA3377",
  "Lincoln 2022"          = "#BBBBBB",
  "Lincoln 2023"          = "#000000"
)

# ---- Loop over each trait ----
for (trait_name in unique(data_long$Trait)) {
  trait_data <- data_long %>%
    filter(Trait == trait_name) %>%
    mutate(
      SNP_Label = gsub("_", ":", SNP)  # Format SNP as chr#:pos
    )
  
  # Order SNPs by average Effect_Size
  snp_order <- trait_data %>%
    group_by(SNP_Label) %>%
    summarise(avg_eff = mean(abs(Effect_Size), na.rm = TRUE)) %>%
    arrange(avg_eff) %>%
    pull(SNP_Label)
  
  # Plot
  p <- ggplot(trait_data, aes(x = Effect_Size, y = factor(SNP_Label, levels = snp_order))) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
    
    geom_point(aes(color = Env_Name), shape = 19, size = 4) +
    
    geom_text(
      data = trait_data %>% filter(!is.na(Significance)),
      aes(label = Significance),
      vjust = -1.2,
      fontface = "bold",
      size = 3.5,
      color = "black",
      show.legend = FALSE
    ) +
    
    annotate(
      "text",
      x = max(trait_data$Effect_Size, na.rm = TRUE) + 5,
      y = length(unique(trait_data$SNP_Label)) + 1,
      label = trait_name,
      hjust = 1,
      vjust = 0,
      fontface = "bold",
      size = 5
    ) +
    
    scale_color_manual(values = env_colors) +
    
    labs(
      x = "Effect Size (% Change)",
      y = "SNP"
    ) +
    
    theme_classic(base_size = 18) +
    theme(
      panel.grid.major.y = element_blank(),
      axis.line.y = element_line(color = "black", linewidth = 0.7),
      axis.text.y = element_text(
        family = "Courier",  # monospace font
        face = "bold",       # bold text
        size = 14,           # 🔼 increase label size
        color = "black"      # 🔲 force label color to black
      ),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = 11)
    ) +
    guides(color = guide_legend(override.aes = list(size = 5, shape = 19)))
  
  # Save the plot
  ggsave(
    filename = paste0("effect_size_", trait_name, "_clean_final.png"),
    plot = p,
    width = 12,
    height = 7,
    dpi = 300
  )
}







library(data.table)
library(dplyr)
library(stringr)
library(ggplot2)

# ---- Load data ----
data <- fread("Effect_Summary_ALL_Envs_Traits.csv")

# ---- Prepare data ----
data_long <- data %>%
  mutate(
    Env_Name = case_when(
      grepl("Ames.*2022", Env) ~ "Ames 2022",
      grepl("Ames.*2023", Env) ~ "Ames 2023",
      grepl("Lincoln.*2022", Env) ~ "Lincoln 2022",
      grepl("Lincoln.*2023", Env) ~ "Lincoln 2023",
      grepl("Missouri.*2022", Env) ~ "Missouri Valley 2022",
      grepl("Missouri.*2023", Env) ~ "Missouri Valley 2023",
      grepl("Crawfordsville.*2022", Env) ~ "Crawfordsville 2022",
      grepl("Crawfordsville.*2023", Env) ~ "Crawfordsville 2023",
      TRUE ~ Env
    ),
    Significance = ifelse(Significance == "", NA, Significance)
  )

# ---- Environment color scheme ----
env_colors <- c(
  "Ames 2022"             = "#4477AA",
  "Ames 2023"             = "#66CCEE",
  "Missouri Valley 2022"  = "#228833",
  "Missouri Valley 2023"  = "#CCBB44",
  "Crawfordsville 2022"   = "#EE6677",
  "Crawfordsville 2023"   = "#AA3377",
  "Lincoln 2022"          = "#BBBBBB",
  "Lincoln 2023"          = "#000000"
)

# ---- Loop over each trait ----
for (trait_name in unique(data_long$Trait)) {
  trait_data <- data_long %>%
    filter(Trait == trait_name) %>%
    mutate(
      SNP_Label = gsub("_", ":", SNP)  # Format SNP like chr3:123456
    )
  
  # ---- Order SNPs by absolute effect size ----
  snp_order <- trait_data %>%
    group_by(SNP_Label) %>%
    summarise(avg_eff = mean(abs(Effect_Size), na.rm = TRUE)) %>%
    arrange(avg_eff) %>%
    pull(SNP_Label)
  
  # ---- Plot ----
  p <- ggplot(trait_data, aes(x = Effect_Size, y = factor(SNP_Label, levels = snp_order))) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
    
    geom_point(aes(color = Env_Name), shape = 19, size = 4) +
    
    geom_text(
      data = trait_data %>% filter(!is.na(Significance)),
      aes(label = Significance),
      vjust = -1.2,
      fontface = "bold",
      size = 3.5,
      color = "black",
      show.legend = FALSE
    ) +
    
    theme(
      plot.margin = margin(t = 10, r = 10, b = 50, l = 10)  # bottom margin expanded
    ) +
    annotate(
      "text",
      x = 0, y = -Inf,
      label = trait_name,
      hjust = 0, vjust = -6,
      fontface = "bold",
      family = "Courier",
      size = 6,
      color = "black"
    ) +
    
    scale_color_manual(values = env_colors) +
    
    scale_x_continuous(
      breaks = seq(-20, 20, by = 10),
      limits = c(-25, 25)
    ) +
    
    labs(
      x = "Effect Size (% Change)",
      y = "SNP"
    ) +
    
    theme_classic(base_size = 14) +
    theme(
      panel.grid.major.y = element_blank(),
      axis.line.y = element_line(color = "black", linewidth = 0.7),
      
      # Axis labels
      axis.text.y = element_text(family = "Courier", face = "bold", size = 14, color = "black"),
      axis.text.x = element_text(family = "Courier", face = "bold", size = 14, color = "black"),
      
      # Axis titles
      axis.title.y = element_text(family = "Courier", face = "bold", size = 16, color = "black"),
      axis.title.x = element_text(family = "Courier", face = "bold", size = 16, color = "black"),
      
      # Legend
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(family = "Courier", face = "bold", size = 11, color = "black")
    ) +
    guides(color = guide_legend(override.aes = list(size = 5, shape = 19)))
  
  # ---- Save the plot ----
  ggsave(
    filename = paste0("effect_size_", trait_name, "_clean_final.png"),
    plot = p,
    width = 12,
    height = 7,
    dpi = 300
  )
}







library(data.table)
library(dplyr)

# Load your data
df <- fread("filtered_data_withsnp.csv")  # Replace with the correct path

# Make sure the SNP column exists
target_snp <- "chr2_74004678"

# Check if column exists
if (!(target_snp %in% colnames(df))) {
  stop(paste("SNP", target_snp, "not found in the data."))
}

# Count genotype types for the specific SNP, grouped by environment
genotype_counts <- df %>%
  select(year, location, nitrogenTreatment, all_of(target_snp)) %>%
  rename(Genotype = !!target_snp) %>%
  mutate(
    allele_type = case_when(
      Genotype == "0|0" ~ "Reference",
      Genotype == "1|1" ~ "Alternate",
      Genotype %in% c("0|1", "1|0") ~ "Heterozygote",
      is.na(Genotype) | Genotype == "" ~ "Missing",
      TRUE ~ "Other"
    )
  ) %>%
  group_by(year, location, nitrogenTreatment, allele_type) %>%
  summarise(count = n(), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = allele_type, values_from = count, values_fill = 0)

# View results
print(genotype_counts)




library(dplyr)

# Assume your data is in long format with columns: SNP, Environment, effect_size
# Example: effect_data <- read.csv("snp_effects_long.csv")

# Step 1: Compute mean, SD, CV
stability_metrics <- effect_data %>%
  group_by(SNP) %>%
  summarise(
    mean_effect = mean(effect_size, na.rm = TRUE),
    sd_effect = sd(effect_size, na.rm = TRUE),
    cv_effect = sd_effect / abs(mean_effect + 1e-6),  # small offset to avoid division by 0
    .groups = "drop"
  )

# Step 2: Compute sign consistency (directional stability)
sign_consistency <- effect_data %>%
  group_by(SNP) %>%
  summarise(
    median_sign = sign(median(effect_size, na.rm = TRUE)),
    sign_agree = mean(sign(effect_size) == median_sign, na.rm = TRUE),
    .groups = "drop"
  )

# Step 3: Merge and compute stability score
snp_stability <- stability_metrics %>%
  left_join(sign_consistency, by = "SNP") %>%
  mutate(stability_score = sign_agree / (cv_effect + 1))

# Step 4: Optional: Rank SNPs by stability
snp_stability <- snp_stability %>%
  arrange(desc(stability_score))

# View top stable SNPs
head(snp_stability)

# Save output if needed
# write.csv(snp_stability, "snp_stability_metrics.csv", row.names = FALSE)
