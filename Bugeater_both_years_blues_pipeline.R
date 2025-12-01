library(data.table)
library(dplyr)
library(tidyr)
library(SpATS)
library(tibble)
library(ggplot2)

getwd()

setwd('~/Documents/PhD work/Projects/Yield_project/')


df_2021 <- fread('20_21_yield_extreme.csv')


unique(df_2021[year == 2021, genotype])



# Read and clean the data
df <- fread("20_21_yield_extreme.csv")
names(df)[names(df) == "Ear Width CM"] <- "earWidth"
names(df)[names(df) == "Ear Length CM"] <- "earLength"
names(df)[names(df) == "Kernels per row"] <- "kernelsPerRow"
names(df)[names(df) == "Kernel row number"] <- "kernelRowNumber"
names(df)[names(df) == "ear weight (grams)"] <- "earWeight"
names(df)[names(df) == "plot weight (grams)"] <- "plotWeight"
names(df)[names(df) == "100 K weight"] <- "hundredKernelMass"

traits <- c("earWidth", "earLength", "kernelsPerRow", "kernelRowNumber",
            "earWeight", "plotWeight", "hundredKernelMass")
years <- sort(unique(df$year))


df_long <- df %>%
  pivot_longer(cols = all_of(traits), names_to = "trait", values_to = "value")

for (yr in years) {
  p <- ggplot(df_long %>% filter(year == yr), aes(x = value)) +
    geom_histogram(bins = 40, fill = ifelse(yr == 2020, "skyblue", "tomato"), color = "black") +
    facet_wrap(~ trait, scales = "free") +
    theme_bw() +
    labs(title = paste("Raw Data Distribution by Trait (", yr, ")", sep=""))
  print(p)
  # Optionally save: ggsave(paste0("raw_histogram_", yr, ".png"), p)
}


thresholds <- list(
  "2020" = list(kernelsPerRow = c(NA, NA)
                # plotWeight = c(10, 150)
  ),
  "2021" = list(kernelsPerRow = c(10, NA),
                earWeight = c(NA, 98),
                plotWeight = c(NA, 80),
                hundredKernelMass = c(13, NA)
  )
)

df_masked <- copy(df)  # Deep copy for safety

for (yr in names(thresholds)) {
  for (trait in names(thresholds[[yr]])) {
    lower <- thresholds[[yr]][[trait]][1]
    upper <- thresholds[[yr]][[trait]][2]
    idx <- df_masked$year == as.integer(yr) &
      !is.na(df_masked[[trait]]) &
      (df_masked[[trait]] < lower | df_masked[[trait]] > upper)
    df_masked[[trait]][idx] <- NA
  }
}
fwrite(df_masked, "raw_data_masked_outliers.csv")



df_long_masked <- df_masked %>%
  pivot_longer(cols = all_of(traits), names_to = "trait", values_to = "value")

for (yr in years) {
  p <- ggplot(df_long_masked %>% filter(year == yr), aes(x = value)) +
    geom_histogram(bins = 40, fill = ifelse(yr == 2020, "skyblue", "tomato"), color = "black") +
    facet_wrap(~ trait, scales = "free") +
    theme_bw() +
    labs(title = paste("Masked Data Distribution by Trait (", yr, ")", sep=""))
  print(p)
  # Optionally save: ggsave(paste0("masked_histogram_", yr, ".png"), p)
}




all_blups <- data.frame()
all_blues <- data.frame()

for (yr in years) {
  df_year <- df_masked %>% filter(year == yr)
  for (trait in traits) {
    if (!(trait %in% colnames(df_year))) next
    df_trait <- df_year %>% filter(!is.na(.data[[trait]]))
    if (nrow(df_trait) < 2) next
    
    rowKnots <- floor(max(df_trait$Row, na.rm=TRUE) / 2) + 1
    colKnots <- floor(max(df_trait$Column, na.rm=TRUE) / 2) + 1
    
    model <- SpATS(
      response = trait,
      genotype = "PlotID",
      genotype.as.random = TRUE,
      spatial = ~ SAP(Column, Row, nseg = c(colKnots, rowKnots)),
      data = df_trait
    )
    
    blup_df <- as_tibble(model$coeff, rownames = "PlotID") %>%
      mutate(
        year = yr,
        trait = trait,
        PlotID = as.character(PlotID),
        value = value + model$coeff["Intercept"]
      ) %>%
      left_join(
        df_trait %>% mutate(PlotID = as.character(PlotID)) %>% 
          select(PlotID, genotype) %>% distinct(),
        by = "PlotID"
      ) %>%
      select(year, trait, PlotID, genotype, value)
    all_blups <- bind_rows(all_blups, blup_df)
    
    mod <- lm(value ~ genotype, data = blup_df)
    cf <- coef(mod)
    intercept <- cf[1]
    geno_coefs <- cf[-1]
    genonames <- sub("genotype", "", names(geno_coefs))
    blue_values <- intercept + geno_coefs
    ref_genotype <- as.character(sort(unique(blup_df$genotype))[1])
    blues_df <- data.frame(
      year = yr,
      trait = trait,
      genotype = c(ref_genotype, genonames),
      BLUE = c(intercept, as.vector(blue_values))
    )
    all_blues <- bind_rows(all_blues, blues_df)
  }
}



for (yr in years) {
  # BLUPs
  this_blups <- all_blups %>% filter(year == yr)
  this_blups_wide <- this_blups %>%
    pivot_wider(
      id_cols = c(PlotID, genotype),
      names_from = trait,
      values_from = value
    )
  write.csv(this_blups_wide, paste0("SpATS_BLUPs_WIDE_", yr, ".csv"), row.names = FALSE)
  
  # BLUEs
  this_blues <- all_blues %>% filter(year == yr)
  this_blues_wide <- this_blues %>%
    pivot_wider(
      id_cols = genotype,
      names_from = trait,
      values_from = BLUE
    )
  write.csv(this_blues_wide, paste0("BLUEs_WIDE_", yr, ".csv"), row.names = FALSE)
}



library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)

traits <- c("earWidth", "earLength", "kernelsPerRow", "kernelRowNumber",
            "earWeight", "plotWeight", "hundredKernelMass")

# Combine both years into one dataframe, add a year column if not already present
blues_2020 <- fread("BLUEs_WIDE_2020.csv") %>% mutate(year = 2020)
blues_2021 <- fread("BLUEs_WIDE_2021.csv") %>% mutate(year = 2021)
blues_all <- bind_rows(blues_2020, blues_2021)



blues_long <- blues_all %>%
  pivot_longer(cols = all_of(traits), names_to = "trait", values_to = "BLUE")

for (yr in unique(blues_long$year)) {
  p <- ggplot(blues_long %>% filter(year == yr), aes(x = BLUE)) +
    geom_histogram(bins = 40, fill = ifelse(yr == 2020, "skyblue", "tomato"), color = "black") +
    facet_wrap(~ trait, scales = "free") +
    theme_bw() +
    labs(title = paste("BLUEs Distribution by Trait (", yr, ")", sep=""))
  print(p)
  # ggsave(paste0("BLUEs_raw_histogram_", yr, ".png"), p)
}




# Example thresholds, set your own!
thresholds_blues <- list(
  "2020" = list(
    earLength = c(11, NA),
    earWidth = c(3.15, NA),
    kernelRowNumber = c(13.88, 14.1),
    earWeight = c(44, 92),
    hundredKernelMass = c(17, NA),
    kernelsPerRow = c(15.5, 30)
  ),
  "2021" = list(
    earWidth = c(NA, 4),
    kernelsPerRow = c(19.5, 31),
    earWeight = c(26, 66),
    plotWeight = c(NA, 56),
    hundredKernelMass = c(NA, 29)
  )
)

blues_masked <- copy(blues_all)  # Deep copy

for (yr in names(thresholds_blues)) {
  for (trait in names(thresholds_blues[[yr]])) {
    lower <- thresholds_blues[[yr]][[trait]][1]
    upper <- thresholds_blues[[yr]][[trait]][2]
    idx <- blues_masked$year == as.integer(yr) &
      !is.na(blues_masked[[trait]]) &
      (blues_masked[[trait]] < lower | blues_masked[[trait]] > upper)
    blues_masked[[trait]][idx] <- NA
  }
}




blues_long_masked <- blues_masked %>%
  pivot_longer(cols = all_of(traits), names_to = "trait", values_to = "BLUE")

for (yr in unique(blues_long_masked$year)) {
  p <- ggplot(blues_long_masked %>% filter(year == yr), aes(x = BLUE)) +
    geom_histogram(bins = 40, fill = ifelse(yr == 2020, "skyblue", "tomato"), color = "black") +
    facet_wrap(~ trait, scales = "free") +
    theme_bw() +
    labs(title = paste("Masked BLUEs Distribution by Trait (", yr, ")", sep=""))
  print(p)
  # ggsave(paste0("BLUEs_masked_histogram_", yr, ".png"), p)
}


fwrite(blues_masked %>% filter(year == 2020) %>% select(-year), "BLUEs_WIDE_2020_no_outliers.csv")
fwrite(blues_masked %>% filter(year == 2021) %>% select(-year), "BLUEs_WIDE_2021_no_outliers.csv")







library(data.table)

blues_2020 <- fread("BLUEs_WIDE_2020_no_outliers.csv")
blues_2021 <- fread("BLUEs_WIDE_2021_no_outliers.csv")



geno_2020 <- blues_2020$genotype
geno_2021 <- blues_2021$genotype



common_genos <- intersect(geno_2020, geno_2021)
length(common_genos)  # How many?
head(common_genos)    # Preview a few








vcf_ids <- fread("798_ids.csv")
pheno_ids <- fread("BLUEs_WIDE_2020_no_outliers.csv")

#head(fread("sample_ids.csv"))
#head(fread("BLUEs_WIDE_2020_no_outliers.csv"))




# Extract the genotype vectors
vcf_ids <- vcf_ids$genotype
pheno_ids <- pheno_ids$genotype

# Genotypes in both
common_ids <- intersect(vcf_ids, pheno_ids)

# Genotypes in VCF but not in phenotype
vcf_only <- setdiff(vcf_ids, pheno_ids)

# Genotypes in phenotype but not in VCF
pheno_only <- setdiff(pheno_ids, vcf_ids)


fwrite(data.table(genotype = common_ids), "common_genotypes_21.csv")
fwrite(data.table(genotype = vcf_only), "vcf_only_genotypes_21.csv")
fwrite(data.table(genotype = pheno_only), "pheno_only_genotypes_21.csv")

















