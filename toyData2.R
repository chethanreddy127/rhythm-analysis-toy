library(ggplot2)


# Creating ground truth dataset

n_genes <- 1000

geneList <- data.frame(
  gene_name = paste0("gene", seq_len(n_genes)),
  rhythmic = sample(c("y", "n"), n_genes, replace = TRUE),
  meanExpr = rnorm(n_genes, mean = 5, sd = 2),
  amplitude = NA_real_,
  acrophase = NA_real_
)

# Restrict mean expression values to 0 to 10
wrongMean <- geneList$meanExpr < 0 | geneList$meanExpr > 10

while (any(wrongMean)) {
  
  # Regenerate only  invalid values
  geneList$meanExpr[wrongMean] <- rnorm(
    sum(wrongMean),
    mean = 5,
    sd = 2
  )
  
  # Update wrongMean
  wrongMean <- geneList$meanExpr < 0 | geneList$meanExpr > 10
}

# Generate amplitudes only for cycling genes
isRhythmic <- geneList$rhythmic == "y"

geneList$amplitude[isRhythmic] <- rnorm(
  sum(isRhythmic),
  mean = 2.5,
  sd = 1
)

# Restrict cycling amplitudes to the range 0 to 5
wrongAmplitude <- isRhythmic &
  (
    geneList$amplitude < 0 |
      geneList$amplitude > 5
  )

while (any(wrongAmplitude, na.rm = TRUE)) {

  # Regenerate invalid amplitudes
  geneList$amplitude[wrongAmplitude] <- rnorm(
    sum(wrongAmplitude, na.rm = TRUE),
    mean = 2.5,
    sd = 1
  )
  
  # Update wrongAmplitude
  wrongAmplitude <- isRhythmic &
    (
      geneList$amplitude < 0 |
        geneList$amplitude > 5
    )
}

# Generate acrophase values only for cycling genes
geneList$acrophase[isRhythmic] <- runif(
  sum(isRhythmic),
  min = 0,
  max = 24
)
timeSeries = data.frame(gene_name = geneList$gene_name, matrix(nrow = 1000, ncol = 12))
sample_by_4 = seq(4, 48, by = 4)
colnames(timeSeries)[2:13] = paste0("t", sample_by_4)

generateTimeData <- function(gene_index) {
  time_points = sample_by_4
  
  
  MESOR <- geneList$meanExpr[gene_index]
  Amp <- geneList$amplitude[gene_index]
  Phase <- geneList$acrophase[gene_index]
  noise = rnorm(length(time_points), 0, 1)
  
  if (geneList$rhythmic[gene_index] == "y") {
    measurements = MESOR + Amp * cos((2*pi/24) * time_points + Phase) + noise
  } else {
    measurements = MESOR + noise
  }
  
  return(measurements)
}
generated_data <- t(apply(matrix(seq_len(nrow(timeSeries)), ncol = 1), 1, FUN = generateTimeData))
timeSeries[, -1] <- generated_data

rhythmicGenes = data.frame(
  gene_name = geneList$gene_name,
  pValue = rep(0,1000),
  BHpVal = rep(0,1000),
  AmpMESORratio = rep (0,1000),
  acrophase = rep(0,1000),
  rhythmic = rep(NA, 1000)
)

# for each gene, create an array of time and expression values and create reduced and full models, then LRF, store p val in rhythmicGenes table

sigTest <- function(gene_index) {
  
  geneName <- geneList$gene_name[gene_index]
  time <- sample_by_4
  expression <- as.numeric(timeSeries[gene_index, -1])
  
  expressionTrig <- data.frame(
    expression = expression,
    time = time,
    cos24 = cos((2 * pi / 24) * time),
    sin24 = sin((2 * pi / 24) * time)
  )
  
  fullModel <- lm(expression ~ cos24 + sin24, data = expressionTrig)
  reducedModel <- lm(expression ~ 1, data = expressionTrig)
  
  lr_result <- lmtest::lrtest(reducedModel, fullModel)
  
  coefficients <- coef(fullModel)
  
  beta0 <- coefficients["(Intercept)"]
  beta_cos <- coefficients["cos24"]
  beta_sin <- coefficients["sin24"]
  
  amplitude <- sqrt(beta_cos^2 + beta_sin^2)
  phase_rad <- atan2(-beta_sin, beta_cos)
  acrophase <- (phase_rad %% (2 * pi)) / (2 * pi) * 24
  
   data.frame(
    gene_name = geneName,
    pValue = lr_result$`Pr(>Chisq)`[2],
    AmpMESORratio = amplitude / beta0,
    acrophase = acrophase
  )
}

results <- lapply(
  seq_len(nrow(timeSeries)),
  sigTest
)



rhythmicGenes <- do.call(
  rbind,
  results
)

rhythmicGenes$BHpVal <- p.adjust(
  rhythmicGenes$pValue,
  method = "BH"
)

rhythmicGenes$rhythmic <- ifelse(
  rhythmicGenes$BHpVal < 0.1 & rhythmicGenes$AmpMESORratio > 0.2,
  "y",
  "n"
)

compare <- merge(
  rhythmicGenes[, c("gene_name", "rhythmic")],
  geneList[, c("gene_name", "rhythmic")],
  by = "gene_name",
  all = TRUE
)

colnames(compare) <- c(
  "gene_name",
  "empirical_rhythmic",
  "ground_truth_rhythmic"
)

agree = subset(compare, empirical_rhythmic == ground_truth_rhythmic)
disagree = subset(compare, empirical_rhythmic != ground_truth_rhythmic)

# ------------------------------------------------------------
# 5. Plot gene groups
# ------------------------------------------------------------

first_rhythmic <- geneList$gene_name[
  geneList$rhythmic == "y"
][seq_len(min(5, sum(geneList$rhythmic == "y")))]

first_non_rhythmic <- geneList$gene_name[
  geneList$rhythmic == "n"
][seq_len(min(5, sum(geneList$rhythmic == "n")))]

disagreement_genes <- compare$gene_name[
  compare$empirical_rhythmic != compare$ground_truth_rhythmic
]

first_disagreements <- disagreement_genes[
  seq_len(min(3, length(disagreement_genes)))
]

create_plot_data <- function(selected_genes) {
  plot_data <- lapply(selected_genes, function(gene_name) {
    
    time_series_index <- match(gene_name, timeSeries$gene_name)
    gene_list_index <- match(gene_name, geneList$gene_name)
    rhythmic_index <- match(gene_name, rhythmicGenes$gene_name)
    
    observed_time <- sample_by_4
    observed_expression <- as.numeric(timeSeries[time_series_index, -1])
    
    expression_data <- data.frame(
      time = observed_time,
      expression = observed_expression,
      cos24 = cos((2 * pi / 24) * observed_time),
      sin24 = sin((2 * pi / 24) * observed_time)
    )
    
    full_model <- lm(expression ~ cos24 + sin24, data = expression_data)
    reduced_model <- lm(expression ~ 1, data = expression_data)
    
    empirical_classification <- rhythmicGenes$rhythmic[rhythmic_index]
    
    fitted_time <- seq(
      min(sample_by_4),
      max(sample_by_4),
      length.out = 200
    )
    
    if (empirical_classification == "y") {
      full_coefficients <- coef(full_model)
      beta0 <- unname(full_coefficients["(Intercept)"])
      beta_cos <- unname(full_coefficients["cos24"])
      beta_sin <- unname(full_coefficients["sin24"])
      
      empirical_curve <- beta0 +
        beta_cos * cos((2 * pi / 24) * fitted_time) +
        beta_sin * sin((2 * pi / 24) * fitted_time)
      
      empirical_label <- "Empirical full cosinor model"
    } else {
      reduced_mean <- unname(coef(reduced_model)["(Intercept)"])
      empirical_curve <- rep(reduced_mean, length(fitted_time))
      empirical_label <- "Empirical reduced model"
    }
    
    true_mesor <- geneList$meanExpr[gene_list_index]
    true_amplitude <- geneList$amplitude[gene_list_index]
    true_acrophase <- geneList$acrophase[gene_list_index]
    
    if (geneList$rhythmic[gene_list_index] == "y") {
      ground_truth_curve <- true_mesor +
        true_amplitude *
        cos((2 * pi / 24) * fitted_time + true_acrophase)
    } else {
      ground_truth_curve <- rep(true_mesor, length(fitted_time))
    }
    
    observed_data <- data.frame(
      gene_name = gene_name,
      time = observed_time,
      value = observed_expression,
      curve = "Observed expression"
    )
    
    ground_truth_data <- data.frame(
      gene_name = gene_name,
      time = fitted_time,
      value = ground_truth_curve,
      curve = "Ground-truth model"
    )
    
    empirical_data_plot <- data.frame(
      gene_name = gene_name,
      time = fitted_time,
      value = empirical_curve,
      curve = empirical_label
    )
    
    rbind(observed_data, ground_truth_data, empirical_data_plot)
  })
  
  do.call(rbind, plot_data)
}

plot_gene_group <- function(selected_genes, plot_title) {
  if (length(selected_genes) == 0) {
    return(
      ggplot() +
        theme_void() +
        labs(title = paste(plot_title, "- no genes available"))
    )
  }
  
  plot_data <- create_plot_data(selected_genes)
  
  ggplot(
    plot_data,
    aes(x = time, y = value, color = curve, group = curve)
  ) +
    geom_point(
      data = subset(plot_data, curve == "Observed expression"),
      size = 2,
      alpha = 0.8
    ) +
    geom_line(
      data = subset(plot_data, curve != "Observed expression"),
      linewidth = 1
    ) +
    facet_wrap(~ gene_name, scales = "free_y") +
    scale_color_manual(
      values = c(
        "Observed expression" = "black",
        "Ground-truth model" = "red",
        "Empirical full cosinor model" = "blue",
        "Empirical reduced model" = "blue"
      )
    ) +
    labs(
      title = plot_title,
      x = "Time",
      y = "Expression",
      color = NULL
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = "bottom"
    )
}

# Plot 1: first 5 ground-truth rhythmic genes
rhythmic_plot <- plot_gene_group(
  first_rhythmic,
  "First Five Ground-Truth Rhythmic Genes"
)
rhythmic_plot

# Plot 2: first 5 ground-truth non-rhythmic genes
non_rhythmic_plot <- plot_gene_group(
  first_non_rhythmic,
  "First Five Ground-Truth Non-Rhythmic Genes"
)
non_rhythmic_plot

# Plot 3: first 3 genes where empirical and ground truth disagree
disagreement_plot <- plot_gene_group(
  first_disagreements,
  "First Three Empirical/Ground-Truth Disagreement Genes"
)
disagreement_plot

