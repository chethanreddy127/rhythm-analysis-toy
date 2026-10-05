library(ggplot2)
library(dplyr)
devtools::install_github('gangwug/MetaCycle')
library(MetaCycle)
install.packages('rain', repos = c('https://bioc.r-universe.dev', 'https://cloud.r-project.org'))
library(rain)

### Creating ground truth dataset ----

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

#Generate acrophase values only for cycling genes
geneList$acrophase[isRhythmic] <- runif(
  sum(isRhythmic),
  min = 0,
  max = 24
)
# Creating Time Series data ----

timeSeries = data.frame(gene_name = geneList$gene_name, matrix(nrow = 1000, ncol = 12))
x = 4
sample_by_x = seq(4, 48, by = x)
colnames(timeSeries)[2:13] = paste0("t", sample_by_x)

generateTimeData <- function(gene_index) {
  time_points = sample_by_x
  
  
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

### Creating nested models and performing significance testing ----


metaTest <- meta2d(infile = "MetaCycle Analysis", 
                   filestyle = "csv", 
                   timepoints = sample_by_x, 
                   outputFile = FALSE, 
                   inDF = timeSeries
                   )
JTK.results <- metaTest$JTK[order(metaTest$JTK$BH.Q),
                            ,
                            drop = FALSE
                            ]
rainTest <- rain(t(timeSeries[, -1]), 
                 deltat = x, period = 24, 
                 adjp.method = "BH", 
                 peak.border = c(0.2, 0.8), 
                 verbose = TRUE
                 )


# Create a table comparing the cycling transcripts from MetaCycle to the ground truth dataset
JTK.results$rhythmic[JTK.results$BH.Q < 0.3] = "y" 
JTK.results$rhythmic[JTK.results$BH.Q >= 0.3] = "n"

compare.meta <- merge(JTK.results, geneList, by.x = "CycID", by.y = "gene_name", all.x = TRUE) %>%
  dplyr::select(CycID, BH.Q, rhythmic.x, rhythmic.y, PER, LAG, AMP, meanExpr, amplitude, acrophase) %>%
  dplyr::rename(JTK.rhythmic = rhythmic.x, ground_truth_rhythmic = rhythmic.y)

compare.meta <- compare.meta %>%
  dplyr::mutate(
    JTK.rhythmic = ifelse(BH.Q < 0.3, "y", "n"),
    correct = dplyr::case_when(
      JTK.rhythmic == "y" & ground_truth_rhythmic == "y" ~ "true positive",
      JTK.rhythmic == "y" & ground_truth_rhythmic == "n" ~ "false positive",
      JTK.rhythmic == "n" & ground_truth_rhythmic == "y" ~ "false negative",
      JTK.rhythmic == "n" & ground_truth_rhythmic == "n" ~ "true negative"
    )
  )

sum(compare.meta$correct == "true positive")
sum(compare.meta$correct == "false positive")
sum(compare.meta$correct == "false negative")
sum(compare.meta$correct == "true negative")

