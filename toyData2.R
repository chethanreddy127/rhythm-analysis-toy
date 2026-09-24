# Creating ground truth dataset

n_genes <- 1000

geneList <- data.frame(
  gene_name = paste0("gene", seq_len(n_genes)),
  cycling = sample(c("y", "n"), n_genes, replace = TRUE),
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
isCycling <- geneList$cycling == "y"

geneList$amplitude[isCycling] <- rnorm(
  sum(isCycling),
  mean = 2.5,
  sd = 1
)

# Restrict cycling amplitudes to the range 0 to 5
wrongAmplitude <- isCycling &
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
  wrongAmplitude <- isCycling &
    (
      geneList$amplitude < 0 |
        geneList$amplitude > 5
    )
}

# Generate acrophase values only for cycling genes
geneList$acrophase[isCycling] <- runif(
  sum(isCycling),
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
  
  if (geneList$cycling[gene_index] == "y") {
    measurements = MESOR + Amp * cos((2*pi/24) * time_points + Phase) + noise
  } else {
    measurements = MESOR + noise
  }
  
  return(measurements)
}
generated_data <- t(apply(matrix(seq_len(nrow(timeSeries)), ncol = 1), 1, FUN = generateTimeData))
timeSeries[, -1] <- generated_data

rhythmicGenes = data.frame(
  gene_name = geneList$gene,
  pValue = rep(0,1000),
  BHpVal = rep(0,1000),
  AmpMESORratio = rep (0,1000),
  acrophase = rep(0,1000),
  rhythmic = rep(NA, 1000)
)








