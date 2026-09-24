library(ggplot2)
library(lmtest)

# Creating ground truth dataset ----

geneList = data.frame(gene = rep(NA, 1000), cycling = rep(NA, 1000), meanExpr = rep(NA, 1000), amplitude = rep(NA, 1000), acrophase = rep(NA, 1000))
cycling = c("y", "n")
for (i in 1:nrow(geneList)) {
  geneList[i,1] = paste0("gene", i)
  geneList[i,2] = sample(cycling, 1, replace = TRUE, prob = NULL)
  geneList[i,3] = rnorm(1, 5, 2)
  while (geneList[i,3] < 0 | geneList[i,3] > 10){
    geneList[i,3] = rnorm(1, 5, 2)
  }
  if (geneList[i,2] == "y") {
    geneList[i,4] = rnorm(1, 2.5, 1)
    while (geneList[i,4] < 0 | geneList [i,4] > 5) {
      geneList[i,4] = rnorm(1, 2.5, 1)
    }
    geneList[i,5] = runif(1, 0, 24)
  }
}
timeSeries = data.frame(gene = geneList$gene, matrix(NA, nrow = 1000, ncol = 48))
colnames(timeSeries)[2:49] = paste0("t", 0:47)


for (i in 1:nrow(timeSeries)) {
  isCycling = geneList[i, 2] == "y"
  meanExpr = geneList[i, 3]
  amplitude = geneList[i, 4]
  acrophase = geneList[i, 5]
  
  for (j in 2:50) {
    t = j - 2
    noise = rnorm(1, 0, 1)
    
    if (isCycling) {
      timeSeries[i, j] = meanExpr + amplitude * cos((2*pi/24) * t + acrophase) + noise
    } else {
      timeSeries[i, j] = meanExpr + noise
    }
  }
}

print(geneList)
print(timeSeries)

#Plotting the first few rhythmic and non-rhythmic genes----

# Pick first 3 cycling and first 3 non-cycling genes
cycling_gene_idx <- which(geneList$cycling == "y")[1:3]
non_cycling_gene_idx <- which(geneList$cycling == "n")[1:3]

# Function to prepare data for plotting
prepare_plot_data <- function(gene_idx, geneList, timeSeries) {
  times <- 0:47
  measured <- as.numeric(timeSeries[gene_idx, 2:49])
  
  meanExpr <- geneList[gene_idx, 3]
  amplitude <- geneList[gene_idx, 4]
  acrophase <- geneList[gene_idx, 5]
  isCycling <- geneList[gene_idx, 2] == "y"
  
  # Calculate ideal values
  if (isCycling) {
    ideal <- meanExpr + amplitude * cos((2*pi/24) * times + acrophase)
  } else {
    ideal <- rep(meanExpr, length(times))
  }
  
  data.frame(
    time = times,
    measured = measured,
    ideal = ideal,
    gene = geneList[gene_idx, 1],
    cycling = geneList[gene_idx, 2]
  )
}

# Prepare data for all 6 genes
plot_data <- data.frame()
for (idx in c(cycling_gene_idx, non_cycling_gene_idx)) {
  plot_data <- rbind(plot_data, prepare_plot_data(idx, geneList, timeSeries))
}

# Create plots
ggplot(plot_data, aes(x = time, y = measured)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_line(aes(y = ideal), color = "red", size = 1) +
  facet_wrap(~gene, scales = "free_y") +
  labs(x = "Time (hours)", y = "Expression Level", title = "Gene Expression Time Series") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# Assuming I know nothing about the cycling nature of genes, I want to model each one with a linear model----

# Screen for rhythmicity
geneData = list()
cosinorModels = list()
nullModels = list()
lrTest = list()
rhythmicGenes = data.frame(
  gene_name = geneList$gene,
  pValue = rep(0,1000),
  BHpVal = rep(0,1000),
  AmpMESORratio = rep (0,1000),
  acrophase = rep(0,1000),
  rhythmic = rep(NA, 1000)
  )
for (i in 1:nrow(timeSeries)) {
  expression_values = as.numeric(timeSeries[i, 2:49])
  gene_name = paste0("gene", i)
  geneData[[gene_name]] = data.frame(
    time = 0:47,
    expression =  expression_values,
    cos24 = 0:47,
    sin24 = 0:47
  )
  geneData[[gene_name]]$cos24 <- cos((2 * pi / 24) * geneData[[gene_name]]$time)
  geneData[[gene_name]]$sin24 <- sin((2 * pi / 24) * geneData[[gene_name]]$time)
  cosinorModels[[gene_name]] = lm(expression ~ cos24 + sin24, data = geneData[[gene_name]])
  
  ## ~ 1 or ~ time?
  nullModels[[gene_name]] = lm(expression ~ 1, data = geneData[[gene_name]])
  
  #Significance testing and coefficient derivation
  lrTest[[gene_name]] = lrtest(nullModels[[gene_name]], cosinorModels[[gene_name]])
  rhythmicGenes$pValue[i] = lrTest[[gene_name]]$"Pr(>Chisq)"[2] 
  beta0 = coef(cosinorModels[[gene_name]])[1]
  beta_cos = coef(cosinorModels[[gene_name]])["cos24"]
  beta_sin = coef(cosinorModels[[gene_name]])["sin24"]
  Amp = sqrt((beta_cos)^2 +(beta_sin)^2)
  Phase = atan2(-beta_sin, beta_cos)
  MESOR = beta0
  rhythmicGenes$AmpMESORratio[i] = Amp/MESOR
  rhythmicGenes$acrophase[i] = Phase
}
q = 0.1

rhythmicGenes[order(rhythmicGenes$pValue), ]
rhythmicGenes$BHpVal <- p.adjust(
  rhythmicGenes$pValue,
  method = "BH"
)
for (i in 1:nrow(rhythmicGenes)) {
  if (rhythmicGenes$BHpVal[i] < q && rhythmicGenes$AmpMESORratio[i] > 0.2) {
    rhythmicGenes$rhythmic[i] = "y"
  } else {
    rhythmicGenes$rhythmic[i] = "n"
    rhythmicGenes$AmpMESORratio[i] = NA
    rhythmicGenes$acrophase[i] = NA
    
  }
}

#Compare my empirically determined rhythmic genes vs ground truth
compare <- merge(
  rhythmicGenes[, c("gene_name", "rhythmic")],
  geneList[, c("gene", "cycling")],
  by.x = "gene_name",
  by.y = "gene",
  all = TRUE
)

colnames(compare) <- c(
  "gene",
  "empirical_rhythmic",
  "ground_truth_rhythmic"
)

agree = subset(compare, empirical_rhythmic == ground_truth_rhythmic)
disagree = subset(compare, empirical_rhythmic != ground_truth_rhythmic)
