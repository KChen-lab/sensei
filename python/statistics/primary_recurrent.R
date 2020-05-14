library(readr)
library(corrplot)

source("corrplot_customized.R")

lgg <- read_csv("lgg.csv")
cell_types <- lgg$X1
rownames(lgg) <- cell_types
lgg <- lgg[-1]

gbm <- read_csv("gbm.csv")
cell_types <- gbm$X1
rownames(gbm) <- cell_types
gbm <- gbm[-1]

ov <- read_csv("ov.csv")
cell_types <- ov$X1
rownames(ov) <- cell_types
ov <- ov[-1]

corr_matrix = matrix(nrow = length(cell_types), ncol = 3)

rownames(corr_matrix) <- cell_types
colnames(corr_matrix) <- c("LGG", "GBM", "OV")

low_matrix = corr_matrix
upp_matrix = corr_matrix

corr_matrix[, "LGG"] = lgg$Correlation
low_matrix[, "LGG"] = lgg$`95% CI Lower`
upp_matrix[, "LGG"] = lgg$`95% CI Upper`

corr_matrix[, "GBM"] = gbm$Correlation
low_matrix[, "GBM"] = gbm$`95% CI Lower`
upp_matrix[, "GBM"] = gbm$`95% CI Upper`

corr_matrix[, "OV"] = ov$Correlation
low_matrix[, "OV"] = ov$`95% CI Lower`
upp_matrix[, "OV"] = ov$`95% CI Upper`

png("primary_recurrent.png", width = 4, height=10, unit = 'in', res = 400)

corrplot(corr_matrix, 
         low = low_matrix, 
         upp = upp_matrix,
         #p.mat = as.matrix(normal_primary_p),
         plotCI = 'rect', cl.pos="n",
         col=colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                                "#F0F0F0", 
                                "#92C5DE", "#4393C3", "#2166AC", "#053061"))(500)
         )

dev.off()
