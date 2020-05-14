library(readr)
library(corrplot)

source("corrplot_customized.R")

res <- read_csv("GSE120575-res.csv")
cell_types <- res$X1
rownames(res) <- cell_types
res <- res[-1]

corr_matrix = matrix(nrow = length(cell_types), ncol = 1)

rownames(corr_matrix) <- cell_types
colnames(corr_matrix) <- c("GSE120575")

low_matrix = corr_matrix
upp_matrix = corr_matrix

corr_matrix[, "GSE120575"] = res$Correlation
low_matrix[, "GSE120575"] = res$`95% CI Lower`
upp_matrix[, "GSE120575"] = res$`95% CI Upper`

png("melanoma.png", width = 4, height=10, unit = 'in', res = 400)

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
