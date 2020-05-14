library(readr)
library(corrplot)

source("corrplot_customized.R")

normal_primary <- read_csv("normal_primary.csv")
cell_types <- normal_primary$X1
rownames(normal_primary) <- cell_types
normal_primary <- normal_primary[c(-1, -23, -22, -21, -20)]

normal_primary_lower <- read_csv("normal_primary_lower.csv")
rownames(normal_primary_lower) <- cell_types
normal_primary_lower <- normal_primary_lower[c(-1, -23, -22, -21, -20)]
#normal_primary_lower[is.na(normal_primary_lower)] <- -1

normal_primary_upper <- read_csv("normal_primary_upper.csv")
rownames(normal_primary_upper) <- cell_types
normal_primary_upper <- normal_primary_upper[c(-1, -23, -22, -21, -20)]
#normal_primary_upper[is.na(normal_primary_upper)] <- 1

normal_primary_p <- read_csv("normal_primary_p.csv")
rownames(normal_primary_p) <- cell_types
normal_primary_p <- normal_primary_p[c(-1, -23, -22, -21, -20)]

#normal_primary[is.na(normal_primary)] = -1
#normal_primary_lower[is.na(normal_primary_lower)] = -1
#normal_primary_upper[is.na(normal_primary_upper)] = -1

corr_matrix <- as.matrix(normal_primary)
rownames(corr_matrix) <- cell_types

png("normal_primary.png", width = 8, height=10, unit = 'in', res = 400)

corrplot(corr_matrix, 
         low = as.matrix(normal_primary_lower), 
         upp = as.matrix(normal_primary_upper),
         #p.mat = as.matrix(normal_primary_p),
         plotCI = 'rect',
         col=colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                                "#F0F0F0", 
                                "#92C5DE", "#4393C3", "#2166AC", "#053061"))(500)
         )

dev.off()
