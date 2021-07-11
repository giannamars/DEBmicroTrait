country <- c('Canada','France','Germany', 'Italy', 'Japan', 'Russia', 'UK', 'USA')
chemistry <- c(4,8,24,1,6,4,23,51)
economics <- c(3,3,1,1,0,3,6,43)
literature <- c(2,11,8,6,2,5,7,8)
medicine <- c(4,12,18,5,3,2,26,70)
peace <- c(1,10,5,1,1,3,11,19)
physics <- c(4,9,24,5,11,10,20,66)

nobel_data <- data.frame(chemistry, economics, literature, medicine, peace, physics)
rownames(nobel_data) <- country
nobel_data$rsum <- rowSums(nobel_data[,])
nobel_data <- rbind(nobel_data, colSums(nobel_data[,]))
rownames(nobel_data)[9] <- 'csum'

nobel_data[, 1:6]
nobel_data_rprofile <- nobel_data[, 1:6]/nobel_data[,7]*100
barplot(t(as.matrix(nobel_data_rprofile)), hori=TRUE)


nobel_data[1:8,1:6]/nobel_data[9,1:6]

nobel_data_cprofile <- nobel_data[1:8,1:6]/nobel_data[9,1:6]


library(FactoMineR)

res <- CA(nobel_data[,-1])

library(readr)
exudation_data <- read_csv("/Users/glmarschmann/.julia/dev/JuliaDEB/files/exudation_rate.csv")
exudation_frame <- data.frame(exudation_data)
cown <- c('sugars', 'organics', 'aminos', 'fattys', 'nucleos', 'auxins')
colnames(exudation_frame) <- cown

res <- CA(exudation_frame)
