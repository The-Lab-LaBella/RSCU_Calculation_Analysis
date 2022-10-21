#generating data

library(heatmap3)

RSCU_data <- list.files(path="", pattern = "^Edited")

#Empty vector
file_matrix <- c()

file.create("Yeast_Matrix.csv")

#print(file_matrix)


for (file in RSCU_data)
{
  loaded_data <- read.csv(file)
  loaded_data <- loaded_data[c("Mean")]
  file_matrix <- c(file_matrix, loaded_data)
  
}

#View(t(as.matrix(file_matrix)))
#print(file_matrix)

write.csv(file_matrix, file= "Yeast_Matrix.csv", append=FALSE)
