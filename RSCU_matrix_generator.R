
RSCU_data <- list.files(path="PATH/TO/DIRECTORY", pattern = "^Edited")

#Empty vector
file_matrix <- c()

file.create("Yeast_Matrix.csv")

#Goes through RSCU average files and converts them into a combined matrix file for further analysis
for (file in RSCU_data)
{
  loaded_data <- read.csv(file)
  loaded_data <- loaded_data[c("Mean")]
  file_matrix <- c(file_matrix, loaded_data)
  
}

write.csv(file_matrix, file= "Yeast_Matrix.csv", append=FALSE)
