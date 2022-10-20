RSCU_data <- list.files(path="/projects/labella_lab/all_files_y1000_plus/RSCU_yeast_files", pattern = "RSCU.csv")


for (file in RSCU_data)
{
  data <- read.csv(file)
  data <- t(data)
  file_header <- "Transposed_"
  filename <- paste(file_header,file, sep ="")
  
  file.create(filename)
  
  write.csv(results, file = filename, append = TRUE)
}