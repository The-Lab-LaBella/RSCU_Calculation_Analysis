
#Reading files to calculate basic stats of RSCU average of every codon in each genome

RSCU_data <- list.files(path="/projects/labella_lab/all_files_y1000_plus/RSCU_yeast_files", pattern = "RSCU.csv")

  
for (data in RSCU_data)
{
    loaded_data <- read.csv(data)
  	X <- loaded_data[['Sequence']]
  	Y <- loaded_data[['RSCU']]
  
  	mean(Y, na.rm = TRUE)
  
  	results <- aggregate(x=Y, by=list(loaded_data[['Codon']]), FUN=mean, na.rm=TRUE, na.action=NULL)
  	  	
  	file_header <- "Mean_summary_"
  	
  	filename <- paste(file_header,data, sep ="")
  	
  	file.create(filename)
  	
  	write.table(results, file = filename, append = TRUE)
}
