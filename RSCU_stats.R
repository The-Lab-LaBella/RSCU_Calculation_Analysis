
#Reading files to calculate basic stats of RSCU average of every codon in each genome

#File Path
RSCU_data <- list.files(path="/projects/labella_lab/all_files_y1000_plus/RSCU_yeast_files", pattern = "RSCU.csv")

#Loop through directory of files
for (data in RSCU_data)
{
    #read in files
    loaded_data <- read.csv(data)
    #Load the sequence name and RSCU values
  	X <- loaded_data[['Sequence']]
  	Y <- loaded_data[['RSCU']]
    
    #calculate the mean
  	mean(Y, na.rm = TRUE)
  
    #calculate the mean of each codon and its RSCU value
  	results <- aggregate(x=Y, by=list(loaded_data[['Codon']]), FUN=mean, na.rm=TRUE, na.action=NULL)
  	  	
    #File name to differentiate
  	file_header <- "Mean_summary_"
  	
    #Append the file name and append it to the regular file name
  	filename <- paste(file_header,data, sep ="")
  	
    #Create the file
  	file.create(filename)
  	
    #Add contents to the created file (RSCU average of each genome for the species)
  	write.table(results, file = filename, append = TRUE)
}
