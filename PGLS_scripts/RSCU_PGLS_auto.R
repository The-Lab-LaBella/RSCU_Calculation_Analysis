
#Install Packages
if (!require("ggplot2")) install.packages("ggplot2"); library(ggplot2)

if (!require("reshape2")) install.packages("reshape2"); library(reshape2)

theme_set(theme_classic())

if (!require("dplyr")) install.packages("dplyr"); library(dplyr)

if (!require("readxl")) install.packages("readxl"); library(readxl)

if (!require("readr")) install.packages("readr"); library(readr)

if (!require("caper")) install.packages("caper"); library(caper)

if (!require("nlme")) install.packages("nlme"); library(nlme)

if (!require("purrr")) install.packages("purrr"); library(purrr)

#Inputting Tree and csv files

#joining the tips to the assembly_IDs
yeast_key<-read_excel("1154yeasts_21outrgoups_info_20220408.xlsx")
yeast_key<-yeast_key[,c("assembly_fullID_updated","Tip_ID_on_tree","clade")]

all_tree<-read.nexus("1175taxa_1403OGs_iqtree_ML_rooted.nex")
fixed_env<-read.delim("updated_env_plant_vict_animal_env_abbe_10_5.txt")
yeast_data <- read.csv("Yeast_Matrix.csv")


#Adding the TipIDs to yeast dataframe
tree.tip <- yeast_key$Tip_ID_on_tree
yeast_data['Tip_ID_on_tree'] <- c(yeast_key$Tip_ID_on_tree)
yeast_data <- as.data.frame(yeast_data)



#Empty values for the variables
summ <- c() #empty list for linear model summary
pvalue <- c() #empty list for p values
coeff <- c() #empty list for slope
lam <- c() #empty list for lambda values
codvs <- c() #empty list for codon labels (AAA vs AAC)
rSqr <- c() #empty list for r squared
adSqr <- c() #empty list for adjusted r squared
intercpt <- c() #empty list for intercept

#seed for reproducibility
set.seed(123)


#create a safe version of pgls
pgls.possible=possibly(.f=pgls, otherwise=NULL)


index<-3

#When doing this for real change to 67
while(index<=66){
  
  column<-3
  #When doing this for real change to 67 
  while (column<=66){
    if(column != index){
      #Y=AAA, X=AAT
      
      summ <- pgls.possible( yeast_data[,c(index)] ~ yeast_data[,c(column)], comparative.data(all_tree,yeast_data,"Tip_ID_on_tree"), lambda = "ML")
      codvs <- append(codvs, print(paste0(names(yeast_data)[index]," vs ", names(yeast_data)[column])))
      if(is.null(summ)==T){
        coeff <- append(coeff, "NULL") 
        
        pvalue <- append(pvalue, "NULL")
        
        lam <- append(lam, "NULL")

        rSqr <- append(rSqr, "NULL")

        adSqr <- append(adSqr, "NULL")

        intercpt <- append(intercpt, "NULL")
      }else{
          coeff <- append(coeff, signif(summary(summ)$coefficients[2,1], 5)) 
    
          pvalue <- append(pvalue, signif(summary(summ)$coefficients[2,4], 5))
      
          lam <- append(lam, summ$param[2])

          rSqr <- append(rSqr, signif(summary(summ)$r.squared, 5))

          adSqr <- append(adSqr, signif(summary(summ)$adj.r.squared, 5))

          intercpt <- append(intercpt, coef(summ)[["(Intercept)"]])

      }
      
      print("column is")
      print(column)
    }
    
    
    column<-column+1
  }
  
  index<-index+1
  print("index is")
  print(index)
}


PGLS_df <- data.frame(codons=codvs, slope=coeff, intercept=intercpt, p_value=pvalue, lambda=lam, r_squared=rSqr, adusted_r_squared=adSqr)

write.csv(PGLS_df,"RSCU_Untransformed_PGLS_profile.csv", row.names=FALSE)
#write.csv(PGLS_df,"/projects/labella_lab/all_files_y1000_plus/RSCU_PGLS/RSCU_PGLS_profile.csv", row.names=FALSE)
