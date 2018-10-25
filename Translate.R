##Trevor Penix
##Translate

source("https://bioconductor.org/biocLite.R")

for(package in c("Biostrings","annotate")) {
  if(!require(package, character.only=T, quietly=T)){
    biocLite(package)
    library(package, character.only=T, quietyly=TRUE)
  }
  
}

translateDNA <- function(DNA){
  
  #Taking in sequence and translating it.
  sequence <- DNAString(DNA)
  protein <- translate(sequence, genetic.code = GENETIC_CODE)
  
  #Calculating ratios of amino acids
  aas <- letterFrequency(protein, "*")
  aas <- aas/length(protein)
  
  aminoAcids <- c("A","I","L","M","F","W","Y","V","R","H","K","D","E","S","T","N","Q","C","U","G","P")
  for(letter in aminoAcids){
    aas <- c(aas,letterFrequency(protein, letter)/length(protein) )
  }
  
  #Print out of information
  print(paste("Protein Length:", length(protein)))
  print("Ratios of Amino Acids")
  print(aas)
  print("Protein Sequence")
  print(protein)
  
  return(protein)
}


translateRNA <- function(RNA){
  
  #Taking in sequence and translating it.
  sequence <- RNAString(RNA)
  protein <- translate(sequence, genetic.code = GENETIC_CODE)
  
  #Calculating ratios of amino acids
  aas <- letterFrequency(protein, "*")
  aas <- aas/length(protein)
  
  aminoAcids <- c("A","I","L","M","F","W","Y","V","R","H","K","D","E","S","T","N","Q","C","U","G","P")
  for(letter in aminoAcids){
    aas <- c(aas,letterFrequency(protein, letter)/length(protein) )
  }
  
  #Print out of information
  print(paste("Protein Length:", length(protein)))
  print("Ratios of Amino Acids")
  print(aas)
  print("Protein Sequence")
  print(protein)
  
  return(protein)
}


