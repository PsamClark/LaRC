#load library
library(LaRC)

#upload data
data<-read_tsv("/Users/mqbppsc3/functions/LaRC/Demo/fitness_data.txt")

#calculate roughness of data
LaRC(data)
