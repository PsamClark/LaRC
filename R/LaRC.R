#version: 0.9.5.1
library(stringdist)
#version: 1.2.1
library(tidyverse)
#version: 1.2-2
library(splus2R)
#version: 3.7-4
library(plotrix)

"Landscape Roughness Calculator: Function for calculating roughness 
by sampling across a fitness landscape based on sequence data"

LaRC<-function(data_set,dims=2,nsize=16,runs=100,meth="lcs"){
  
  #create distance matrix based on sequences
  dist<-stringdistmatrix(data_set$sequence,method = meth)
  
  #scale distance matrix onto provided dimensions
  cmds<-cmdscale(dist,k=dims)
  
  #combine scaled distances and fit data
  mds<-as.tibble(cmds)
  
  new_set<-
    bind_cols(data_set,mds)
  
  rs_set<-vector('numeric')
  
  #iterate sampling runs
  for (i in 1:runs){
    
    #extract sample set
    samp_set<-
      new_set%>%
      select(-sequence)%>%
      sample_n(nsize)
  #fit  linear model to landscape  
  mod<-lm(data = samp_set,formula = fit~.)
  
  #calculate mean slope of model fit
  ms<-mean(abs(mod$coefficients[2:length(mod$coefficients)]),na.rm = TRUE)
  
  #calculate roughness (square of residuals) over mean slope 
  rs<-sqrt(sum((mod$residuals)^2)/nsize)/ms
  
  rs_set<-c(rs_set,rs)
  }
  
  #calculate mean and standard error of roughness over slope
  mean_rs<-mean(rs_set)
  std_rs<-std.error(rs_set)
  
  #return values
  return(list(mean=mean_rs,se=std_rs))
  
}