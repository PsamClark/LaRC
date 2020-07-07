#' Landscape Roughness Calculator: 
#' 
#' Function for calculating roughness by sampling 
#' across a fitness landscape based on sequence data. 
#' Function initially calculates string distances between 
#' sequences and and condenses the resultant distance
#' matrix to the number of dimensions requested
#' a linear fit between the fitness value and the dimensions
#' is then calculated and the average slope for this fit is c
#' alculated.The final step is to calculate the sum of the 
#' residuals to the fit squared divided by the mean slope
#' to give the roughness/slope.
#' This calculation is repeated several times with bootstrap 
#' sampling on the full set and the mean and standard error
#' of the roughness/slope values is calculated.
#' 
#' 
#' 
#' @param data_set tibble containing paired sequences and fitness values
#' @param dims number of dimensions requested for landscape
#' @param nsize number of data points extracted in each sampling
#' @param runs number of times sampling is performed
#' @param meth string distance calculation method 
#' 
#' @return Returns a list object containing the mean 
#' and standard error of the roughness/slope values obtained
#' from sampling.  
#' 
#' @export 

LaRC<-function(data_set,dims=2,nsize=16,runs=100,meth="lcs"){
  
  #create distance matrix based on sequences
  dist<-stringdistmatrix(data_set$sequence,method = meth)
  
  #scale distance matrix onto provided dimensions
  cmds<-cmdscale(dist,k=dims)
  
  #combine scaled distances and fit data
  mds<-as_tibble(cmds)
  
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