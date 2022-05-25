#############################################
####### Identifying alert concentrations ###
####### by a parametric bootstrap approach ###
####### Author: K.Moellenhoff ###############
##############################################

library(matrixStats)
library(DoseFinding)
library(parallel)
library(foreach)
library(doParallel)

### Please note: You have to specify the data set by loading the file in line 20 
### The model has to be specified (up to now you can choose sigEmax or betaMod)
### You can adapt the code easily to other scenarios by changing the doses, the number of samples and the threshold 
### Further you can change the bootstrap repetitions and the significance level
### In summary: You can modify a lot of things between line 20 and 44. Further, in line 91 you can activate/deactivate the parallelization.

# load the dataset of interest
load("DataSet.ScenarioIV.VarLarge.RData")

doses <- c(0, 25, 150, 350, 450, 550, 800, 1000) #concentrations
samples <- 3 #number of samples per concentration
lambda <- log2(1.5) #threshold
scal <- max(doses)*1.2
alpha <- 0.05 #significance level of the test
B1 <- 500 # "outer" bootstrap for critical value (=quantile)
B2 <- 25 # "inner" bootstrap for standard error
# here we construct a grid for searching the LEC (can be modified if required)
epsilon <- 1
grid <- seq(1,max(doses),epsilon)
Nsim <- 1000

# Specify the model: choose just one of these two options: sigEmax or betaMod 

# Option 1
model <- sigEmax
model_name <- "sigEmax"

# Option 2
model <- function(dose, e0, eMax, delta1, delta2) {
  betaMod(dose, e0, eMax, delta1, delta2,scal=scal)
}
model_name <- "betaMod"

### from here on there are no changes/modifications necessary except for turning on/off the paralellization 

# this is needed for parallelizing 
numCores <- detectCores()

# vectors for storing the results
lo.cf <- vector()
crit_val <- vector()
lec <- vector()
result <- vector()
conf_all <- matrix(NA,nrow=Nsim,ncol=length(grid))

# data simulating function
simul <- function(theta,sd)
{
  b <- theta[1]
  c <- theta[2]
  d <- theta[3]
  e <- theta[4]
  resp <- model(rep(doses,each=3),b, c, d, e)
  resp + rnorm(length(resp), mean = 0, sd = sd)
}

# bootstrap function returns replicates of abs(\hat Delta)
bootstrap <- function(theta,B){
  # store bootstrap results
  boot <- matrix(NA,nrow=B,ncol=length(grid))
  theta_boot_vec<- matrix(NA,nrow=B,ncol=length(theta))
  for(m in 1:B){
    data_boot <- simul(theta,sd=sd_est)
    try(mod_boot <- suppressMessages(fitMod(dose=rep(doses, each=3), resp=data_boot, model=model_name)))
    theta_boot <- unname(mod_boot$coef)
    diff_boot <- model(grid,theta_boot[1],theta_boot[2],theta_boot[3],theta_boot[4])-model(0,theta_boot[1],theta_boot[2],theta_boot[3],theta_boot[4])
    boot[m,] <- abs(diff_boot)
    theta_boot_vec[m,] <- theta_boot  
    }
  return(list(boot=boot,theta_boot=theta_boot_vec))
}

# this function returns the simulated values D^* 
# NOTE: if you change do to dopar (line 91), the parallelization is activated. This only works for Linux and MacOS!
crit.val <- function(theta,B){
  D <- vector()
  Dmax <- vector()
  registerDoParallel(numCores)
  res <- foreach (m=1:B,.combine=rbind) %do% {
    inner_boot <- bootstrap(bootst$theta_boot[m,],B=B2)$boot # inner bootstrap for estimating the SE*
    for(l in 1:length(grid)){
      fx <- model(grid[l], theta[1], theta[2], theta[3], theta[4])
      f0 <- model(0, theta[1], theta[2], theta[3], theta[4])
    D[l] <- (boot[m,l]-abs(fx-f0))/sd(inner_boot[,l])
    }
D
  }
  stopImplicitCluster()
  Dmax <- rowMaxs(res)
   return(Dmax)
  }


conf.bands <- function(theta){
  for(l in 1:length(grid)){
    fx <- model(dose = grid[l], theta[1], theta[2], theta[3], theta[4])
    f0 <- model(dose = 0, theta[1], theta[2], theta[3], theta[4])
    #pointwise estimated variance (each l is one point of the grid)
    var.FC <- var(boot[,l])
    lo.cf[l] <- abs(fx-f0) - c*sqrt(var.FC)
  }
  return(list(low=lo.cf))
}


#####################################################
####### Here starts the evaluation of the data #######
#####################################################

for(i in 1:Nsim){
  # load the i-th dataset
  data <- mat.expression[i,]
  # fit the model and the standard deviation
  mod_est <- suppressMessages(fitMod(dose=rep(doses, each=3), resp=data, model=model_name))
  sd_est <- sqrt(mod_est$RSS/mod_est$df)
  theta <- unname(coef(mod_est))
  bootst <- bootstrap(theta=theta,B=B1) # outer bootstrap with B1 repetitions of theta and the test statistic 
  boot <- bootst$boot
  crit <- crit.val(theta=theta,B=B1) # here hat D^* is simulated by performing another bootstrap (nested)
  if(all(is.na(crit)) | all(crit==-Inf)){result[i]=NA;next} # in case of numeric errors just skip this dataset 
  c <- quantile(crit,1-alpha,na.rm=TRUE) # take the quantile for the critical value
  crit_val[i] <- c # save it
  conf <- conf.bands(theta=theta) # calculate the confidence band
  conf_all[i,] <- conf$low
  #has an alert been detected?
  if(any(conf$low>=lambda)){result[i]=1;lec[i]=grid[min(which(conf$low>=lambda))]}else{result[i]=0;lec[i]=NA} 
  print(i)
}

# saving the results
saveRDS(list(n=sum(result,na.rm=TRUE),sum.nas=sum(is.na(result)),lec=lec,conf.bands=conf_all),file="Results.rds")

