
# Packages ----------------------------------------------------------------

rm(list=ls())
packages = c("evgam","evd","ismev","quantreg")
package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)

# Accounting for marginal non-stationarity ------------------------------------------------------------

#Example computing lambda for bivariate gaussian

rho <- seq(0, 1, 0.1)
w_lam <- seq(0, 1, 0.01)
lambda_ex <- matrix(0, nrow = length(rho), ncol = length(w_lam))
for(i in 1:length(rho)){
  for(j in 1:length(w_lam)){
    if(rho[i]^2 < min(w_lam[j] / (1 - w_lam[j]), (1 - w_lam[j]) / w_lam[j])){
      lambda_ex[i, j] <- (1 - 2 * rho[i] * sqrt(w_lam[j] * (1 - w_lam[j]))) / (1 - rho[i]^2)
    }
    else{
      lambda_ex[i, j] <- max(w_lam[j], 1 - w_lam[j])
    }
  }
}
colfunc = colorRampPalette(c("red","blue"))
cols = colfunc(length(rho))

pdf('lambda_gauss.pdf')
plot(w,lambda_ex[1,],type="l",lwd=2,col=cols[1],xlab="w",ylab='lambda(w)',
     xlim=c(0,1),ylim=c(0.4,1.2),cex.main=1.5,cex.axis=1.5,cex.lab=1.5,
     main = 'Angular Dependence of Bivariate Gaussian')
for(i in 2:length(rho)){
  lines(w,lambda_ex[i,],lwd=2,col=cols[i])
}                                  
dev.off()


##Dat cleaning and testing various BC datasets

#load(file="UKCP18_data.RData")#file containing marginal datasets, with time covariates. Both marginal dataset can be obtained from https://ukclimateprojections-ui.metoffice.gov.uk/products

setwd("/Users/JZ/Desktop/Qualifying-Paper-2/data/EC_raw 2")

files <- c('sechelt.csv', 'callaghanvalley.csv', 'summerland.csv', 
           'pointatkinson.csv', 'powellriver.csv', 'fannyisland.csv', 'deaselake.csv')

#for(file in files){
file <- 'summerland.csv'
  data_raw <- read.csv(paste('/Users/JZ/Desktop/Qualifying-Paper-2/data/EC_raw 2/', file, sep = ''), skip = 1)
  data_raw <- data_raw[,c('air_temperature', 'air_temperature_yesterday_high', 'relative_humidity', 'time', 'air_temperature_yesterday_low')]
  colnames(data_raw) <- c('temp', 'temp_max', 'rel_humid', 'time', 'temp_min')
  data_raw$time <- trimws(data_raw$time)
  data_raw$time <- trimws(data_raw$time)
  data_raw$rel_humid <- trimws(data_raw$rel_humid)
  data_raw <- data_raw[data_raw[,'rel_humid'] != 'None', ]
  
  data_raw$Year <- as.numeric(substr(data_raw[,'time'], 1, 4))
  data_raw$Month <- as.numeric(substr(data_raw[,'time'], 6, 7))
  data_raw$Day <- as.numeric(substr(data_raw[,'time'], 9, 10))
  data_raw$Dryness <- 100 - as.numeric(data_raw$rel_humid)
  data_raw$Temperature <- as.numeric(data_raw$temp_max)
  data_raw <- data_raw[(data_raw[,'temp_max'] != ' None'), ]
  
  data_raw <- data_raw[(data_raw[,'Month'] %in% c(6, 7, 8)), ]
  data_raw$DayIndex <- data_raw$Day + 30*(data_raw$Month == 7) + 61*(data_raw$Month == 8)
  data_raw$TimeIndex <- data_raw$Day + 30*(data_raw$Month == 7) + 61*(data_raw$Month == 8) + 92*(data_raw$Year - 2012)
  
  dryness_data <- data_raw[,c('TimeIndex', 'DayIndex', 'Month', 'Year', 'Dryness')]
  temperature_data <- data_raw[,c('TimeIndex', 'DayIndex', 'Month', 'Year', 'Temperature')]
  
  
  colfunc = colorRampPalette(c("red","blue"))
  cols = colfunc(dim(temperature_data)[1])
  
  #printing diagnostics
  
  pdf('climate_diagnostic.pdf', width = 10, height = 5)
  par(mfrow = c(1, 2))
  plot(temperature_data$Temperature, dryness_data$Dryness, col = cols, 
       xlab = 'Temperature', ylab = 'Dryness', main = 'Summerland Climate', 
       cex = 1, pch = 16)
  
  plot(seq(1, dim(temperature_data)[1], 1), temperature_data$Temperature, type = 'l',
       xlab = 'Time', ylab = 'Temperature', main= 'Summerland Temperature over Time')
  dev.off()
#}




data = cbind(temperature_data$Temperature,dryness_data$Dryness) #the bivariate dataset in question
n = dim(dryness_data)[1] #number of observations

#Fitting GAMs to remove non-stationarity in the body of data
#Temperature Data

m_temp <- gam( formula = list(Temperature~s(DayIndex, bs = "cc", k = 90) + s(TimeIndex),~s(TimeIndex)),family = gaulss(b=0),data = temperature_data) #this gam fits long term trends to both location and scale, along with a seasonal trend to the location
summary(m_temp) #Each function is significant
mu_temp = m_temp$fitted.values[,1] #obtaining fitted location function
sig_temp = 1/m_temp$fitted.values[,2] #obtaining fitted scale function

sdata = matrix(0,ncol=2,nrow=n) #matrix for pre-processed data
sdata[,1] = (temperature_data$Temperature-mu_temp)/(sig_temp) #temperature data after pre-processing

#Dryness Data

m_dry <- gam( formula = list(Dryness~s(DayIndex, bs = "cc", k = 90) + s(TimeIndex),~s(TimeIndex)),family = gaulss(b=0),data = dryness_data) #this gam fits long term trends to both location and scale, along with a seasonal trend to the location
summary(m_dry) #Each function is significant
mu_dry = m_dry$fitted.values[,1] #obtaining fitted location function
sig_dry = 1/m_dry$fitted.values[,2] #obtaining fitted scale function

sdata[,2] = (dryness_data$Dryness-mu_dry)/(sig_dry) #dryness data after pre-processing

#fitting GPDs to the upper tail, likelihood ratio tests are used to test signficance of fitted trends 

nyear = 90 #number of points in each year

lin_trend = as.matrix(c(1:n)/(n+1)) #for evaluating linear trends in scale parameter
cyc_trend = as.matrix(cbind(cos(2*pi*c(1:n)/nyear),sin(2*pi*c(1:n)/nyear))) #for evaluating seasonal trends in scale parameter

thresh1 = quantile(sdata[,1],0.96) #high threshold of pre-processed temperature data
fit1_stat = gpd.fit(xdat = sdata[,1],threshold = thresh1,show=F) #stationary GPD fit
fit1_linear = gpd.fit(xdat = sdata[,1],threshold = thresh1,ydat=lin_trend,siglink = exp,sigl = 1,show=F) #GPD fit with linear trend in scale
dev = 2*(fit1_stat$nllh-fit1_linear$nllh);print(dev) #non significant deviance value at 5% level
fit1_seasonal = gpd.fit(xdat = sdata[,1],threshold = thresh1,ydat=cyc_trend,siglink = exp,sigl = 1:2,show=F) #GPD fit with linear trend in scale
dev = 2*(fit1_stat$nllh-fit1_seasonal$nllh);print(dev) #significant deviance value at 5% level
#use non stationary GPD fit with seasonal trends

thresh2 = quantile(sdata[,2],0.96)
fit2_stat = gpd.fit(xdat = sdata[,2],threshold = thresh2,show=F) #stationary GPD fit
fit2_linear = gpd.fit(xdat = sdata[,2],threshold = thresh2,ydat=lin_trend,siglink = exp,sigl = 1,show=F) #GPD fit with linear trend in scale
dev = 2*(fit2_stat$nllh-fit2_linear$nllh);print(dev) #non significant deviance value for linear trend at 5% level
fit2_seasonal = gpd.fit(xdat = sdata[,2],threshold = thresh2,ydat=cyc_trend,siglink = exp,sigl = 1:2,show=F) #GPD fit with linear trend in scale
dev = 2*(fit2_stat$nllh-fit2_seasonal$nllh);print(dev) non #significant deviance value for seasonal trend at 5% level
#use stationary GPD

unif_data = matrix(0,ncol=2,nrow=n) #matrix for storing the data transformed to uniform margins using the Probability integral transform 

#transforming the data to uniform margins using fitted GPDs

q1_gpd = 0.96
unif_data[sdata[,1]<=thresh1,1] = (rank(sdata[,1])[sdata[,1]<=thresh1])/(n+1)
unif_data[sdata[,1]>thresh1,1] = 1-(1-q1_gpd)*pgpd((sdata[,1])[sdata[,1]>thresh1],loc=thresh1,scale=fit1_seasonal$vals[,1],shape=fit1_seasonal$mle[4],lower.tail = FALSE)

q2_gpd = 0.96
unif_data[sdata[,2]<=thresh2,2] = (rank(sdata[,2])[sdata[,2]<=thresh2])/(n+1)
unif_data[sdata[,2]>thresh2,2] = 1-(1-q2_gpd)*pgpd((sdata[,2])[sdata[,2]>thresh2],loc=thresh2,scale=fit2_stat$mle[1],shape=fit2_stat$mle[2],lower.tail = FALSE)

data_exp = apply(unif_data, 2, FUN = qexp) #matrix of data on exponential margins


# Supplementary plots  ----------------------------------------------------

#Fitted mu and sigma functions against empirical estimates

#To compute standard error for yearly averages, we simulate from the fitted GAM multivariate normal distributions for parameters and evaluate the resulting posterior distribution

newd = cbind(temperature_data$DayIndex,temperature_data$TimeIndex) #data for prediction matrices
newd = as.data.frame(newd)
names(newd) = c("DayIndex","TimeIndex")

rmvn <- function(n,mu,sig) { ## MVN random deviates
  L <- mroot(sig);m <- ncol(L);
  t(mu + L%*%matrix(rnorm(m*n),m,n)) 
}

mu_index_temp = 1:(which(names(coef(m_temp))=="(Intercept).1")-1) #which smooth terms relate to location
sig_index_temp = (which(names(coef(m_temp))=="(Intercept).1")):length(coef(m_temp))#which smooth terms relate to scale

mu_index_dry = 1:(which(names(coef(m_dry))=="(Intercept).1")-1) #which smooth terms relate to location
sig_index_dry = (which(names(coef(m_dry))=="(Intercept).1")):length(coef(m_dry))#which smooth terms relate to scale

Xp_temp <- predict(m_temp,newd,type="lpmatrix") #matrix of predictions
Xp_dry <- predict(m_dry,newd,type="lpmatrix") #matrix of predictions

#empty vectors for storing

mu_x = c()
mu_y = c()
sig_x = c()
sig_y = c()
est_mu_x = c()
est_mu_x_upper = c()
est_mu_x_lower = c()
est_mu_y = c()
est_mu_y_upper = c()
est_mu_y_lower = c()
est_sig_x = c()
est_sig_x_upper = c()
est_sig_x_lower = c()
est_sig_y = c()
est_sig_y_upper = c()
est_sig_y_lower = c()
count = 1
years = unique(temperature_data$Year) #vector of years in the observation period

for(i in years){
  window_days = temperature_data$Year == i #subsetting the data for each year
  
  data_ind = data[window_days,]
  
  #empirical estimates of mean and sd
  
  mu_x[count] = mean(data_ind[,1])
  mu_y[count] = mean(data_ind[,2])
  
  sig_x[count] = sd(data_ind[,1])
  sig_y[count] = sd(data_ind[,2])
  
  #averaged GAM estimates of location and scale
  
  est_mu_x[count] = mean(mu_temp[window_days])
  est_sig_x[count] = mean(sig_temp[window_days])
  
  br <- rmvn(1000,coef(m_temp),m_temp$Vp) ## 1000 replicate param vectors
  
  res_mu <- rep(0,1000)
  res_sig <- rep(0,1000)
  for (j in 1:1000){ 
    pr_mu = Xp_temp[window_days,mu_index_temp] %*% br[j,mu_index_temp] ## replicate yearly predictions of location
    
    pr_sig = Xp_temp[window_days,sig_index_temp] %*% br[j,sig_index_temp] ## replicate yearly predictions of scale
    
    res_mu[j] <- mean(pr_mu) 
    
    res_sig[j] <- mean(exp(pr_sig))
  }
  
  est_mu_x_lower[count] = quantile(res_mu,0.025) #quantiles from resulting posterior distributions. `x' here is temperature
  est_mu_x_upper[count] = quantile(res_mu,0.975)
  
  est_sig_x_lower[count] = quantile(res_sig,0.025)
  est_sig_x_upper[count] = quantile(res_sig,0.975)
  
  est_mu_y[count] = mean(mu_dry[window_days])
  est_sig_y[count] = mean(sig_dry[window_days])
  
  br <- rmvn(1000,coef(m_dry),m_dry$Vp) ## 1000 replicate param. vectors
  
  res_mu <- rep(0,1000)
  res_sig <- rep(0,1000)
  for (j in 1:1000){ 
    pr_mu = Xp_dry[window_days,mu_index_dry] %*% br[j,mu_index_dry] ## replicate yearly predictions
    
    pr_sig = Xp_dry[window_days,sig_index_dry] %*% br[j,sig_index_dry] ## replicate yearly predictions
    
    res_mu[j] <- mean(pr_mu) 
    
    res_sig[j] <- mean(exp(pr_sig))
  }
  
  est_mu_y_lower[count] = quantile(res_mu,0.025)
  est_mu_y_upper[count] = quantile(res_mu,0.975)
  
  est_sig_y_lower[count] = quantile(res_sig,0.025)
  est_sig_y_upper[count] = quantile(res_sig,0.975)
  
  count = count + 1
}

par(mfrow=c(1,2),mgp=c(2.5,1,0),mar=c(5,4,4,2)+0.1)
plot(years,mu_x,type="l",lwd=3,xlab=paste("Year"),ylab="",ylim=c(range(mu_x,mu_y,mu_temp,mu_dry,est_mu_x_lower,est_mu_x_upper,est_mu_y_lower,est_mu_y_upper)[1],range(mu_x,mu_y,mu_temp,mu_dry,est_mu_x_lower,est_mu_x_upper,est_mu_y_lower,est_mu_y_upper)[2]+10),cex.lab=1.5, cex.axis=1.5)
polygon(x = c(years, rev(years)),y = c(est_mu_x_lower, rev(est_mu_x_upper)),col =  adjustcolor("red", alpha.f = 0.2), border = NA)
lines(years,est_mu_x,col=2,lwd=3)
lines(years,mu_y,col=1,lty=2,lwd=3)
polygon(x = c(years, rev(years)),y = c(est_mu_y_lower, rev(est_mu_y_upper)),col =  adjustcolor("green", alpha.f = 0.2), border = NA)
lines(years,est_mu_y,lty=2,col=3,lwd=3)
legend(2020, 90,bg="white",legend=c("Empirical mean temperature","Estimated location temperature","Empirical mean dryness","Estimated location dryness"),col=c(1:2,1,3),lty=c(1,1,2,2),lwd=3,cex=1.2)

plot(years,sig_x,type="l",lwd=3,xlab=paste("Year"),ylab="",ylim=c(range(sig_x,sig_y,sig_temp,sig_dry,est_sig_x_upper,est_sig_x_lower,est_sig_y_lower,est_sig_y_upper)[1],range(sig_x,sig_y,sig_temp,sig_dry,est_sig_x_upper,est_sig_x_lower,est_sig_y_lower,est_sig_y_upper)[2]+10),cex.lab=1.5, cex.axis=1.5)
polygon(x = c(years, rev(years)),y = c(est_sig_x_lower, rev(est_sig_x_upper)),col =  adjustcolor("red", alpha.f = 0.2), border = NA)
lines(years,est_sig_x,col=2,lwd=3)
lines(years,sig_y,col=1,lty=2,lwd=3)
polygon(x = c(years, rev(years)),y = c(est_sig_y_lower, rev(est_sig_y_upper)),col =  adjustcolor("green", alpha.f = 0.2), border = NA)
lines(years,est_sig_y,lty=2,col=3,lwd=3)
legend(c(2020, 2024), c(30, 22),bg="white",legend=c("Empirical s.d. \n temperature","Estimated scale \n temperature","Empirical s.d. \n dryness","Estimated scale \n dryness"),col=c(1:2,1,3),lty=c(1,1,2,2),lwd=3,cex=0.5)


window_length_gap = 2
years2 = years[-c(1:window_length_gap,(length(years)-window_length_gap+1):(length(years)))] #years for constructing +/- 15 year intervals

rate_x = c()
upper_x = c()
lower_x = c()
rate_y = c()
upper_y = c()
lower_y = c()
count = 1
for(i in years2){
  window_days = temperature_data$Year >= i-window_length_gap & temperature_data$Year <= i+window_length_gap #subset of indices +/- 15 year interval for year i 
  dataexp = data_exp[window_days,] #subset of data for window
  
  #MLE estimates of exponential rate parameters, with asymptotic 95% confidence intervals
  
  N = dim(dataexp)[1]
  rate_x[count] = 1/mean(dataexp[,1])
  upper_x[count] = rate_x[count] + 1.96*sqrt((rate_x[count]^2)/N)
  lower_x[count] = rate_x[count] - 1.96*sqrt((rate_x[count]^2)/N)
  rate_y[count] = 1/mean(dataexp[,2])
  upper_y[count] = rate_y[count] + 1.96*sqrt((rate_y[count]^2)/N)
  lower_y[count] = rate_y[count] - 1.96*sqrt((rate_y[count]^2)/N)
  
  count = count + 1
}

par(mfrow=c(1,2),mgp=c(2.5,1,0),mar=c(5,4,4,2)+0.1)
plot(years2,rate_x,type="l",lwd=2,xlab=paste("Reference Year (+/-",window_length_gap,"year windows)"),ylab="",main=expression(paste("Temperature")),ylim=range(rate_x,rate_y,upper_x,upper_y,lower_x,lower_y),cex.lab=1.5, cex.axis=1.5,cex.main=2)
lines(years2,upper_x,lty=2,col=4,lwd=2)
lines(years2,lower_x,lty=2,col=4,lwd=2)
abline(h=1,lwd=3,col=2)
plot(years2,rate_y,type="l",lwd=2,xlab=paste("Reference Year (+/-",window_length_gap,"year windows)"),ylab="",main=expression(paste("Dryness")),ylim=range(rate_x,rate_y,upper_x,upper_y,lower_x,lower_y),cex.lab=1.5, cex.axis=1.5,cex.main=2)
lines(years2,upper_y,lty=2,col=4,lwd=2)
lines(years2,lower_y,lty=2,col=4,lwd=2)
abline(h=1,lwd=3,col=2)

# Fitting BP estimator with constrained coefficients  ---------------------------------------------

bp2d <- function(w, k)
{
  q <- length(w) # number of points of the grid
  basis <- array(0, dim=c(q, k+1))
  
  for(j in 0:k){
    basis[,j+1] <- choose(k,j) * w^j * (1-w)^(k-j)} #BBP basis function
  
  return(basis)
}

para_estimation = function(para_vec,av_emp_vec,n,k,w){
  lam_vec = c()
  beta_matrix = matrix(0,ncol=k+1,nrow = n)
  beta_matrix[,1] = 1
  beta_matrix[,k+1] = 1
  
  mat = matrix(para_vec,ncol=2,byrow=T) #theta
  mat2 = rbind(rep(1,n),1:n) #z_t = (1, t)
  new_mat = t(mat%*%mat2)
  beta_matrix[,2:k] = exp(new_mat) #log link
  
  lam_vec = as.vector(t(bp2d(w=w,k=k)%*%t(beta_matrix)))
  
  minsum = (1/(n*length(w)))*sum(abs(lam_vec + av_emp_vec)) #S (objective)
  return(minsum)
}

minimisation_function = function(w,k,n,av_emp_vec){
  para_init  = rep(c(log(0.5),0),k-1)
  tol = 0.0001
  results = tryCatch(optim(par=para_init,fn=para_estimation,av_emp_vec=av_emp_vec,n=n,k=k,w=w,method="BFGS",control = list(maxit=100000)),error = function(e){1})
  if(is.list(results)){ 
    optim_output = results
  } else {
    optim_output = optim(par=para_init,fn=para_estimation,av_emp_vec=av_emp_vec,n=n,k=k,w=w,method="Nelder-Mead",control = list(maxit=100000))
  }
  results = tryCatch(optim(par=optim_output$par,fn=para_estimation,av_emp_vec=av_emp_vec,n=n,k=k,w=w,method="Nelder-Mead",control = list(maxit=100000)),error = function(e){1})
  if(is.list(results)){ 
    optim_output2 = results
  } else {
    optim_output2 = optim(par=optim_output$par,fn=para_estimation,av_emp_vec=av_emp_vec,n=n,k=k,w=w,method="BFGS",control = list(maxit=100000))
  }
  while(abs(optim_output2$val - optim_output$val)>=tol){
    optim_output = optim_output2
    results = tryCatch(optim(par=optim_output$par,fn=para_estimation,av_emp_vec=av_emp_vec,n=n,k=k,w=w,method="Nelder-Mead",control = list(maxit=100000)),error = function(e){1})
    if(is.list(results)){ 
      optim_output2 = results
    } else {
      optim_output2 = optim(par=optim_output$par,fn=para_estimation,av_emp_vec=av_emp_vec,n=n,k=k,w=w,method="BFGS",control = list(maxit=100000))
    } 
  }
  return(optim_output2$par)
}

w = seq(0,1,length.out = 101) #rays to consider in [0,1]
k = 5 #degree for Bernstein polynomial
time = 1:n #time covariate 

gap = 0.04
q1_init = 0.9 #Base quantile level for fitting non-stationary GPD model 
q1 = seq(0.90,0.95,length.out=30) #q1 quantile levels
q2 = q1 + gap #q2 quantile levels

u_vals = list() #this list is used for saving the q1 quantiles from the min projection. Will be required for computing return curve estimates

eps_vals = list() #This list is used for storing the differences between min projection quantiles (v - u)

for(j in 1:length(q1)){
  eps_vals[[j]] = vector() 
}

min_proj = list()
for(i in 1:length(w)){
  min_proj[[i]] = pmin(data_exp[,1]/w[i],data_exp[,2]/(1-w[i])) #min-projection (K_w,t) at each ray w[i]
}

for(i in 1:length(w)){
  mp_df = as.data.frame(cbind(min_proj[[i]],time,time^2,time^3))
  names(mp_df) = c("MP","t","t2","t3")
  fits_quantreg_q1_init <- rq(MP ~ t + t2 + t3, tau = q1_init, data = mp_df)$fitted.values #Computing base threshold level using standard quantile regression
  
  mp_df_exc = data.frame(MP = mp_df$MP, t=mp_df$t, MP_exc = mp_df$MP - fits_quantreg_q1_init)
  
  mp_df_exc = mp_df_exc[mp_df_exc$MP_exc>0,] #Computing threshold exceedances for min-projection
  
  fmla_gpd = 'MP_exc ~ s(t, k=15, bs="cr")' #Specifying GAM model for the GPD scale parameter
  
  m_gpd = evgam(list(as.formula(fmla_gpd), ~1) , data=mp_df_exc, family="gpd") #Fitting non-stationary GPD to threshold exceedances using evgam.
  
  predicted_q1 = predict(m_gpd,newdata=data.frame(t=time),type="quantile",prob=(q1-q1_init)/(1-q1_init)) #min-projection quantile estimates for q1 (minus the threshold)
  predicted_q2 = predict(m_gpd,newdata=data.frame(t=time),type="quantile",prob=(q2-q1_init)/(1-q1_init)) #min-projection quantile estimates for q2 (minus the threshold)
  u_vals[[i]] = predicted_q1 + fits_quantreg_q1_init #min-projection quantile estimates for q1
  
  rm(m_gpd)
  
  for(j in 1:length(q1)){
    eps_vals[[j]][((i-1)*n+1):(i*n)] = predicted_q2[,j] - predicted_q1[,j] #Computing differences between quantile estimates
  }
}

len_vec = n*length(w)
av_emp_vec = c()

for(r in 1:len_vec){ #This for loop generates a vector of pointwise ADF estimates used for optimising the Bernstein Polynomial coefficients
  temp = c()
  for(j in 1:length(q1)){
    temp[j] = (1/eps_vals[[j]][r])*log((1-q2[j])/(1-q1[j])) #lambda_qr
  }
  av_emp_vec[r] = mean(temp) #hat(lambda_qr)
}

rm(eps_vals) #removed to save memory

para_estimates = minimisation_function(w=w,k=k,n=n,av_emp_vec=av_emp_vec) #Estimating the Bernstein Polynomial parameters via pointwise ADF estimates 

# Estimating non-stationary ADFs and return curves ------------------------

time_indices = which(temperature_data$DayIndex == 45)[floor(seq(1,length(which(temperature_data$DayIndex == 45)),length=11))] #time indices for July 15th at a subset of years over observation period

temperature_data$DayIndex[time_indices] #should all be 45
temperature_data$Month[time_indices] #should all be 7
temperature_data$Year[time_indices] #should be increasing subset of years from 1981-2080

est_lam_bp = list() #list for storing estimated ADFs for time points corresponding to July 15th 

for(i in 1:length(time_indices)){
  est_lam_bp[[i]] = vector()
  beta_t = c()
  beta_t[1] = 1
  beta_t[k+1] = 1
  for(j in 2:k){
    temp = para_estimates[(j-1)*2-1] + para_estimates[(j-1)*2]*((1:n)[time_indices[i]])
    beta_t[j] = exp(temp) #estimate of beta cofficient at time indices 'i', exponential link function
  }
  est_lam_bp[[i]] = bp2d(w=w,k=k)%*%beta_t #get lambda as BBP
  
  #Imposing the lower bound condition. 
  for(j in 1:length(w)){
    if(est_lam_bp[[i]][j] < max(w[j],1-w[j])){
      est_lam_bp[[i]][j] = max(w[j],1-w[j])
    }
  }
  
  #Imposing the shape constraints on lambda
  for(j in 50:1){
    if((w/est_lam_bp[[i]])[j+1]<(w/est_lam_bp[[i]])[j] ){
      est_lam_bp[[i]][j] = (w[j])*(est_lam_bp[[i]]/w)[j+1] 
    }
    if(((1-w)/est_lam_bp[[i]])[j+1]>((1-w)/est_lam_bp[[i]])[j] ){
      est_lam_bp[[i]][j] = ((1-w)[j])*(est_lam_bp[[i]]/(1-w))[j+1] 
    }
  }
  for(j in 52:(length(w))){
    if((w/est_lam_bp[[i]])[j]<(w/est_lam_bp[[i]])[j-1] ){
      est_lam_bp[[i]][j] = (w[j])*(est_lam_bp[[i]]/w)[j-1] 
    }
    
    if(((1-w)/est_lam_bp[[i]])[j]>((1-w)/est_lam_bp[[i]])[j-1] ){
      est_lam_bp[[i]][j] = ((1-w)[j])*(est_lam_bp[[i]]/(1-w))[j-1] 
    }
  }

}

#My GAM estimates for lambda


## attempt to get thresholds over all z_t and w

# all_rays <- rep(0, length(w) * length(time))
# for(i in 1:length(w)){
#   start <- (i - 1) * length(time) + 1
#   end <- i * length(time)
#   all_rays[start:end] <- min_proj[[i]]
# }
# all_rays_df <- data.frame(mp = all_rays, t = rep(time, length(w)), w = rep(w, rep(length(time), length(w))))
# ray_quant <- rq(all_rays ~ w + t + 1, tau = q1[15], data = all_rays_df)$fitted.values

gam_lam <- matrix(0, ncol = length(w), nrow = length(time_indices))
for(i in 1:10){
  # ray_mp <- rep(0, length(w))
  # ray_quant
  # for(j in 1:length(w)){
  #   ray_mp[j] <- min_proj[[j]][time_indices[i]]
  # }
  # ray_quant <- rq(ray_mp ~ w + 1, tau = q1[15], data = data.frame(ray_mp = ray_mp, w = w))$fitted.values
  # 
  # ray_u_df <- data.frame(w = w, ray_mp = ray_mp, ray_exc = ray_mp - ray_quant)
  # ray_u_df <- ray_u_df[ray_u_df$ray_exc > 0, ]
  # print(ray_u_df)
  #evgam_fit <- evgam(list(as.formula('ray_exc ~ s(w, k=4, bs="cr")')) , data=ray_u_df, family = 'exponential')
  #pred <- predict(evgam_fit,newdata=data.frame(w=w),type="link")
  #print(pred)
  
  #print(evgam_fit$lograte)
  
  ray_mp <- min_proj[[i]]
  ray_u <- u_vals[[i]][,15]
  ray_df <- data.frame(t = time, exc = ray_mp - ray_u)
  ray_df <- ray_df[ray_df$exc > 0,]

  evgam_fit <- evgam(list(as.formula('exc ~ s(t, k=15, bs="cr") + 1')) , data=ray_df, family = 'exponential')
  
  pred <- predict(evgam_fit,newdata=data.frame(t=time),type="link") #min-projection quantile estimates for q1 (minus the threshold)
  print(exp(pred[time_indices,1]))
  gam_lam[,i] <- pred[time_indices,1]
}

#Imposing the lower bound condition.
for(i in 1:length(time_indices)){
  for(j in 1:length(w)){
    if(gam_lam[i, j] < max(w[j],1-w[j])){
      gam_lam[i, j] = max(w[j],1-w[j])
    }
  }
  for(j in 50:1){
    if((w/gam_lam[i,])[j+1]<(w/gam_lam[i,])[j] ){
      gam_lam[i,j] = (w[j])*(gam_lam[i,]/w)[j+1]
    }
    if(((1-w)/gam_lam[i,])[j+1]>((1-w)/gam_lam[i,])[j] ){
      gam_lam[i,j] = ((1-w)[j])*(gam_lam[i,]/(1-w))[j+1]
    }
  }
  for(j in 52:(length(w))){
    if((w/gam_lam[i,])[j]<(w/gam_lam[i,])[j-1] ){
      gam_lam[i,j] = (w[j])*(gam_lam[i,]/w)[j-1]
    }
    
    if(((1-w)/gam_lam[i,])[j]>((1-w)/gam_lam[i,])[j-1] ){
      gam_lam[i,j] = ((1-w)[j])*(gam_lam[i,]/(1-w))[j-1]
    }
  }
}


gam_lam

#Estimating return curves on exponential margins 

est_curves_bp = list() #list for storing return curve estimates on exponential margins
prob = 1/(90*10000) #probability considered, corresponds to 10000 year return period in stationary setting

for(i in 1:length(time_indices)){
  est_curves_bp[[i]] = matrix(0,ncol=2,nrow=length(w))
  
  for(j in 1:length(w)){
    x_vec = c()
    y_vec = c()
    for(q in 1:length(q1)){ #see section 3.6 of paper
      r = -(1/est_lam_bp[[i]][j])*log(prob/(1-q1[q]))
      u = (u_vals[[j]])[time_indices[i],q]
      x_vec[q] = w[j]*(u+r)
      y_vec[q] = (1-w[j])*(u+r)
    }
    est_curves_bp[[i]][j,] = c(mean(x_vec),mean(y_vec)) #average curve over quantile levels
  }
}

#Transforming curves to uniform margins

unif_curves = list()

for(i in 1:length(time_indices)){
  unif_curves[[i]] = matrix(0,ncol=2,nrow=length(w))
  unif_curves[[i]] = apply(est_curves_bp[[i]],2,pexp)
}

#Transforming curves to margins of pre-processed data

norm_curves = list()

for(i in 1:length(time_indices)){
  norm_curves[[i]] = matrix(0,ncol=2,nrow=length(w))
  for(j in 1:length(w)){
    if(unif_curves[[i]][j,1]>q1_gpd){
      norm_curves[[i]][j,1] = qgpd((unif_curves[[i]][j,1]-q1_gpd)/(1-q1_gpd),loc=thresh1, scale = exp(t(c(1,cos(2*pi*time_indices[i]/nyear),sin(2*pi*time_indices[i]/nyear)))%*%fit1_seasonal$mle[1:3]), shape = fit1_seasonal$mle[4])
      #norm_curves[[i]][j,1] = qgpd((unif_curves[[i]][j,1]-q1_gpd)/(1-q1_gpd),loc=thresh1, scale = fit1_seasonal$vals[,1], shape = fit1_seasonal$mle[4])
    } else {
      norm_curves[[i]][j,1] = quantile(sdata[,1],unif_curves[[i]][j,1])
    }
    if(unif_curves[[i]][j,2]>q2_gpd){
      norm_curves[[i]][j,2] = qgpd((unif_curves[[i]][j,2]-q2_gpd)/(1-q2_gpd),loc=thresh2, scale = fit2_stat$mle[1], shape = fit2_stat$mle[2])
      #norm_curves[[i]][j,2] = qgpd((unif_curves[[i]][j,2]-q2_gpd)/(1-q2_gpd),loc=thresh2, scale = exp(t(c(1,time_indices[i]/(n+1),cos(2*pi*time_indices[i]/nyear),sin(2*pi*time_indices[i]/nyear)))%*%fit2_nonstat$mle[1:4]), shape = fit2_nonstat$mle[5])
    } else {
      norm_curves[[i]][j,2] = quantile(sdata[,2],unif_curves[[i]][j,2])
    }
  }
  
}

#Transforming curves to original margins of both time series; for this, we perform the reverse transformation of the pre-processing

ns_curves = list()

for(i in 1:length(time_indices)){
  ns_curves[[i]] = matrix(0,ncol=2,nrow=length(w))
  
  ns_curves[[i]][,1] = norm_curves[[i]][,1]*sig_temp[time_indices[i]] + mu_temp[time_indices[i]]
  ns_curves[[i]][,2] = norm_curves[[i]][,2]*sig_dry[time_indices[i]] + mu_dry[time_indices[i]]
  
  #Here, we correct curves using theoretical properties derived in Murphy-Barltrop et al. (2023). 
  
  for(k in 1:(length(ns_curves[[i]][,1])-2)){#This ensures the curves have the correct shape
    if(ns_curves[[i]][k+1,1]<ns_curves[[i]][k,1]){
      ns_curves[[i]][k+1,1] = ns_curves[[i]][k,1]
    }
    if(ns_curves[[i]][k+1,2]>ns_curves[[i]][k,2]){
      ns_curves[[i]][k+1,2] = ns_curves[[i]][k,2]
    }
  }
  
  ns_curves[[i]] = rbind(c(-10,ns_curves[[i]][1,2]),ns_curves[[i]],c(ns_curves[[i]][length(w),1],-10))   #This is used to inform the values of the curves at the margins
}

#comparing model estimates of eta to empirical 

time_indices2 = 1:n

k = 5
eta_bp = vector()

for(i in 1:length(time_indices2)){ #We compute the ADF for each time point, impose the shape constraints, then extract the corresponding estimate of eta
  est_lam_bp2 = vector()
  beta_t = c()
  beta_t[1] = 1
  beta_t[k+1] = 1
  for(j in 2:k){
    temp = para_estimates[(j-1)*2-1] + para_estimates[(j-1)*2]*((1:n)[time_indices2[i]])
    beta_t[j] = exp(temp)
  }
  est_lam_bp2 = bp2d(w=w,k=k)%*%beta_t
  
  if(sum( est_lam_bp2 < pmax(w,1-w) ) > 0){
    min_ind = min(which( est_lam_bp2 < pmax(w,1-w) ))
    max_ind = max(which( est_lam_bp2 < pmax(w,1-w) ))
    if(min_ind < 52){
      w1_ind = max(which( (est_lam_bp2 < pmax(w,1-w))[1:51] ))
      est_lam_bp2[1:w1_ind] = pmax(w,1-w)[1:w1_ind]
    } 
    if(max_ind >= 52){
      w2_ind = min(which( (est_lam_bp2 < pmax(w,1-w))[52:length(w)] ))+51
      est_lam_bp2[w2_ind:length(w)] = pmax(w,1-w)[w2_ind:length(w)]
    }
  }
  
  for(j in 50:1){
    if((w/est_lam_bp2)[j+1]<(w/est_lam_bp2)[j] ){
      est_lam_bp2[j] = (w[j])*(est_lam_bp2/w)[j+1] 
    }
    if(((1-w)/est_lam_bp2)[j+1]>((1-w)/est_lam_bp2)[j] ){
      est_lam_bp2[j] = ((1-w)[j])*(est_lam_bp2/(1-w))[j+1] 
    }
  }
  for(j in 52:(length(w))){
    if((w/est_lam_bp2)[j]<(w/est_lam_bp2)[j-1] ){
      est_lam_bp2[j] = (w[j])*(est_lam_bp2/w)[j-1] 
    }
    
    if(((1-w)/est_lam_bp2)[j]>((1-w)/est_lam_bp2)[j-1] ){
      est_lam_bp2[j] = ((1-w)[j])*(est_lam_bp2/(1-w))[j-1] 
    }
  }
  
  #Extracting estimating of eta via the ADF w = 0.5
  eta_bp[i] = 1/(2*est_lam_bp2[51])
}

years = unique(temperature_data$Year) #vector of years in the observation period

window_length_gap = 2 #how many years to look above/beyond reference
years = years[-c(1:window_length_gap,(length(years)-window_length_gap+1):(length(years)))] #years for constructing +/- 15 year intervals

eta_est = c()
upper_eta = c()
lower_eta = c()
mean_eta =c() # vector for storing mean eta model estimate for each time window
count = 1
for(i in years){
  window_days = temperature_data$Year >= i-window_length_gap & temperature_data$Year <= i+window_length_gap #+/- year time window from year i
  
  dataexp = data_exp[window_days,] #subset of data for time window
  
  min_proj = pmin(dataexp[,1]/.5,dataexp[,2]/.5) #min-projection at w = 0.5
  
  thresh=quantile(min_proj,0.95)
  eta_est[count] = (0.5)*mean(min_proj[min_proj>thresh] - thresh);if(eta_est[count]>1){eta_est[count] = 1} #hill estimator of eta at 0.5. Cannot exceed 1, hence we bound above 
  
  #95% pointwise confidence intervals for eta
  
  N = length(min_proj[min_proj>thresh])
  upper_eta[count] = eta_est[count] + 1.96*sqrt((eta_est[count]^2)/N);if(upper_eta[count]>1){upper_eta[count] = 1}
  lower_eta[count] = eta_est[count] - 1.96*sqrt((eta_est[count]^2)/N)
  
  mean_eta[count] = mean(eta_bp[window_days]) #mean model estimate in time window
  
  count = count + 1
}

par(mfrow=c(1,2),mgp=c(2.3,1,0),mar=c(5,4,4,2)+0.1)

#averaged model eta estimated against empirical

plot(years,eta_est,type="l",lwd=2,col=1,xlab=paste("Reference Year (+/-",window_length_gap,"year windows)"),ylab=expression(paste(eta," Estimates")),main=expression(bold(paste("Rolling window ",eta, " estimates"))),ylim=range(eta_est,upper_eta,lower_eta),cex.main=1.5,cex.axis=1.5,cex.lab=1.5)
lines(years,upper_eta,lty=2,col=4,lwd=2)
lines(years,lower_eta,lty=2,col=4,lwd=2)
lines(years,mean_eta,lwd=3,col=3)
#legend(2010,0.59,legend=c("Empirical","95% Confidence Region","Averaged Model Estimates"),col=c(1,4,3),lwd=c(2,2,3),lty=c(1,2,1),bg="white",cex=1.2)

#non-stationary return curve estimates over time frame

plot(ns_curves[[1]],type="l",lwd=2,col=cols[1],xlab=expression("Heysham Temperature ("*degree*C*")"),ylab=expression("Heysham Dryness"),main = expression(bold(paste("Non-stationary Return Curve Estimates"))),xlim=c(2,50),ylim=c(2,range(ns_curves)[2]),cex.main=1.5,cex.axis=1.5,cex.lab=1.5)
for(i in 2:(length(time_indices))){
  lines(ns_curves[[i]],lwd=2,col=cols[i])
}
legend(10,30,legend = c("1981","2080"),col=c(2,4),lwd=2,cex=1.5)
point_1981 = ns_curves[[1]][which(w==0.5)+1,] #point from 1981 curve corresponding to w = 0.5
points(point_1981[1],point_1981[2],pch=16,col="green",cex=1.5) 

#assessing point from 1981 curve in the year 2080

#removing non-stationarity in body (relative to 2080)
point_1981[1] = (point_1981[1]- mu_temp[time_indices[length(time_indices)]])/sig_temp[time_indices[length(time_indices)]]
point_1981[2] = (point_1981[2]- mu_dry[time_indices[length(time_indices)]])/sig_dry[time_indices[length(time_indices)]]

#transforming to uniform margins (relative to 2080)
print(point_1981[1] > thresh1) #not greater than gpd threshold for temperature
point_1981[2] = sum(sdata[,2]<=point_1981[2])/(n+1) #estimate using empirical below threshold
print(point_1981[2] > thresh2) #greater than gpd threshold for dryness, so estimate using fitted gpd
point_1981[1] = 1 - (1-q2_gpd)*pgpd(point_1981[1],loc=thresh1, scale = exp(t(c(1,cos(2*pi*time_indices[length(time_indices)]/nyear),sin(2*pi*time_indices[length(time_indices)]/nyear)))%*%fit1_seasonal$mle[1:3]), shape = fit1_seasonal$mle[4],lower.tail = F)

#calculating marginal return periods
RPs = 1/(nyear*(1 - point_1981))
RPs

#transforming to exponential margins
point_1981 = qexp(point_1981)

#empirical value of ray at this point
w_1981 = point_1981[1]/sum(point_1981)
#closest ray in values we have obtained point estimates
w_ind = which.min(abs(w-w_1981))

prob = c() 
#computing joint survival probabilities at each quantile level used for obtained lambda estimator
for(q in 1:length(q1)){
  u = (u_vals[[w_ind]])[time_indices[length(time_indices)],q]
  r = point_1981[1]/w_1981
  v = r - u
  prob[q] = (1-q1[q])*exp(-est_lam_bp[[length(time_indices)]][w_ind]*(v))
}
#computing joint return period
RP2 = 1/(mean(prob)*nyear)
RP2



## printing lambda and return curve plots

colfunc = colorRampPalette(c("red","blue"))
cols = colfunc(length(time_indices))


pdf('results_bbp.pdf', width = 10, height = 5)
par(mfrow = c(1, 2))
plot(w,est_lam_bp[[1]],type="l",lwd=2,col=cols[1],
     xlab="w",ylab = 'lambda(w)', main = 'ADF estimates using BBP',
     xlim=c(0,1),ylim=c(0.4500000,1.5),cex.main=1.5,cex.axis=1.5,cex.lab=1.5)
for(i in 2:length(time_indices)){
  lines(w,est_lam_bp[[i]],lwd=2,col=cols[i])
}
lines(w, pmax(w, 1 - w), lwd = 2, col = 'black')



plot(ns_curves[[1]],type="l",lwd=2,col=cols[1],
     xlab='Temperature',ylab='Dryness',main = "Return Curve Estimates",
     xlim=c(2,50),ylim=c(2,range(ns_curves)[2]),cex.main=1.5,cex.axis=1.5,cex.lab=1.5)
for(i in 2:(length(time_indices))){
  lines(ns_curves[[i]],lwd=2,col=cols[i])
}
points(point_1981[1],point_1981[2],pch=16,col="green",cex=1.5) 
dev.off()


##printin gam
colfunc = colorRampPalette(c("red","blue"))
cols = colfunc(length(time_indices))


pdf('results_gam.pdf')
plot(w,gam_lam[1,],type="l",lwd=2,col=cols[1],
     xlab="w",ylab = 'lambda(w)', main = 'ADF estimates using GAM',
     xlim=c(0,1),ylim=c(0.4500000,1.5),cex.main=1.5,cex.axis=1.5,cex.lab=1.5)
for(i in 2:length(time_indices)){
  lines(w,gam_lam[i,],lwd=2,col=cols[i])
}
lines(w, pmax(w, 1 - w), lwd = 2, col = 'black')
dev.off()










