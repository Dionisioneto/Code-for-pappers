## -----
## Code of Simualtion studies
## Monte Carlo replicates
## A NEW APPROACH FOR JOINT CURE MODELLING USING THE
## DEFECTIVE GOMPERTZ DISTRIBUTION 
## -----

## [Libraries]
rm(list = ls())
library(INLA)
library(INLAjoint)
library(MASS)


## 1. [Function to generate data from defective Gompertz data]
gen.cure.gompertz = function(n,a,b,p) {
  if (length(b) == 1) b <- rep(b, n)
  if (length(p) == 1) p <- rep(p, n)
  rm = rbinom(n=n,size=1,prob=1-p)
  t=rep(NA,n)
  for(i in 1:n){
    t[i]=ifelse(rm[i]==0,
                Inf,
                log((-(a/b[i])*log(1-runif(n=1,min=0,max=1-p[i]))) + 1)*(1/a))
  }
  t_finite = ifelse(t==Inf,0,t)
  u2 = runif(n=n,0,max(t_finite))
  t2 = pmin(t,u2)
  delta = ifelse(t<u2,1,0)
  return(cbind(t2,delta))
}

## 1. [Function to generate longitudinal and survival data from defective Gompertz data]
gen.joint.cure.dgompertz = function(n,alpha,gammas,betas.fixed,
                                    sd.int,sd.slope,rho,
                                    gap, followup){
  
  nmesindiv=followup/gap+1 # number of individual measurement
  cov=sd.int*sd.slope*rho
  Sigma=matrix(c(sd.int^2,cov,
                 cov,sd.slope^2),ncol=2,nrow=2)
  
  MVnorm <- mvrnorm(n = n, mu = rep(0,2), Sigma = Sigma)
  
  b_int = MVnorm[,1]
  b_intY <- rep(b_int, each=nmesindiv) # b0
  b_slo = MVnorm[,2]
  b_sloY <- rep(b_slo, each=nmesindiv) # b1
  
  binXs=rbinom(n,1, 0.8) # binary covariate
  #binX=rbinom(n,1, 0.5) # binary covariate
  
  bin=rep(binXs, each=nmesindiv)
  ctsX=round(runif(n,0, 1),5) # continuous covariate
  cts=rep(ctsX, each=nmesindiv)
  mestime=seq(0,followup,gap) # measurement times
  timej=rep(mestime, n) # time column
  nmesy= nmesindiv*n# number of longi measurements
  
  # linear predictor
  b_0=betas.fixed[1];b_1=betas.fixed[2];b_2=betas.fixed[3];b_3=betas.fixed[4]
  
  linPred <- b_0+b_intY+(b_1+b_sloY)*timej+b_2*cts+b_3*bin
  # count longitudinal outcome
  biomarker <- rpois(nmesy, exp(linPred))
  idY<-as.factor(rep(1:n, each=nmesindiv)) # id
  #longitudinal dataset
  longDat <- data.frame(id=as.integer(idY), timej, biomarker, cts, bin)
  
  # now our covariates are b_int and b_slo (from longi)
  eta = cbind(1, b_int, b_slo,binXs) %*% c(gamma0,gamma1,gamma2,gamma3)
  mu = exp(eta)
  p_par = exp(mu/alpha)
  
  survDat = gen.cure.gompertz(n=n,a=alpha,b=mu,p=p_par)
  colnames(survDat) = c("time", "delta")
  survDat = data.frame(id=1:n,survDat,binXs,b_int,b_slo)
  
  return(list(survdata=survDat,longdata=longDat))
}

## (fisrt scenario of simualtion) 
# beta0=2.5# intercept
# beta1=-0.2 # slope
# beta2=-0.01 # continuous covariate
# beta3=0.1 # binary covariate
# 
# betas.fixed0=c(beta0, beta1, beta2, beta3)
# sd.int0 = 0.25
# sd.slo0 = 0.25
# rho0 = -0.05
# 
# alpha0=-0.66 ##defective
# gamma0=-0.68;gamma1=0.68;gamma2=0.17;gamma3=-0.37 # gamma3:=psi1
# gammas0 = c(gamma0,gamma1,gamma2,gamma3)
# 
# parameters = c(alpha0,gammas0,betas.fixed0,sd.int0,sd.slo0,rho0)
# names(parameters)=c("alpha","gamma0","gamma1", "gamma2","gamma3",
#                     "beta0","beta1","beta2","beta3",
#                     "sd.int","sd.slo","rho")
# 

## (second scenario of simualtion)

## biomarkers
beta0=1# intercept
beta1=-1 # slope
beta2=-0.1 # continuous covariate
beta3=0.5 # binary covariate

betas.fixed0=c(beta0, beta1, beta2, beta3)
sd.int0 = 0.3
sd.slo0 = 0.3
rho0 = 0.1

alpha0=-1 ##defective
gamma0=0.8;gamma1=2;gamma2=1.2;gamma3=0.5
gammas0 = c(gamma0,gamma1,gamma2,gamma3)

parameters = c(alpha0,gammas0,betas.fixed0,sd.int0,sd.slo0,rho0)
names(parameters)=c("alpha","gamma0","gamma1", "gamma2","gamma3",
                    "beta0","beta1","beta2","beta3",
                    "sd.int","sd.slo","rho")

nsamples = c(100,500,1000)
niter = 1000

## the order for estimates is always (alpha, gamma0, gamma1, gamma2, beta0, beta1, beta2, beta3, sd.int0, sd.slope0, rho)
res.matrix = array(data=NA, dim = c(niter, length(parameters)*4, length(nsamples)),
                   dimnames=list(1:niter, c(rep("mean",length(parameters)), rep("sd",length(parameters)),
                                            rep("bias",length(parameters)), rep("cp",length(parameters))),nsamples))

censoring_values = array(data=NA, dim = c(niter, 1, length(nsamples)),
                         dimnames=list(1:niter, "censoring",nsamples))

for(j in 1:length(nsamples)){
  for (i in 1:niter) {
    set.seed(546+i)
    print(paste("Monte Carlo realization: ", i, ", sample size:", nsamples[j]))
    
    ## Gerar os dados
    data = gen.joint.cure.dgompertz(n=nsamples[j], alpha=alpha0,
                                    gammas=gammas0, betas.fixed=betas.fixed0,
                                    sd.int=sd.int0, sd.slope=sd.slo0, rho=rho0,
                                    gap=0.3, followup=0.9) #gap=0.6,followup=1.8
    
    survdata = data$survdata
    longdata = data$longdata
    
    
    # Ajuste do modelo de sobrevivência
    modS = inla(inla.surv(time, delta) ~ binXs + b_int + b_slo,
                data = survdata,
                family = "dgompertzsurv",
                control.inla = list(cmin = 0),
                control.fixed = list(prec.intercept = 1),
                control.family = list(hyper = list(alpha = list(param = c(0, 0.001)))))
    #summary(modS)
    
    
    # Ajuste do modelo conjunto
    modLS_0 = joint(formSurv = inla.surv(time, delta) ~ binXs,
                     dataSurv = survdata, basRisk = "weibullsurv", NbasRisk = 50,
                     formLong = biomarker ~ timej + cts + bin + (1 + timej | id),
                     family = c("poisson"), dataLong = longdata, id = "id",
                     timeVar = "time", assoc = c("SRE_ind"),
                     run = FALSE, control = list(int.strategy = "eb",
                                                 priorSRE_ind=list(mean=0,prec=1),
                                                 priorFixed=list(mean.intercept=0,prec.intercept=0.1,
                                                                 mean=0, prec=0.1)))
    
    modLS_0$.args$family = c("poisson", "dgompertzsurv")
    modLS_0$.args$control.family[[2]] = modS$.args$control.family
    modLS_0$.args$control.inla$cmin = 0
    modLS = joint.run(modLS_0, class = "inla")

    
    # sample of random effects
    MC_samples = inla.iidkd.sample(10^4, modLS, "IDIntercept_L1", return.cov = TRUE)
    ref_vals = t(sapply(MC_samples, function(Sigma) {
      cov_xy = Sigma[1, 2]
      sd_x = sqrt(Sigma[1, 1])
      sd_y = sqrt(Sigma[2, 2])
      cor = cov_xy / (sd_x * sd_y)
      c(sd_int = sd_x, sd_slo = sd_y, rho = cor)
    }))
    
    sd_int_est = ref_vals[,"sd_int"]
    sd_slo_est = ref_vals[,"sd_slo"]
    sd_rho_est = ref_vals[,"rho"]
    
    res = rbind(
      modLS$summary.hyperpar[1, c(1, 2, 3, 5)],
      modLS$summary.fixed[1, c(1, 2, 3, 5)],
      modLS$summary.hyperpar[5:6, c(1, 2, 3, 5)],
      modLS$summary.fixed[2, c(1, 2, 3, 5)],
      modLS$summary.fixed[-c(1,2), c(1, 2, 3, 5)],
      c(mean(sd_int_est), sd(sd_int_est), quantile(sd_int_est, probs = c(0.025, 0.975))),
      c(mean(sd_slo_est), sd(sd_slo_est), quantile(sd_slo_est, probs = c(0.025, 0.975))),
      c(mean(sd_rho_est), sd(sd_rho_est), quantile(sd_rho_est, probs = c(0.025, 0.975)))
    )
    
    check_ci = 1 * ((res$`0.025quant` <= parameters) & (res$`0.975quant` >= parameters))
    res$cp = check_ci
    res = cbind(parameters, res)
    bias = (res$mean - parameters)
    
    # Saving results
    res.matrix[i, 1:12, j] = res$mean
    res.matrix[i, 13:24, j] = res$sd
    res.matrix[i, 25:36, j] = bias
    res.matrix[i, 37:48, j] = res$cp
    censoring_values[i, , j] = 1 - mean(survdata$delta)
    
    if (i == niter) {
      write.csv2(res.matrix[,,j], file = paste("joint_model_res33_", "rep_", niter, "_n_", nsamples[j], ".csv", sep = ""))
      write.csv2(censoring_values[,,j], file = paste("joint_model_cens33_", "rep_", niter, "_n_", nsamples[j], ".csv", sep = ""))
    }
  }
}


resumo = as.matrix(colMeans(res.matrix[,,2]))

parameters
resumo[rownames(resumo)=="mean"]
resumo[rownames(resumo)=="sd"]
resumo[rownames(resumo)=="bias"]
resumo[rownames(resumo)=="cp"]


#res500=read.csv2("joint_model_res22_cen2_rep_1000_n_500.csv"); res500=res500[,-1]
#res100=read.csv2("joint_model_res22_cen2_rep_1000_n_100.csv"); res100=res100[,-1]
#res1000=read.csv2("joint_model_res22_cen2_rep_1000_n_1000.csv"); res1000=res1000[,-1]





