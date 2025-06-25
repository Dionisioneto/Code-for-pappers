## application to data with shared random effects
## in defective gompertz distribution

library(survival)
library(joineRML)
library(INLA)
library(INLAjoint)

set.seed(1)

data=joineRML::epileptic.qol
data$time <- data$time/365.25
data$with.time <- data$with.time/365.25

epileptic.l <- data[, -c(2,4,9)]
epileptic.s = data[!duplicated(data$id), c(1:3, 9)] # filtering the survival time and censoring
head(epileptic.l)
head(epileptic.s)

## Analyzing the characteristics of the biomarkers
## for the poisson, mean=variance
mean(epileptic.l$anxiety,na.rm=T);var(epileptic.l$anxiety,na.rm=T)
mean(epileptic.l$aep,na.rm=T);var(epileptic.l$aep,na.rm=T) # breaks the assumption
mean(epileptic.l$depress,na.rm=T);var(epileptic.l$depress,na.rm=T)


# Inla estimation of the survival
survpart = inla(inla.surv(with.time, with.status2) ~ 1,
         data = epileptic.s,
         family = "dgompertzsurv",
         control.inla = list(cmin = 0),
         control.family=list(hyper=list(alpha=list(param=c(0, 0.01)))),
         verbose = TRUE)
#summary(survpart) ## defective model


## Phase 1: fitting just one longitudinal biomarker

## (1L) anxiety
shared_long1 = joint(formSurv=inla.surv(with.time, with.status2) ~ 1,
                    dataSurv=epileptic.s, basRisk="weibullsurv", NbasRisk=50,
                    formLong=anxiety ~ time + trt + (1+time|id),
                    family=c("poisson"), dataLong=epileptic.l, id="id",
                    timeVar="time", assoc=c("SRE_ind"),
                    run=F, control=list(int.strategy="eb"))

shared_long1$.args$family <- c("poisson", "dgompertzsurv")
shared_long1$.args$control.family[[2]] <- survpart$.args$control.family
shared_long1$.args$control.inla$cmin <- 0
modLS1 <- joint.run(shared_long1, class="inla")
sum1 = summary(modLS1)

# estimation of the random effects
MC_samples1 <- inla.iidkd.sample(10^4, modLS1,
                                 "IDIntercept_L1", return.cov=TRUE) # for random-effects covariance terms
VarCov1 <- matrix(unlist(MC_samples1), nrow = 2^2)

VarCovMeans1 <- matrix(rowMeans(VarCov1),2,2);VarCovMeans1
# sd and quantiles over the estimate:
VarCovSD_1 <- matrix(apply(VarCov1,1,sd),2,2);VarCovSD_1
VarCov025_1 <- matrix(apply(VarCov1,1,function(x) quantile(x, 0.025)),2,2);VarCov025_1
VarCov975_1 <- matrix(apply(VarCov1,1,function(x) quantile(x, 0.975)),2,2);VarCov975_1

## extrair os vlores das amostras mcmc da correlação
cor_vals1 <- sapply(MC_samples1, function(Sigma) {
  cov_xy <- Sigma[1, 2]
  sd_x <- sqrt(Sigma[1, 1])
  sd_y <- sqrt(Sigma[2, 2])
  cov_xy / (sd_x * sd_y)
})

# mean(cor_vals1)
# sd(cor_vals1)
# quantile(cor_vals1, probs = c(0.025,0.975))

## parameters = c(alpha,gamma0,gamma1,gamma2,
##              beta0, beta1, beta2, sd_int, sd_slo, rho)

## posterior means
posterior_means1 = c(sum1$hyperpar[1,"mean"],
                    c(sum1$fixed[1,"mean"], sum1$hyperpar[c(5,6),"mean"]),
                    as.vector(sum1$fixed[-1,"mean"]),
                    sqrt(VarCovMeans1[1,1]), sqrt(VarCovMeans1[2,2]),
                    mean(cor_vals1))

## posterior sds
post_sds = c(sum1$hyperpar[1,"sd"],
             c(sum1$fixed[1,"sd"], sum1$hyperpar[c(5,6),"sd"]),
             as.vector(sum1$fixed[-1,"sd"]),
             sqrt(VarCovSD_1[1,1]),
             sqrt(VarCovSD_1[2,2]), sd(cor_vals1))

## credibility interval
cimatrix1 = as.matrix(rbind(sum1$hyperpar[1,c("0.025quant","0.975quant")],
                            sum1$fixed[1,c("0.025quant","0.975quant")],
                            sum1$hyperpar[c(5,6),c("0.025quant","0.975quant")],
                            sum1$fixed[-1,c("0.025quant","0.975quant")],
                            c(sqrt(VarCov025_1[1,1]),sqrt(VarCov975_1[1,1])),
                            c(sqrt(VarCov025_1[2,2]),sqrt(VarCov975_1[2,2])),
                            quantile(cor_vals1, probs = c(0.025,0.975))))

rownames(cimatrix1)=c("alpha","gamma0", "gamma1", "gamma2",
                     "beta0","beta1","beta2",
                     "sd.int","sd.slope","rho")

cimatrix1

## (1L) depress
shared_long2 = joint(formSurv=inla.surv(with.time, with.status2) ~ 1,
                    dataSurv=epileptic.s, basRisk="weibullsurv", NbasRisk=50,
                    formLong=depress ~ time + trt + (1+time|id),
                    family=c("poisson"), dataLong=epileptic.l, id="id",
                    timeVar="time", assoc=c("SRE_ind"),
                    run=F, control=list(int.strategy="eb"))

shared_long2$.args$family <- c("poisson", "dgompertzsurv")
shared_long2$.args$control.family[[2]] <- survpart$.args$control.family
shared_long2$.args$control.inla$cmin <- 0
modLS2 <- joint.run(shared_long2, class="inla")

sum2 = summary(modLS2)

# estimation of the random effects
MC_samples2 <- inla.iidkd.sample(10^4, modLS2,
                                 "IDIntercept_L1", return.cov=TRUE) # for random-effects covariance terms
VarCov2 <- matrix(unlist(MC_samples2), nrow = 2^2)

VarCovMeans2 <- matrix(rowMeans(VarCov2),2,2);VarCovMeans2
# sd and quantiles over the estimate:
VarCovSD_2 <- matrix(apply(VarCov2,1,sd),2,2);VarCovSD_2
VarCov025_2 <- matrix(apply(VarCov2,1,function(x) quantile(x, 0.025)),2,2);VarCov025_2
VarCov975_2 <- matrix(apply(VarCov2,1,function(x) quantile(x, 0.975)),2,2);VarCov975_2

## extrair os vlores das amostras mcmc da correlação
cor_vals2 <- sapply(MC_samples2, function(Sigma) {
  cov_xy <- Sigma[1, 2]
  sd_x <- sqrt(Sigma[1, 1])
  sd_y <- sqrt(Sigma[2, 2])
  cov_xy / (sd_x * sd_y)
})



## posterior means
posterior_means2 = c(sum2$hyperpar[1,"mean"],
                     c(sum2$fixed[1,"mean"], sum2$hyperpar[c(5,6),"mean"]),
                     as.vector(sum2$fixed[-1,"mean"]),
                     sqrt(VarCovMeans2[1,1]), sqrt(VarCovMeans2[2,2]),
                     mean(cor_vals2))

## posterior sds
post_sds = c(sum2$hyperpar[1,"sd"],
             c(sum2$fixed[1,"sd"], sum2$hyperpar[c(5,6),"sd"]),
             as.vector(sum1$fixed[-1,"sd"]),
             sqrt(VarCovSD_2[1,1]),
             sqrt(VarCovSD_2[2,2]), sd(cor_vals2))

## credibility interval
cimatrix2 = as.matrix(rbind(sum2$hyperpar[1,c("0.025quant","0.975quant")],
                            sum2$fixed[1,c("0.025quant","0.975quant")],
                            sum2$hyperpar[c(5,6),c("0.025quant","0.975quant")],
                            sum2$fixed[-1,c("0.025quant","0.975quant")],
                            c(sqrt(VarCov025_2[1,1]),sqrt(VarCov975_2[1,1])),
                            c(sqrt(VarCov025_2[2,2]),sqrt(VarCov975_2[2,2])),
                            quantile(cor_vals2, probs = c(0.025,0.975))))

rownames(cimatrix2)=c("alpha","gamma0", "gamma1", "gamma2",
                      "beta0","beta1","beta2",
                      "sd.int","sd.slope","rho")

cimatrix2



## Phase 2: fitting two longitudinal biomarkers
## (2L) depress & aep

shared_long3 = joint(formSurv=inla.surv(with.time, with.status2) ~ 1,
                    dataSurv=epileptic.s, basRisk="weibullsurv", NbasRisk=50,
                    formLong=list(depress ~ time + trt + (1+time|id),
                                  aep ~ time + trt + (1+time|id)),
                    family=c("poisson","poisson"), dataLong=epileptic.l, id="id",
                    timeVar="time", assoc=c("SRE_ind","SRE_ind"), corLong=T,
                    run=F, control=list(int.strategy="eb"))

shared_long3$.args$family <- c("poisson", "poisson", "dgompertzsurv")
shared_long3$.args$control.family[[3]] <- survpart$.args$control.family
shared_long3$.args$control.inla$cmin <- 0
modLS3 <- joint.run(shared_long3)
SDC <- summary(modLS3, sdcor=T)$ReffList[[1]];SDC # <- aqui
class(modLS3) <- "inla"

sum3=summary(modLS3)
sum3


# estimation of the random effects
# MC_samples3 <- inla.iidkd.sample(10^4, modLS3,
#                                  "IDIntercept_L1", return.cov=TRUE) # for random-effects covariance terms
# VarCov3 <- matrix(unlist(MC_samples3), nrow = 2^2)
#
# VarCovMeans3 <- matrix(rowMeans(VarCov3),2,2);VarCovMeans3
# # sd and quantiles over the estimate:
# VarCovSD_3 <- matrix(apply(VarCov3,1,sd),2,2);VarCovSD_3
# VarCov025_3 <- matrix(apply(VarCov3,1,function(x) quantile(x, 0.025)),2,2);VarCov025_3
# VarCov975_3 <- matrix(apply(VarCov3,1,function(x) quantile(x, 0.975)),2,2);VarCov975_3
#
# ## extrair os vlores das amostras mcmc da correlação
# cor_vals3 <- sapply(MC_samples3, function(Sigma) {
#   cov_xy <- Sigma[1, 2]
#   sd_x <- sqrt(Sigma[1, 1])
#   sd_y <- sqrt(Sigma[2, 2])
#   cov_xy / (sd_x * sd_y)
# })
#
# sd_vals3 <- sapply(MC_samples3, function(Sigma) {
#   sd_x <- sqrt(Sigma[1, 1])
# })
#
# sd_vals3 <- sapply(MC_samples3, function(Sigma) {
#   sd_x <- sqrt(Sigma[1, 1])
# })



## parameters = c(alpha,gamma0,gamma1,gamma2,gamma3,gamma4,
##              betay0, betay1, betay2,
##              betaz0, betaz1, betaz2,
##              sd_int, sd_slo, rho)


## posterior means
# posterior_means3 = c(sum3$hyperpar[1,"mean"],
#                     c(sum3$fixed[1,"mean"],sum3$hyperpar[12:15,"mean"]),
#                     as.vector(sum3$fixed[-1,"mean"]),
#                     sqrt(VarCovMeans3[1,1]), sqrt(VarCovMeans3[2,2]),
#                     mean(cor_vals3))
#
# ## posterior sds
# post_sds3 = c(sum3$hyperpar[1,"sd"],
#             c(sum3$fixed[1,"sd"],sum3$hyperpar[12:15,"sd"]),
#             as.vector(sum3$fixed[-1,"sd"]),
#             sqrt(VarCovSD_3[1,1]), sqrt(VarCovSD_3[2,2]),
#             sd(cor_vals3))
#
# ## credibility interval
#
# cimatrix3 = as.matrix(rbind(sum3$hyperpar[1,c("0.025quant","0.975quant")],
#                             sum3$fixed[1,c("0.025quant","0.975quant")],
#                             sum3$hyperpar[12:15,c("0.025quant","0.975quant")],
#                             sum3$fixed[-1,c("0.025quant","0.975quant")],
#                             c(sqrt(VarCov025_3[1,1]),sqrt(VarCov975_3[1,1])),
#                             c(sqrt(VarCov025_3[2,2]),sqrt(VarCov975_3[2,2])),
#                             quantile(cor_vals3, probs = c(0.025,0.975))))
#
# rownames(cimatrix2)=c("alpha","gamma0", "gamma1", "gamma2",
#                       "beta0","beta1","beta2",
#                       "sd.int","sd.slope","rho")
#
#
#









