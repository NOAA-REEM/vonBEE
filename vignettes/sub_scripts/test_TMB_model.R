#'
#'
#'
#'test_TMB_model.R
#'
#'This code tests the TMB model
#'
#'

out_fn <- file.path(prefn,"data","out","test_EBS")
if(!dir.exists(out_fn)) dir.create(out_fn)

# select subset for EBS pollock
dat <- dat.all%>%filter(Sp==species[1],REGION==regions[1])%>%
  mutate(cohortYr   = YEAR-age,
         cohortYrfac = (factor(YEAR-age)))

cov_dat <- dat%>%select(Temp, YEAR)%>%group_by(YEAR)%>%
  summarize(mnTemp = mean(Temp,na.rm=T),
            seTemp = se (Temp,na.rm=T))%>%
  mutate(factYr = factor(YEAR),
         YrIndex = as.numeric(factYr))%>%
  relocate(YrIndex,YEAR, mnTemp, seTemp,factYr)%>%arrange(YrIndex)%>%
  ungroup()

cohort_cov_dat <- dat%>%select(Temp, cohortYr)%>%
  group_by(cohortYr)%>%
  summarize(chort_mnTemp = mean(Temp,na.rm=T),
            chort_seTemp = se (Temp,na.rm=T))%>%
  ungroup()

dat     <- dat%>%left_join(cov_dat)%>%left_join(cohort_cov_dat)
covs    <- as.matrix(data.frame(blank=0,mnTemp=dat$mnTemp,mnTemp2=dat$mnTemp^2,chort_mnTemp=dat$chort_mnTemp))
covs_se <- as.matrix(data.frame(blank=0,seTemp=dat$seTemp,seTemp2 =dat$seTemp^2,chort_seTemp = dat$chort_seTemp))
ncov    <- dim(covs)[2]
nyrs    <- dim(cov_dat)[1]
nobs    <- dim(dat)[1]
nyrs_cohort<- length(levels(dat$cohortYrfac))

datIN <- list(YrIndex    = dat$YrIndex-1, 
              age        = dat$age, 
              weight     = dat$weight, 
              cohortYr   = as.numeric(dat$cohortYrfac)-1, 
              ncov       = ncov,
              nobs       = nobs,
              nyrs       = nyrs,
              d_offset   = 0 ,#0.01,  # prevents WInf --> 0
              #covYr      = cov_dat$YrIndex,
              vonb_cov   = matrix(t(covs),ncov,nobs),
              sdvonb_cov = matrix(t(covs_se),ncov,nobs))

save(dat,file = file.path(out_fn,"dat.Rdata"))
save(covs,file = file.path(out_fn,"covs.Rdata"))
save(covs_se,file = file.path(out_fn,"covs_se.Rdata"))

# run sim
# ---------------------------------------
parIN <- parIN2<- list(
  logH   = log(12.0),
  logK   = log(0.15),
  t0  = -0.25,
  mu_d = -.005,
  log_sigma_proc   = -.6,
  log_sigma_obs   = -.3,
  beta_c     = (c(0,0,0,0)),
  beta_1     = (0),
  beta_0     = (rep(0,nyrs_cohort)),
  u_y        = (rep(0,nyrs)))

parIN2$beta_c <- (c(0,0.03,-0.0074,0))

simd    <- sim_dat(par = parIN, data = datIN)
simd_T2 <- sim_dat(par = parIN2, data = datIN)
names(simd_T2$out)[-(1:3)]<-paste0("cov",names(simd_T2$out)[-(1:3)])

x<-1:20  
y <-0.063*x+-0.0094*x^2
plot(x,y)
p1<- ggplot()+
  geom_point(data=simd$out,aes(x=age,y=W_obs,color="observed"))+
  geom_point(data=simd_T2$out,aes(x=age,y=W_hat,color=factor(round(cov2))))+
  geom_point(data=simd$out,aes(x=age,y=W_hat,color="predicted"))+theme_minimal()+ scale_color_brewer(palette = "Spectral")
p1
sclr<-1.5
jpeg(filename = file.path(out_fn,"sim_plot.jpg"), res = 350, width = sclr*6, height = sclr*4,units = "in")
print(p1)
dev.off()

# Now set up VonBT models
#---------------------------------------
## TEST the model with simulated data
dataIN_test        <- datIN
dataIN_test$weight <- simd$What

# define starting par
truepar <- testpar <- parIN

# define different starting points for fitting the model
testpar$logH   = log(2.0)
testpar$logK   = log(0.25)
testpar$mu_d   = 1.01

# map it
tmp_map <- makeMap(param=truepar,estpar=list(
  logH   = TRUE,
  logK   = TRUE,
  t0     = TRUE,
  mu_d   = TRUE,
  log_sigma_proc  = FALSE,
  log_sigma_obs   = TRUE,
  beta_c     = FALSE,
  beta_1     = FALSE,
  beta_0     = FALSE,
  u_y        = FALSE))

for(i in 1:length(dataIN_test)){
  if(any(is.na(dataIN_test[[i]])))
    stop(cat("problem ! NA found in dataIN_test: ",names(dataIN_test)[i]))
}
dlistIN_true <- list(data     = dataIN_test, 
                     parameters  = truepar, 
                     random      = NULL,  
                     map         = tmp_map)
dlistIN_test <- list(data    = dataIN_test, 
                     parameters  = testpar, 
                     random      = NULL,  
                     map         = tmp_map)

true_model  <- MakeADFun(
  data                 =  dlistIN_true$data,
  parameters           =  dlistIN_true$parameters,
  random               =  dlistIN_true$random,
  DLL                  =  "VonBTv1",
  checkParameterOrder  =  TRUE,
  hessian              =  TRUE,
  map                  =  dlistIN_true$map,
  silent               =  TRUE)

test_model  <- MakeADFun(
  data                 =  dlistIN_test$data,
  parameters           =  dlistIN_test$parameters,
  random               =  dlistIN_test$random,
  DLL                  =  "VonBTv1",
  checkParameterOrder  =  TRUE,
  hessian              =  TRUE,
  map                  =  dlistIN_test$map,
  silent               =  TRUE)

set.seed(1) ## optional
tt_true  <- true_model$simulate()
tt_test  <- test_model$simulate()

p <- ggplot()+
  geom_point(data=data.frame(age=dataIN_test$age,weight=dataIN_test$weight),aes(x=age,y = weight,color="a obs"),size=6)+
  geom_point(data=data.frame(age=tt_true$age,weight=tt_true$weight,What=tt_true$What),aes(x=age,y = weight,color="b TMB_Whobs"),size=4)+
  geom_point(data=data.frame(age=tt_true$age,weight=tt_true$weight,What=tt_true$What),aes(x=age,y = What,color="c TMB_What"),size=2)+
  geom_point(data=data.frame(age=tt_test$age,weight=tt_test$weight,What=tt_test$What),aes(x=age,y = What,color="d TMB_test_What_start_val"),size=.5)+
  theme_minimal()+scale_colour_manual( values = c(col_ramp(3),"gray","gray"))
p

sclr<-1.5
jpeg(filename = file.path(out_fn,"test_start.jpg"), res = 350, width = sclr*6, height = sclr*4,units = "in")
print(p)
dev.off()

test_model$gr(test_model$par)
test_model$fn(test_model$par)
test_model$fn(true_model$par)
test_model$gr(true_model$par)

tmpmod <- list()
true_model$par
tmpmod$fit      <-   nlminb(start = test_model$par, obj = test_model$fn, gr = test_model$gr, 
                            control=list(iter.max=10000,eval.max=10000))
print(test_model$env$last.par.best)

for(i in 1:3){
  tmpmod$fit      <-  nlminb(start = test_model$env$last.par.best,
                             obj = test_model$fn,
                             gr = test_model$gr, 
                             control=list(iter.max=10000,eval.max=10000))
  print(test_model$env$last.par.best)
}
true_model$par
test_model$par

summary(sdreport(test_model, getJointPrecision = TRUE)) 
test_model$report()

# p + geom_point(data=data.frame(What= test_model$report()$What, age= test_model$report()$age),
#                aes(x=age,y=What, color="TMB Test"))

m_test0 <- suppressWarnings(
  runMod(dlistIN    = dlistIN_test,
         src_fldr   = getwd(),
         version    = "VonBTv1",
         nrep       = 3,
         hessianIN  = TRUE,
         silentIN   = FALSE,
         se.fit     = TRUE,
         simulate   = FALSE,
         sim_nitr   = 1000,
         recompile  = FALSE,
         maxitr     = 10000,
         maxeval    = 10000))

# plot fitted to true values
p + geom_point(data=m_test0$fitted,
               aes(x=age,y=W_hat, color="TMB Test"),size=.3)
p
sclr<-1.5
jpeg(filename = file.path(out_fn,"test_end.jpg"), res = 350, width = sclr*6, height = sclr*4,units = "in")
print(p + geom_point(data=m_test0$fitted,
                     aes(x=age,y=W_hat, color="TMB Test"),size=.3))
dev.off()


