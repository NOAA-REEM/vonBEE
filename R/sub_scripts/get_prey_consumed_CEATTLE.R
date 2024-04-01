
# ConsumAge_fut g/pred/yr
#ration2Age_fut  = ceattle_2$ration2_1 kg/pred/yr
# avgN in thousands

load(file.path(run_fl,"data/in/dietdata/mnBy_YrBin.Rdata"))
load(file.path(run_fl,"data/in/dietdata/mnBy_Yr_JuvAdult.Rdata"))
load(file.path(run_fl,"data/in/dietdata/diet_p4CEATTLE.Rdata"))

assmnt_path <- file.path(run_fl,"docs/MSM_Assessment")
load(file.path(fl_2,"results/CEATTLE_results.Rdata")); ceattle_2<-tmp
TempC    <- as.numeric(ceattle_2$TempC)
Years    <- as.numeric(ceattle_2$styr)+1:as.numeric(ceattle_2$nyrs)-1
yearset  <- 1982:2012
mean(TempC[Years%in%yearset])
TempC_scaled <- (as.numeric(ceattle_2$TempC)-mean(TempC[Years%in%yearset]))/sqrt(sd(TempC[Years%in%yearset]))
nages    <- c(ceattle_2$nages_1,ceattle_2$nages_2,ceattle_2$nages_3)
nspp     <- as.numeric(ceattle_2$nspp)
spnm     <- c("walleye pollock","Pacific cod","arrowtooth flounder")
fnIN     <- file.path(run_fl,"data/in/dat_input_files/stomach2023.dat")
idsetUSE <- c("YEAR","REGION","SN","CN","juv2","annual_ration")
idsetUSE2 <- c("YEAR","REGION","SN","CN")
# E_est_all_1/1e6

# Biomass_eaten <- rbind(
#   data.frame(Year = yrs,Biomass=rowSums(E_est_all_1)/1e6,Species = "Prey: Walleye pollock"),
#   data.frame(Year = yrs,Biomass=rowSums(E_est_all_2)/1e6,Species = "Prey: Pacific cod"),
#   data.frame(Year = yrs,Biomass=rowSums(E_est_all_3)/1e6,Species = "Prey: Arrowtooth flounder"))



mean(TempC)

A2L_matrix <- list()
A2L_matrix[[1]] <-readxl::read_xlsx(file.path(run_fl,"data/in/dat_input_files/A2L_matrix.xlsx"),sheet="plk")
A2L_matrix[[2]] <-readxl::read_xlsx(file.path(run_fl,"data/in/dat_input_files/A2L_matrix.xlsx"),sheet="pcod")
A2L_matrix[[3]] <-readxl::read_xlsx(file.path(run_fl,"data/in/dat_input_files/A2L_matrix.xlsx"),sheet="atf")

l_bins      <- list()
l_bins[[1]] <- as.double(names(A2L_matrix[[1]]))
l_bins[[2]] <- as.double(names(A2L_matrix[[2]]))
l_bins[[3]] <- as.double(names(A2L_matrix[[3]]))


typeUSE <- "BW"
#tmpDATA <- mnBy_YrBin%>%
tmpDATA <- mnBy_Yr_JuvAdult%>%
  filter(!CN%in%"sablefish",type=="BW",PredL!=0,kind == "mean",REGION=="BS")%>%
  select(-c(kind,mode,DIG_cutoff
            ,C1_TWT,CNT,propDig_mnbyW,
            DIG_mnbyW,DIG,CNT,S0_TWT,W_use
            ,PredW,SST,BT,sum_propB))

tmpDATA$juv2 <- factor(c("Adult","Juvenile")[round(tmpDATA$juv,0)+1],
                       levels=c("Juvenile","Adult"))

#tmpDATA <- tmpDATA%>%mutate(BIN2 = getBIN2(BIN_cm_mid, BIN_IN = l_bins[[1]], type=1 ,divid=1))

sp<-1
#tmp_rat <- data.frame(matrix(NA, dim(ceattle_2$biomassByage_1)[1],lengths(l_bins)[sp]))
tmp_rat <- data.frame(matrix(NA, dim(ceattle_2$AvgN_1)[1],lengths(l_bins)[sp]))

colnames(tmp_rat)<-colnames(A2L_matrix[[sp]])
for(y in 1:dim(ceattle_2$AvgN_1)[1])
  tmp_rat[y,] <- colSums( (1000*ceattle_2$AvgN_1[y,]*ceattle_2$ration2_1[y,])*A2L_matrix[[sp]] )  #kg/yr
tmp_rat$YEAR <- 1:as.numeric(ceattle_2$nyrs)+as.numeric(ceattle_2$styr)-1
tmp_rat$CN   <- spnm[sp]
tmp_rat      <- tmp_rat%>%melt(id =c("CN","YEAR"),variable_name = "BIN2")%>%dplyr::rename(annual_ration = value)

annualrationByBin <- tmp_rat

sp <- 2
#tmp_rat <- data.frame(matrix(NA, dim(ceattle_2$AvgNe_2)[1],lengths(l_bins)[sp]))
tmp_rat <- data.frame(matrix(NA, dim(ceattle_2$AvgN_2)[1],lengths(l_bins)[sp]))

colnames(tmp_rat) <- colnames(A2L_matrix[[sp]])
for(y in 1:dim(ceattle_2$AvgN_2)[1])
  tmp_rat[y,] <- colSums( (1000*ceattle_2$AvgN_2[y,]*ceattle_2$ration2_2[y,])*A2L_matrix[[sp]] )
tmp_rat$YEAR <- 1:as.numeric(ceattle_2$nyrs)+as.numeric(ceattle_2$styr)-1
tmp_rat$CN   <- spnm[sp]
tmp_rat      <- tmp_rat%>%melt(id =c("CN","YEAR"),variable_name = "BIN2")%>%dplyr::rename(annual_ration = value)

annualrationByBin <- rbind(annualrationByBin,tmp_rat)

sp<-3
tmp_rat <- data.frame(matrix(NA, dim(ceattle_2$AvgN_3)[1],lengths(l_bins)[sp]))
colnames(tmp_rat)<-colnames(A2L_matrix[[sp]])
for(y in 1:dim(ceattle_2$AvgN_3)[1])
  tmp_rat[y,] <- colSums( (1000*ceattle_2$AvgN_3[y,]*ceattle_2$ration2_3[y,])*A2L_matrix[[sp]]) #kg/yr
tmp_rat$YEAR <- 1:as.numeric(ceattle_2$nyrs)+as.numeric(ceattle_2$styr)-1
tmp_rat$CN   <- spnm[sp]
tmp_rat      <- tmp_rat%>%melt(id =c("CN","YEAR"),variable_name = "BIN2")%>%dplyr::rename(annual_ration = value)

annualrationByBin <- rbind(annualrationByBin,tmp_rat)
annualrationByBin$BIN2 <- as.numeric(as.character(annualrationByBin$BIN2))
annualrationByBin$juv2 <- factor(c("Adult"),levels=c("Juvenile","Adult"))
annualrationByBin$juv2[annualrationByBin$BIN2<45] <- factor(c("Juvenile"),levels=c("Juvenile","Adult"))
annualrationJA <- annualrationByBin%>%group_by(CN,YEAR,juv2)%>%summarize(annual_ration= sum(annual_ration))

#tmpDATA <- tmpDATA%>%left_join(annualrationByBin)
tmpDATA  <- diet_p$data%>%left_join(annualrationJA)%>%
  mutate( ration =annual_ration* percent/100 )%>%
  dplyr::rename(total_ration_kg = annual_ration)

## now sum across juv and adult
annual_ration_JA <- tmpDATA
annual_ration    <- tmpDATA%>%group_by(across(all_of(c(idsetUSE2, "prey"))))%>%
  summarize(eaten = sum(ration,na.rm=T),total_ration_kg = sum(total_ration_kg,na.rm=T))

plot_prey_core = c("Arrowtooth","W..Pollock","P..Cod","Opilio","Unid.Chion","King.Crab","other","Euphausiid","Copepod","Salmon")

plot_aka_rat <- ggplot(annual_ration%>%filter(eaten!=0, prey%in%plot_prey_core))+
  geom_point(aes(x=YEAR,y=eaten*0.001/1e6,color=CN))+  #mil t / yr
  geom_line(aes(x=YEAR,y=eaten*0.001/1e6,color=CN))+
  facet_wrap(prey~.,scales="free_y",ncol=2)+
  theme_minimal()+ylab("eaten (million t per yr)")+
  expand_limits(y = 0)


#' now calculate the average diet composition for warm 
#' and cold years and expand that to annual rations:

  # block years
  upp  <- mean(TempC[Years%in%yearset])+0.5*sd(TempC[Years%in%yearset])
  dwn  <- mean(TempC[Years%in%yearset])-0.5*sd(TempC[Years%in%yearset])
  
  mn_diet <- diet_p$data%>%mutate(Yr_block=factor("Avg",levels=c("Warm","Avg","Cold")))
  mn_diet$Yr_block[mn_diet$YEAR%in%Years[TempC>=upp]] <- factor("Warm",levels=c("Warm","Avg","Cold"))
  mn_diet$Yr_block[mn_diet$YEAR%in%Years[TempC<=dwn]] <- factor("Cold",levels=c("Warm","Avg","Cold"))
  mn_diet <- mn_diet%>%group_by(REGION,SN,CN,juv2,prey,Yr_block)%>%
    summarize(percent = mean(percent,na.rm=T))
  mn_diet<- mn_diet%>%
    left_join(mn_diet%>%group_by(REGION,SN,CN,juv2,Yr_block)%>%
    summarize(total = sum(percent,na.rm=T)))%>%
    mutate(percent =100*percent/total)%>%select(-total)
  
  tmpDATA  <- mn_diet%>%left_join(annualrationJA,multiple = "all")%>%
    mutate( ration =annual_ration* percent/100 )%>%
    dplyr::rename(total_ration_kg = annual_ration)

## now sum across juv and adult
annual_ration_JA_allyr <- tmpDATA
annual_ration_allyr    <- tmpDATA%>%group_by(across(all_of(c(idsetUSE2, "prey"))))%>%
  summarize(eaten = sum(ration,na.rm=T),total_ration_kg = sum(total_ration_kg,na.rm=T))  

plot_prey_core = c("Arrowtooth","W..Pollock","P..Cod","Opilio","Unid.Chion","King.Crab","other","Euphausiid","Copepod")

plot_rat_allyr <- ggplot(annual_ration_allyr%>%filter(eaten!=0, prey%in%plot_prey_core))+
  geom_point(aes(x=YEAR,y=eaten*0.001/1000,color=CN))+
  geom_line(aes(x=YEAR,y=eaten*0.001/1000,color=CN))+
  facet_wrap(prey~.,scales="free_y",ncol=2)+
  theme_minimal()+ylab("eaten (thousand t per yr)")+
  expand_limits(y = 0)

dat1 <- annual_ration_allyr%>%filter((eaten*0.001/1000)>0.005, prey%in%plot_prey_core)
dat2 <- annual_ration%>%filter((eaten*0.001/1000)>0.005, prey%in%plot_prey_core)


plot_rat_all <- ggplot()+
  geom_point(data=dat1,aes(x=YEAR,y=eaten*0.001/1000,color=CN, shape="mndiet"))+  #eaten = 1000 kg/yr  t/yr
  geom_line(data=dat1,aes(x=YEAR,y=eaten*0.001/1000,color=CN, linetype="mndiet"))+
  geom_point(data=dat2,aes(x=YEAR,y=eaten*0.001/1000,color=CN, shape="obs diet"))+
  geom_line(data=dat2,aes(x=YEAR,y=eaten*0.001/1000,color=CN, linetype="obs diet"))+
  facet_wrap(prey~.,scales="free_y",ncol=2)+
  theme_minimal()+ylab("Biomass eaten ( x 1000 t)")+
  expand_limits(y = 0)+scale_color_viridis_d(end=0,begin=.8)

# jpeg(file = file.path(plot_file,"Prey_eaten_annually.jpg"), width = 7 , height=8,unit="in",res=250)
# print(plot_rat_all)
# dev.off()

jpeg(file = file.path(assmnt_path,"Figs/Prey_eaten_annually.jpg"), width = 7 , height=8,unit="in",res=250)
print(plot_rat_all)
dev.off()

save(annual_ration_allyr,file=file.path(assmnt_path,"Figs/annual_ration_allyr.Rdata"))
save(annual_ration,file=file.path(assmnt_path,"Figs/annual_ration.Rdata"))

