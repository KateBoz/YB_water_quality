##############################################
# R code for Sigleo and Bosley isotope paper
#
# Created by: Katelyn Bosley
# Date: 6/2/2018
##############################################

# read in the data
setwd("E:\\Sigleo_Bosley_MS")

dat<-read.csv("Y_iso_data_fixed.csv")
head(dat)


# Load in packages needed for the analysis
load.libraries<-function(){
  library(ggplot2)
  library(lme4)
  library(dplyr)
  library(GGally)
}
  
load.libraries()

#fix up the station codes
dat$Station_alt<-as.factor(dat$Station_alt)
unique(dat$Station_alt)


#################################################
#make exploratory plots with some basic analyses
#################################################


#setting up the theme for ggplots
cust_theme<-
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.background = element_rect(fill="transparent"),
        panel.border = element_rect(colour = "black"))+
  #theme(legend.title = element_blank())+
  theme(strip.text.x = element_text(size = 10, colour = "black", face="bold"))



#############################################################################################
# Start with a simple matrix of all the relationships and see if anything interesting emerges

mat.dat<-dat[,c(3,8,9,10,11,12,13,14)]

#generate pdf with plots
pdf("YB_pairs.pdf",paper='letter') 
ggpairs(mat.dat)
dev.off()



# 







##########################
#Salinity

#build.df
sal.df<-dat[,c("SampleID","Station_alt","Year","Month","Day","Salinity_ppt")]

#calculate site means by month in each year
sal.ag<-sal.df %>% group_by(Station_alt,Year)%>%summarize(SampleN=length(Salinity_ppt),mean_sal=mean(Salinity_ppt,na.rm=T),sd_sal=sd(Salinity_ppt, na.rm=T))
#head(sal.ag)
sal.ag$Year<-as.factor(sal.ag$Year)
sal.ag<-na.omit(sal.ag)

#make a nice plot!

sal.p<-ggplot(sal.ag, aes(x=Station_alt,y=mean_sal,group=Year))+
  geom_point(aes(col=Year, shape=Year),size = 3)+
  geom_line(aes(col=Year),linetype = 2)+
  scale_color_manual(name="Year", values = c("grey20","grey50","grey70"))+
  ylab("Mean Salinity")+
  xlab("Station_ID")+
  cust_theme+
  ggtitle("Mean Annual Salinity")


# add in the monthly variation
#build.df
sal.ag.mo<-sal.df %>% group_by(Station_alt,Year,Month)%>%summarize(SampleN=length(Salinity_ppt),mean_sal=mean(Salinity_ppt,na.rm=T))

sal.ag.mo$Year<-as.factor(sal.ag.mo$Year)
sal.ag.mo$Month<-as.factor(sal.ag.mo$Month)


sal.mo.p<-ggplot(sal.ag.mo, aes(x=Station_alt,y=mean_sal,group=Month))+
  geom_point(aes(col=Month),size = 3, shape=8)+
  geom_line(aes(col=Month),linetype = 1,size=1.2)+
  #scale_color_manual(name="Year", values = c("grey20","grey50","grey70"))+
  ylab("Mean Salinity")+
  xlab("Station_ID")+
  cust_theme+
  ggtitle("Mean Monthly Salinity")+
  facet_grid(Year~.)


#simple model to look at differences among months and years

sal.mod<-aov(sal.ag.mo$mean_sal~sal.ag.mo$Station_alt+sal.ag.mo$Month+sal.ag.mo$Year)
summary(sal.mod)
# different by station and Month but not year


#########################################
#repeat the above plot with other values
##########################################
#silica

sil.df<-dat[,c("SampleID","Station_alt","Year","Month","Day","Silica")]

# add in the monthly variation

sil.ag.mo<-sil.df %>% group_by(Station_alt,Year,Month)%>%summarize(SampleN=length(Silica),mean_sil=mean(Silica,na.rm=T))

sil.ag.mo$Year<-as.factor(sil.ag.mo$Year)
sil.ag.mo$Month<-as.factor(sil.ag.mo$Month)


sil.mo.p<-ggplot(sil.ag.mo, aes(x=Station_alt,y=mean_sil,group=Month))+
  geom_point(aes(col=Month),size = 3, shape=8)+
  geom_line(aes(col=Month),linetype = 1,size=1.2)+
  #scale_color_manual(name="Year", values = c("grey20","grey50","grey70"))+
  ylab("Mean Silica conc")+
  xlab("Station_ID")+
  cust_theme+
  ggtitle("Mean Monthly Silica")+
  facet_grid(Year~.)


sil.mod<-aov(sil.ag.mo$mean_sil~sil.ag.mo$Station_alt+sil.ag.mo$Month+sil.ag.mo$Year)
summary(sil.mod)
# different by station and Month but not year


#merge sil and sal for combined plot
m.master<-merge(sil.ag.mo,sal.ag.mo)
#head(m.master)

sil.sal.p<-ggplot(m.master, aes(x=mean_sal,y=mean_sil,group=Month))+
  geom_point(aes(col=Month),size = 3, shape=8)+
  geom_line(aes(col=Month),linetype = 1,size=1)+
  #scale_color_manual(name="Year", values = c("grey20","grey50","grey70"))+
  ylab("Mean Silica conc")+
  xlab("Mean Salinity (ppt)")+
  cust_theme+
  ggtitle("Salinity vs. Silica")+
  facet_grid(Year~.)


# model
sil.sal.mod<-aov(m.master$mean_sil~m.master$mean_sal+m.master$Month*m.master$Year)
summary(sil.sal.mod)

#################################
#P
p.df<-dat[,c("SampleID","Station_alt","Year","Month","Day","P")]

# add in the monthly variation

p.ag.mo<-p.df %>% group_by(Station_alt,Year,Month)%>%summarize(SampleN=length(P),mean_P=mean(P,na.rm=T))

p.ag.mo$Year<-as.factor(p.ag.mo$Year)
p.ag.mo$Month<-as.factor(p.ag.mo$Month)


p.mo.p<-ggplot(p.ag.mo, aes(x=Station_alt,y=mean_P,group=Month))+
  geom_point(aes(col=Month),size = 3, shape=8)+
  geom_line(aes(col=Month),linetype = 1,size=1.2)+
  #scale_color_manual(name="Year", values = c("grey20","grey50","grey70"))+
  ylab("Mean Silica conc")+
  xlab("Station_ID")+
  cust_theme+
  ggtitle("Mean Monthly Phosporus")+
  facet_grid(Year~.)


p.mod<-aov(p.ag.mo$mean_P~p.ag.mo$Station_alt+p.ag.mo$Month+p.ag.mo$Year)
summary(p.mod)
# different by station and Month but not year

plot(resid(p.mod)) #looking really good!


#add P to master
m.master<-merge(m.master,p.ag.mo)


# everything is realated to salinity
p.sal.p<-ggplot(m.master, aes(x=mean_sal,y=log(mean_P),group=Month))+
  geom_point(aes(col=Month),size = 3, shape=8)+
  geom_line(aes(col=Month),linetype = 1,size=1)+
  #scale_color_manual(name="Year", values = c("grey20","grey50","grey70"))+
  ylab("Mean P conc")+
  xlab("Mean Salinity (ppt)")+
  cust_theme+
  ggtitle("Salinity vs. Phos")+
  facet_grid(Year~.)


# model
p.sal.mod<-aov(log(m.master$mean_P)~m.master$mean_sal+m.master$Month*m.master$Year)
summary(p.sal.mod)
plot(resid(p.sal.mod))

############################
#Nitrate
############################
#P
n.df<-dat[,c("SampleID","Station_alt","Year","Month","Day","Ni_Nitrate")]

# add in the monthly variation

n.ag.mo<-n.df %>% group_by(Station_alt,Year,Month)%>%summarize(SampleN=length(Ni_Nitrate),mean_Ni=mean(Ni_Nitrate,na.rm=T))

n.ag.mo$Year<-as.factor(n.ag.mo$Year)
n.ag.mo$Month<-as.factor(n.ag.mo$Month)


n.mo.p<-ggplot(n.ag.mo, aes(x=Station_alt,y=mean_Ni,group=Month))+
  geom_point(aes(col=Month),size = 3, shape=8)+
  geom_line(aes(col=Month),linetype = 1,size=1.2)+
  #scale_color_manual(name="Year", values = c("grey20","grey50","grey70"))+
  ylab("Mean Nitrate conc")+
  xlab("Station_ID")+
  cust_theme+
  ggtitle("Mean Monthly Nitrate")+
  facet_grid(Year~.)


n.mod<-aov(n.ag.mo$mean_Ni~n.ag.mo$Station_alt+n.ag.mo$Month+n.ag.mo$Year)
summary(n.mod)
# different by station and Month but not year

plot(resid(n.mod)) #looking really good!


#add P to master
m.master<-merge(m.master,n.ag.mo)


# everything is realated to salinity
n.sal.p<-ggplot(m.master, aes(x=mean_sal,y=log(mean_Ni),group=Month))+
  geom_point(aes(col=Month),size = 3, shape=8)+
  geom_line(aes(col=Month),linetype = 1,size=1)+
  #scale_color_manual(name="Year", values = c("grey20","grey50","grey70"))+
  ylab("Mean Nitrate conc")+
  xlab("Mean Salinity (ppt)")+
  cust_theme+
  ggtitle("Salinity vs. Nitrate")+
  facet_grid(Year~.)


#model
n.sal.mod<-aov(log(m.master$mean_Ni)~m.master$mean_sal+m.master$Month+m.master$Year)
summary(n.sal.mod)
plot(resid(n.sal.mod))


###################################
# ISOTOPES
###################################
#del13C
delC.df<-dat[,c("SampleID","Station_alt","Year","Month","Day","del_13C")]

# add in the monthly variation

delC.ag.mo<-delC.df %>% group_by(Station_alt,Year,Month)%>%summarize(SampleN=length(del_13C),mean_delC=mean(del_13C,na.rm=T))

delC.ag.mo$Year<-as.factor(delC.ag.mo$Year)
delC.ag.mo$Month<-as.factor(delC.ag.mo$Month)

delC.mo.p<-ggplot(delC.ag.mo, aes(x=Station_alt,y=mean_delC,group=Month))+
  geom_point(aes(col=Month),size = 3, shape=8)+
  geom_line(aes(col=Month),linetype = 1,size=1.2)+
  #scale_color_manual(name="Year", values = c("grey20","grey50","grey70"))+
  ylab("Mean del 13 C")+
  xlab("Station_ID")+
  cust_theme+
  ggtitle("Mean Monthly del 13C")+
  facet_grid(Year~.)


delC.mod<-aov(delC.ag.mo$mean_delC~delC.ag.mo$Station_alt+delC.ag.mo$Month+delC.ag.mo$Year)
summary(delC.mod)
# different by station and Month but not year

plot(resid(delC.mod)) #looking really good!


#add P to master
m.master<-merge(m.master,delC.ag.mo)


# everything is realated to salinity
delC.sal.p<-ggplot(m.master, aes(x=mean_sal,y=mean_delC,group=Month))+
  geom_point(aes(col=Month),size = 3, shape=8)+
  geom_line(aes(col=Month),linetype = 1,size=1)+
  #scale_color_manual(name="Year", values = c("grey20","grey50","grey70"))+
  ylab("Mean del 13C")+
  xlab("Mean Salinity (ppt)")+
  cust_theme+
  ggtitle("Salinity vs. del 13C")+
  facet_grid(Year~.)


#model
delC.sal.mod<-aov(m.master$mean_delC~m.master$mean_sal+m.master$Month+m.master$Year)
summary(delC.sal.mod)
plot(resid(delC.sal.mod))

###############
# del 15N
###############

#del15N
delN.df<-dat[,c("SampleID","Station_alt","Year","Month","Day","del_15N")]

# add in the monthly variation

delN.ag.mo<-delN.df %>% group_by(Station_alt,Year,Month)%>%summarize(SampleN=length(del_15N),mean_delN=mean(del_15N,na.rm=T))

delN.ag.mo$Year<-as.factor(delN.ag.mo$Year)
delN.ag.mo$Month<-as.factor(delN.ag.mo$Month)

delN.mo.p<-ggplot(delN.ag.mo, aes(x=Station_alt,y=mean_delN,group=Month))+
  geom_point(aes(col=Month),size = 3, shape=8)+
  geom_line(aes(col=Month),linetype = 1,size=1.2)+
  #scale_color_manual(name="Year", values = c("grey20","grey50","grey70"))+
  ylab("Mean del 15 N")+
  xlab("Station_ID")+
  cust_theme+
  ggtitle("Mean Monthly del 15N")+
  facet_grid(Year~.)


delN.mod<-aov(delN.ag.mo$mean_delN~delN.ag.mo$Station_alt+delN.ag.mo$Month+delN.ag.mo$Year)
summary(delN.mod)

# different by station and Month but not year

plot(resid(delN.mod)) #looking really good!


#add P to master
m.master<-merge(m.master,delN.ag.mo)


# everything is realated to salinity
delN.sal.p<-ggplot(m.master, aes(x=mean_sal,y=mean_delN,group=Month))+
  geom_point(aes(col=Month),size = 3, shape=8)+
  geom_line(aes(col=Month),linetype = 1,size=1)+
  #scale_color_manual(name="Year", values = c("grey20","grey50","grey70"))+
  ylab("Mean del 15N")+
  xlab("Mean Salinity (ppt)")+
  cust_theme+
  ggtitle("Salinity vs. del 15N")+
  facet_grid(Year~.)


#model
delN.sal.mod<-aov(m.master$mean_delN~m.master$mean_sal+m.master$Month+m.master$Year)
summary(delN.sal.mod)
plot(resid(delN.sal.mod))


####################################################################################

#isotope biplots?

#complete the Siber analysis  
plot(dat$del_13C~dat$del_15N)


# present zooplankton data relative to the seston
# traditional biplot
#


######################
#Multivariate Analysis
#####################

