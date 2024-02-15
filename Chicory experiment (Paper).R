######======= Chicory (Cichorium intybus) reduces cyathostomin egg excretion and development in grazing horses  ========
## == Packages ==
library(ggplot2)
library(RColorBrewer)
theme_set(theme_bw())
require(geepack)
require(dplyr)
library(gridExtra)
library(showtext)
font_add_google("Lato")
library(ggridges)
library(ggpubr)
library(rstatix)
require(plyr)
library(reshape2)
library(tidyverse)
library(purrr)
require(drc)
require(ade4)
require(vegan)
require(ggrepel)
require(circlize)
library(tidyr)
require(scales)
library(nlme)
library("knitr")
library("BiocStyle")
require('dada2')
require('phyloseq')
library(lme4)
library(data.table)
require(stringi)
require(stringr)
library(readr)
library(DECIPHER)
library(ShortRead)
library(Biostrings)
library(Hmisc)
require(car)
library(vegan)
library(patchwork)

data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)}

get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

## === Work file ===
path<-"~/Scripts/Chicory paper"
setwd(path)

########## Herbage characteristics and dietary choices ##########
# === Body weight data ====
BW=read.csv(file='Weight_XP_Chicory.csv',header=T,sep=';',dec=',')
BW$Group<- factor(BW$Group, levels = c('Control','Chicory'))

BW_summary=data_summary(BW, varname="Weight",groupnames=c("Day",'Group'))
BW_summary
#   Day   Group  Weight       sd
# 1 d-5 Control 457.558 18.18193
# 2 d-5 Chicory 460.561 20.21653
# 3  d0 Control 456.565 17.92421
# 4  d0 Chicory 461.523 20.20675
# 5 d45 Control 471.242 20.45362
# 6 d45 Chicory 458.896 17.55711

# === Statistical analysis ====
## BW to create each group
BW_5=BW[BW$Day=='d-5',]
#Identify outliers
BW_5 %>% 
  group_by(Group) %>%
  identify_outliers(Weight)
# Group   Name             Horses Day   Weight is.outlier is.extreme
# 1 Chicory KASSCROUT DE SPI      2 d-5     436. TRUE       FALSE     
# 2 Chicory KECHAZULA DE SPI      3 d-5     486. TRUE       FALSE     
# 3 Chicory KELELAPOM DE SPI      6 d-5     421. TRUE       TRUE      
# 4 Chicory KONTENTE DE SPI       9 d-5     487. TRUE       FALSE 

#Normality
shapiro_test(residuals(lm(Weight ~ Group, data = BW_5)))
#   variable                                 statistic p.value
# 1 residuals(lm(Weight ~ Group, data = BW_5))      0.937   0.212

BW_5 %>%
  group_by(Group) %>%
  shapiro_test(Weight)
#    Group   variable statistic     p
# 1 Control Weight       0.856 0.0676
# 2 Chicory Weight       0.906 0.255 

#Equality of variance
BW_5 %>% levene_test(Weight ~ Group)
#    df1   df2 statistic     p
#   1    18   0.00172 0.967

#ANOVA
res.aov <- BW_5 %>% anova_test(Weight ~ Group)
res.aov
# ANOVA Table (type II tests)
#   Effect DFn DFd     F     p p<.05   ges
# 1  Group   1  18 0.122 0.731       0.007

## BW during the experiment
BW_test=BW[BW$Day!="d-5",]

BW_test %>% 
  group_by(Group) %>%
  identify_outliers(Weight)

#Normality
shapiro_test(residuals(lm(Weight ~ Group, data = BW_test)))
#   variable                                      statistic p.value
# 1 residuals(lm(Weight ~ Group, data = BW_test))     0.984   0.840

BW_test %>%
  group_by(Group) %>%
  shapiro_test(Weight)
#  Group   variable statistic     p
# 1 Control Weight       0.951 0.376
# 2 Chicory Weight       0.976 0.865

#Equality of variance
BW_test %>% levene_test(Weight ~ Group)
#   df1   df2 statistic     p
#   1     1    38     0.133 0.717

#ANOVA
BW_cov=BW[BW$Day!="d-5",]
BW_cov=BW_cov %>% spread(Day,Weight)
res.aov_cov <- BW_cov %>% anova_test(d45 ~ Group+d0)
res.aov_cov
# ANOVA Table (type II tests)
# 
# Effect DFn DFd       F        p p<.05   ges
# 1  Group   1  17  40.567 6.99e-06     * 0.705
# 2     d0   1  17 168.668 2.97e-10     * 0.908

########## Chicory effect on FEC ##########
# === Parasite data ====
Chico = read.csv(file='Parasite_data_XP_Chicory.csv',header=T,sep=';',dec=',')
Chico$Group<- factor(Chico$Group, levels = c('Control','Chicory'))

#Mean after log+1 transformation
Chico$logEPG_plus1=log(Chico$EPG+1)
FEC_summary=data_summary(Chico, varname="logEPG_plus1",groupnames=c("Day",'Group'))
FEC_summary
#   Day   Group logEPG_plus1        sd
# 1  d0 Control     7.607597 0.4195966
# 2  d0 Chicory     7.612341 0.4131267
# 3 d16 Control     6.931347 0.3663331
# 4 d16 Chicory     5.329088 1.0678714
# 5 d31 Control     7.000914 0.5090608
# 6 d31 Chicory     4.859432 1.8811817
# 7 d45 Control     6.695623 0.3997935
# 8 d45 Chicory     3.907118 1.8811513

#Mean after back-transformation
Chico$backEPG_plus1=exp(Chico$logEPG_plus1)-1
FEC_summary=data_summary(Chico, varname="backEPG_plus1",groupnames=c("Day",'Group'))
FEC_summary
#   Day   Group FEC_summary$Low       sd
# 1  d0 Control        2169.0 853.5573
# 2  d0 Chicory        2169.0 794.7390
# 3 d16 Control        1084.5 383.0176
# 4 d16 Chicory         285.0 175.7840
# 5 d31 Control        1233.0 654.8206
# 6 d31 Chicory         282.0 307.9610
# 7 d45 Control         864.0 316.9543
# 8 d45 Chicory         153.0 215.3834

FEC_summary$error= qnorm(0.975)*FEC_summary$sd/sqrt(10)
FEC_summary$Low = FEC_summary$backEPG_plus1 - FEC_summary$error
FEC_summary$High =FEC_summary$backEPG_plus1 + FEC_summary$error
FEC_summary
#   Day   Group backEPG_plus1       sd    error        Low      High
# 1  d0 Control        2169.0 853.5573 529.0306 1639.96945 2698.0306
# 2  d0 Chicory        2169.0 794.7390 492.5753 1676.42475 2661.5753
# 3 d16 Control        1084.5 383.0176 237.3924  847.10758 1321.8924
# 4 d16 Chicory         285.0 175.7840 108.9500  176.04998  393.9500
# 5 d31 Control        1233.0 654.8206 405.8545  827.14545 1638.8545
# 6 d31 Chicory         282.0 307.9610 190.8727   91.12728  472.8727
# 7 d45 Control         864.0 316.9543 196.4467  667.55333 1060.4467
# 8 d45 Chicory         153.0 215.3834 133.4935   19.50645  286.4935

#Plot
ggplot(Chico, aes(x=Day, y=EPG, fill=Group))+
  geom_boxplot(alpha=0.4)+
  geom_point(aes(x = Day,y= EPG, group = Group), size = 1.5, shape = 1,position = position_jitterdodge(0))+
  labs(title=, y='Fecal egg count (Eggs per gram)', x='Day')+
  theme(axis.line = element_line(size = 1, linetype = "solid"),
        legend.title = element_blank(),legend.position="bottom",
        legend.text = element_text(size=33, family = "Lato"),
        axis.title.x = element_text(size=35, family = "Lato", margin = margin(t = 0.4, unit="cm")), 
        axis.title.y = element_text(size=35, family = "Lato", margin = margin(r = 0.4, unit="cm")),
        axis.text.x = element_text(size=25, face = "bold", family = "Lato"),
        axis.text.y = element_text(size=25, face='bold', family = "Lato"))+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))+
  scale_fill_manual(values=c('#737373','#006d2c'))


# === Statistical analysis ====
gee.fit.EPG <-glmer.nb(EPG ~ Group * Day + (1|Horses), data = Chico)
summary(gee.fit.EPG)
# Fixed effects:
#                     Estimate Std. Error z value Pr(>|z|)    
# (Intercept)           7.6550     0.2462  31.094  < 2e-16 ***
# GroupChicory          0.0746     0.3517   0.212 0.832006    
# Dayd16               -0.7033     0.2862  -2.457 0.013996 *  
# Dayd31               -0.5903     0.2882  -2.048 0.040533 *  
# Dayd45               -0.9277     0.2870  -3.232 0.001230 ** 
# GroupChicory:Dayd16  -1.5133     0.4124  -3.670 0.000243 ***
# GroupChicory:Dayd31  -1.8378     0.4204  -4.372 1.23e-05 ***
# GroupChicory:Dayd45  -2.2261     0.4309  -5.166 2.39e-07 ***
#   ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

D0Con=Chico$EPG[Chico$Day=='d0' & Chico$Group=="Control"]
D0Chico=Chico$EPG[Chico$Day=='d0' & Chico$Group=="Chicory"]
wilcox.test(D0Con,D0Chico)
# data:  D0Con and D0Chico
# W = 49, p-value = 0.9698
# alternative hypothesis: true location shift is not equal to 0


D16Con=Chico$EPG[Chico$Day=='d16' & Chico$Group=="Control"]
D16Chico=Chico$EPG[Chico$Day=='d16' & Chico$Group=="Chicory"]
wilcox.test(D16Con,D16Chico)
#data:  D16Con and D16Chico
# W = 100, p-value = 0.0001817
# alternative hypothesis: true location shift is not equal to 0 

D31Con=Chico$EPG[Chico$Day=='d31' & Chico$Group=="Control"]
D31Chico=Chico$EPG[Chico$Day=='d31' & Chico$Group=="Chicory"]
wilcox.test(D31Con,D31Chico)
#data:  D31Con and D31Chico
# W = 94, p-value = 0.0003248
# alternative hypothesis: true location shift is not equal to 0

D45Con=Chico$EPG[Chico$Day=='d45' & Chico$Group=="Control"]
D45Chico=Chico$EPG[Chico$Day=='d45' & Chico$Group=="Chicory"]
wilcox.test(D45Con,D45Chico)
# data:  D45Con and D45Chico
# W = 97, p-value = 0.0004353
# alternative hypothesis: true location shift is not equal to 0


# gee.fit.EPG <- geeglm(EPG ~ Group * Day,id = Horses, data = Chico, family = poisson,
#                       corstr = "ar1", scale.fix = TRUE, std.err = "san.se")
# summary(gee.fit.EPG)
# # Coefficients:
# #                         Estimate   Std.err    Wald Pr(>|W|)    
# # (Intercept)          7.68e+00  1.18e-01 4234.11  < 2e-16 ***
# # GroupChicory        -5.77e-18  1.61e-01    0.00  1.00000    
# # DayD16              -6.93e-01  1.59e-01   19.09  1.2e-05 ***
# # DayD31              -5.65e-01  1.98e-01    8.11  0.00439 ** 
# # DayD45              -9.20e-01  1.61e-01   32.52  1.2e-08 ***
# # GroupChicory:DayD16 -1.34e+00  2.67e-01   24.98  5.8e-07 ***
# # GroupChicory:DayD31 -1.48e+00  3.98e-01   13.71  0.00021 ***
# # GroupChicory:DayD45 -1.73e+00  4.65e-01   13.84  0.00020 ***
# # ---
# # Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#==== Chicory efficacy ====
###----- Bayesian hierarchical models
fec_ctl0 = Chico$EPG[Chico$Group=='Control' & Chico$Day=='d0']
fec_ctl16 = Chico$EPG[Chico$Group=='Control' & Chico$Day=='d16']
fec_ctl31 = Chico$EPG[Chico$Group=='Control' & Chico$Day=='d31']
fec_ctl45 = Chico$EPG[Chico$Group=='Control' & Chico$Day=='d45']
fec_chico0 = Chico$EPG[Chico$Group=='Chicory' & Chico$Day=='d0']
fec_chico16 = Chico$EPG[Chico$Group=='Chicory' & Chico$Day=='d16']
fec_chico31 = Chico$EPG[Chico$Group=='Chicory' & Chico$Day=='d31']
fec_chico45 = Chico$EPG[Chico$Group=='Chicory' & Chico$Day=='d45']

mod_0 <- eggCounts::fecr_stan(fec_ctl0, fec_chico0, rawCounts = TRUE, preCF=50, postCF=50,
                              paired = TRUE, indEfficacy = TRUE)
mod_0$posterior.summary
#                         mean        sd      2.5%       50%      97.5%  HPDLow95       mode  HPDHigh95
# FECR                  0.2068    0.0791    0.0953     0.193     0.3968    0.0839     0.1724     0.3622
# meanEPG.untreated 12741.1135 3649.9504 6766.5243 12412.380 20667.9531 6664.5898 11202.4075 20221.1567
# meanEPG.treated   10104.5014 3079.6644 5100.3281  9716.054 16956.6893 4985.4587  8969.9611 16695.1345

mod_16 <- eggCounts::fecr_stan(fec_ctl16, fec_chico16, rawCounts = TRUE, preCF=50, postCF=50,
                               paired = TRUE, indEfficacy = TRUE)
mod_16$posterior.summary
#                         mean       sd      2.5%        50%      97.5% HPDLow95       mode  HPDHigh95
# FECR                  0.7288    0.078    0.5396     0.7366     0.8575    0.559     0.7568     0.8665
# meanEPG.untreated 14310.2089 4102.247 7529.6886 13877.2577 23310.6459 6865.438 12947.4600 22383.9747
# meanEPG.treated    3867.3025 1541.893 1591.5577  3628.3201  7515.5530 1285.703  3057.3578  7035.7192

mod_31 <- eggCounts::fecr_stan(fec_ctl31, fec_chico31, rawCounts = TRUE, preCF=50, postCF=50,
                               paired = TRUE, indEfficacy = TRUE)
mod_31$posterior.summary
#                         mean        sd      2.5%        50%      97.5%  HPDLow95       mode  HPDHigh95
# FECR                  0.7853    0.0965    0.5464     0.8014     0.9315    0.5771     0.8337     0.9472
# meanEPG.untreated 13384.1253 3867.7946 6870.4659 13060.4077 21439.5923 6532.4142 12527.7072 20753.9412
# meanEPG.treated    2877.4335 1587.3500  775.0548  2547.6081  6866.3950  519.7044  1994.3112  6093.5680

mod_45 <- eggCounts::fecr_stan(fec_ctl45, fec_chico45, rawCounts = TRUE, preCF=50, postCF=50,
                               paired = TRUE, indEfficacy = TRUE)
mod_45$posterior.summary
#                         mean        sd      2.5%        50%      97.5%  HPDLow95       mode  HPDHigh95
# FECR                  0.8548    0.0831    0.6475     0.8706     0.9746    0.6951     0.8997     0.9952
# meanEPG.untreated 14763.1220 4104.1911 7506.4874 14501.1554 23427.1651 7262.2810 13118.9021 23008.4416
# meanEPG.treated    2147.3106 1402.3396  301.0928  1807.4226  5754.2500   71.8540  1448.0010  4769.5084

###----- Cumulative FEC analysis
data_FEC=select(Chico, Horses, Group,Day,logEPG_plus1)
data_FEC=data_FEC %>% spread(Day,logEPG_plus1)

data_FEC$mean_FEC_d0d16=((data_FEC$d0+data_FEC$d16)/2)*16
data_FEC$mean_FEC_d16d31=((data_FEC$d16+data_FEC$d31)/2)*15
data_FEC$mean_FEC_d31d45=((data_FEC$d31+data_FEC$d45)/2)*14

data_FEC$FEC_0plus16=data_FEC$d0+data_FEC$mean_FEC_d0d16
data_FEC$FEC_0plus16plus31=data_FEC$FEC_0plus16+data_FEC$mean_FEC_d16d31
data_FEC$FEC_0plus16plus31plus45=data_FEC$FEC_0plus16plus31+data_FEC$mean_FEC_d31d45

FEC_cum=select(data_FEC, Horses, Group,FEC_0plus16 ,FEC_0plus16plus31,FEC_0plus16plus31plus45)
FEC_cum=FEC_cum %>% pivot_longer(cols = c("FEC_0plus16","FEC_0plus16plus31","FEC_0plus16plus31plus45"), names_to = "day", values_to = "cumul")
mean_FEC_cum=data_summary(FEC_cum, varname="cumul",groupnames=c('Group','day'))
mean_FEC_cum
#     Group                     day    cumul        sd
# 1 Control             FEC_0plus16 28197.00  9554.488
# 2 Control       FEC_0plus16plus31 45578.25 14740.692
# 3 Control FEC_0plus16plus31plus45 60257.25 19843.162
# 4 Chicory             FEC_0plus16 21801.00  7880.940
# 5 Chicory       FEC_0plus16plus31 26053.50 10336.284
# 6 Chicory FEC_0plus16plus31plus45 29098.50 13076.311

cumul_ctl45 = FEC_cum$cumul[FEC_cum$Group=='Control' & FEC_cum$day=='FEC_0plus16plus31plus45']
cumul_chi45 = FEC_cum$cumul[FEC_cum$Group=='Chicory' & FEC_cum$day=='FEC_0plus16plus31plus45']

wilcox.test(cumul_ctl45,cumul_chi45, alternative = "two.sided")
# Wilcoxon rank sum exact test
# 
# data:  cumul_ctl45 and cumul_chi45
# W = 96, p-value = 0.0001299
# alternative hypothesis: true location shift is not equal to 0

########## Chicory effect on larval development rate ##########
Chico$id=row.names(Chico)
OUT=Chico %>% 
  group_by(Group, Day) %>%
  identify_outliers(Dev)
dtout=data.frame(Group=OUT$Group,Day=OUT$Day, Name=OUT$Name, Horses=OUT$Horses,
                 Date=OUT$Date, EPG=OUT$EPG, EPG_corrected=OUT$EPG_corrected, 
                 Dev=OUT$Dev,MF_g=OUT$MF_g, MS_g=OUT$MS_g,id=OUT$id)
dtout=dtout[dtout$Day!="d0",]
Chico=anti_join(Chico,dtout, by = "id")

Dev_summary_byDay=data_summary(Chico, varname="Dev",groupnames=c("Day"))
Dev_summary_byDay
# Day   Dev    sd
# 1  D0  1.01  2.45
# 2 D16 12.77  6.80
# 3 D31 16.34 14.91
# 4 D45 14.54 10.53

Dev_summary_byDay_Group=data_summary(Chico, varname="Dev",groupnames=c("Day","Group"))
Dev_summary_byDay_Group
#   Day   Group   Dev    sd
# 1  D0 Control  0.67  1.07
# 2  D0 Chicory  1.35  3.35
# 3 D16 Control 13.17  7.52
# 4 D16 Chicory 12.38  6.37
# 5 D31 Control 21.06 11.08
# 6 D31 Chicory  7.69  6.25
# 7 D45 Control 21.60  9.21
# 8 D45 Chicory  7.48  6.20

#Plot
ggplot(Chico, aes(x=Day, y=Dev, fill=Group))+
  geom_boxplot(alpha=0.4)+
  geom_point(aes(x = Day,y= Dev, group = Group), size = 1.5, shape = 1,position = position_jitterdodge(0))+
  labs(title=, y='Larval development rate (%)', x='Day')+
  theme(axis.line = element_line(size = 1, linetype = "solid"),
        legend.title = element_blank(),legend.position="bottom",
        legend.text = element_text(size=33, family = "Lato"),
        axis.title.x = element_text(size=35, family = "Lato", margin = margin(t = 0.4, unit="cm")), 
        axis.title.y = element_text(size=35, family = "Lato", margin = margin(r = 0.4, unit="cm")),
        axis.text.x = element_text(size=25, face = "bold", family = "Lato"),
        axis.text.y = element_text(size=25, face = "bold", family = "Lato"))+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))+
  scale_fill_manual(values=c('#737373','#006d2c'))

# === Statistical analysis ====
NoD0=Chico[Chico$Day!="d0",]

gee.fit.Dev <- geeglm(Dev ~ Group * Day, id = Horses, data = NoD0, family = poisson,
                      corstr = "ar1", scale.fix = TRUE, std.err = "san.se")
summary(gee.fit.Dev)
# Coefficients:
# Estimate  Std.err    Wald Pr(>|W|)    
# (Intercept)          2.57794  0.17130 226.477  < 2e-16 ***
# GroupChicory        -0.06226  0.23063   0.073  0.78719    
# Dayd31               0.46929  0.23813   3.884  0.04876 *  
# Dayd45               0.49458  0.21382   5.350  0.02072 *  
# GroupChicory:Dayd31 -0.94528  0.37430   6.378  0.01155 *  
# GroupChicory:Dayd45 -0.99777  0.36253   7.575  0.00592 ** 
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

ChicoD16= Chico$Dev[Chico$Group=='Chicory' & Chico$Day=='d16']
CtlD16= Chico$Dev[Chico$Group=='Control' & Chico$Day=='d16']
wilcox.test(ChicoD16,CtlD16)
# data:  ChicoD16 and CtlD16
# W = 48, p-value = 0.9118
# alternative hypothesis: true location shift is not equal to 0

ChicoD31=NoD0[NoD0$Group=="Chicory" & NoD0$Day=="d31",]
ChicoD31=ChicoD31$Dev
CtlD31=NoD0[NoD0$Group=="Control" & NoD0$Day=="d31",]
CtlD31=CtlD31$Dev

wilcox.test(ChicoD31,CtlD31)
# data:  ChicoD31 and CtlD31
# W = 12, p-value = 0.006
# alternative hypothesis: true location shift is not equal to 0

ChicoD45=NoD0[NoD0$Group=="Chicory" & NoD0$Day=="d45",]
ChicoD45=ChicoD45$Dev
CtlD45=NoD0[NoD0$Group=="Control" & NoD0$Day=="d45",]
CtlD45=CtlD45$Dev

wilcox.test(ChicoD45,CtlD45)
# data:  ChicoD45 and CtlD45
# W = 6, p-value = 0.001
# alternative hypothesis: true location shift is not equal to 0

wilcox.test(CtlD16,CtlD31)
# data:  CtlD16 and CtlD31
# W = 25, p-value = 0.1
# alternative hypothesis: true location shift is not equal to 0

wilcox.test(CtlD31,CtlD45)
# data:  CtlD31 and CtlD45
# W = 44, p-value = 1
# alternative hypothesis: true location shift is not equal to 0


########## Chicory effect on equine cyathostomins larval community structure ##########
tr=read.table("track_output_chicory_dada2_R_ci_mxee25_trunc200_BS16.tsv")
head(tr)
summary(tr$filtered)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 820  402444  458629  465366  548932  703166 

###----====== 
set.seed(100) # Initialize random number generator for reproducibility

### Training set - once and for all
taxmeth = 'idtaxa'
train <- readDNAStringSet("curated_nemaca_ITS2_Strongylidae.fasta") 
tax <- read_tsv("idtaxa_03022022.tax") 

trainingSet <- LearnTaxa(train, names(train), tax)

### Read in seqtab
seqtab_nochim = read.table(file = paste0('output_chicory_dada2_R_ci_mxee25_trunc200_BS16.tsv'),
                           header=T, sep='\t')

dim(seqtab_nochim)
#[1] 4047   61

dna <- DNAStringSet(getSequences(seqtab_nochim$OTUID))
seqtab_nochim$OTUID = paste0('ASV_',1:nrow(seqtab_nochim))
rownames(seqtab_nochim) = seqtab_nochim$OTUID
names(dna) = seqtab_nochim$OTUID
seqtab_nochim$OTUID = NULL

colnames(seqtab_nochim) = sapply(stringr::str_split(colnames(seqtab_nochim),'\\.'),
                                 function(x) gsub('WP3','',x[2]))

###-------------- Taxonomy assignment
trainingSet <- LearnTaxa(train, names(train), tax)

ids <- IdTaxa(dna,
              trainingSet,
              strand = "both",
              threshold = 50,
              bootstraps = 100,
              processors = NULL,
              verbose = TRUE,
              type = "extended")

ranks <- c("kingdom","phylum", "class", "order", "family", "genus", "species") # ranks of interest
# Convert the output object of class "Taxa" to a matrix analogous to the output from assignTaxonomy
taxidITS <- t(sapply(ids, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa
}))
taxidITS = gsub('_',' ',taxidITS)
colnames(taxidITS) <- ranks; rownames(taxidITS) <- names(dna) #getSequences(seqtab_nochim)

#### Metadata
metadata = data.frame(sample.id = colnames(seqtab_nochim))
rownames(metadata) = metadata$sample.id
metadata$Group = substr(metadata$sample.id,1,2)
metadata$Group[metadata$Group=='Ch'] = 'Chicory'
metadata$Group[metadata$Group=='Te'] = 'Control'
metadata$Group = factor(metadata$Group, levels = c('Control','Chicory','Mo'))

metadata$Horse = sapply(stringr::str_split(colnames(seqtab_nochim),'D'),
                        function(x) x[1])
metadata$Day = sapply(stringr::str_split(colnames(seqtab_nochim),'D'),
                      function(x) paste0('d',x[2]))
metadata$Day[metadata$Day=='dNA']='Mock'

head(metadata)
#       sample.id   Group Horse Day
# Ch1D0     Ch1D0 Chicory   Ch1   0
# Ch2D0     Ch2D0 Chicory   Ch2   0
# Te4D0     Te4D0 Control   Te4   0
# Te5D0     Te5D0 Control   Te5   0
# Te6D0     Te6D0 Control   Te6   0
# Te7D0     Te7D0 Control   Te7   0

table(metadata$Day,metadata$Group)
#      Control Chicory
# 0          4       2
# 16         8      10
# 31        10       8
# 45         9       5
# Mock       0       0

### A few horses have complete data b/c of missing samples at day0
table(metadata$Horse)
# Ch1  Ch10   Ch2   Ch3   Ch4   Ch5   Ch6   Ch7   Ch8   Ch9 Mock1 Mock2 Mock3 Mock4  Te11  Te12 
#   4     1     3     3     3     3     2     2     2     2     1     1     1     1     3     3 
# Te13  Te14  Te15  Te16  Te17  Te18  Te19  Te20   Te4   Te5   Te6   Te7 
#    3     2     2     3     3     3     2     3     1     1     1     1 

table(metadata$Horse,metadata$Day)

### Phyloseq object
ps = phyloseq(
  otu_table(seqtab_nochim, taxa_are_rows = TRUE),
  tax_table(taxidITS),
  sample_data(metadata))

## Species level only
psITS = tax_glom(ps,taxrank='species',NArm=F)

### Remove samples that did not work
summary(colSums(otu_table(psITS)))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 231  386325  431192  432364  496499  674697 

quantile(colSums(otu_table(psITS)))
plot(density(colSums(otu_table(psITS))))
abline(v=80000)

which(colSums(otu_table(psITS))<100000)
# Te5D0 Te13D16 
# 4      19 
not_worked = names(which(colSums(otu_table(psITS))<100000))

### Remove contaminants according to count data
psITS.f =  filter_taxa(psITS, function(x) sum(x) > 10000, TRUE)
psITS.f
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 15 taxa and 60 samples ]
# sample_data() Sample Data:       [ 60 samples by 4 sample variables ]
# tax_table()   Taxonomy Table:    [ 15 taxa by 7 taxonomic ranks ]

### Remove samples that failed
psITS.f.wk = subset_samples(psITS.f,!(sample.id %in% not_worked))
psITS.f.wk
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 15 taxa and 58 samples ]
# sample_data() Sample Data:       [ 58 samples by 4 sample variables ]
# tax_table()   Taxonomy Table:    [ 15 taxa by 7 taxonomic ranks ]

### Transform counts
psITS.f.wk.t = transform_sample_counts(psITS.f.wk,
                                       function(OTU) OTU/sum(OTU))
psITS.f.wk.t
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 15 taxa and 58 samples ]
# sample_data() Sample Data:       [ 58 samples by 4 sample variables ]
# tax_table()   Taxonomy Table:    [ 15 taxa by 7 taxonomic ranks ]

### Remove contaminants from relative abundances & based on mock
conta=names(which(rowSums(otu_table(psITS.f.wk.t))<= 0.2))
conta
# "ASV_137"

psITS.f.wk.t.f =  filter_taxa(psITS.f.wk.t, function(x) sum(x) > 0.2, TRUE)
psITS.f.wk.t.f
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 14 taxa and 58 samples ]
# sample_data() Sample Data:       [ 58 samples by 4 sample variables ]
# tax_table()   Taxonomy Table:    [ 14 taxa by 7 taxonomic ranks ]

### Missing rate ; species level only
nalist_sp = rownames(tax_table(psITS.f.wk.t.f)[which(is.na(tax_table(psITS.f.wk.t.f)[,7]))])
nalist_sp
#[1] "ASV_12" "ASV_25"
which(is.na(tax_table(psITS.f.wk.t.f)[,6]))
#integer(0)

df_na = otu_table(psITS.f.wk.t.f)[which(rownames(otu_table(psITS.f.wk.t.f)) %in% nalist_sp),]
df_na = df_na[,which(colSums(df_na)>0)]
df_na = reshape2::melt(df_na)
summary(df_na$value)
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0000000 0.0000801 0.0009602 0.0344833 0.0195691 0.4523214 

### Assign missing taxonomy at the genus level
tax_table(psITS.f.wk.t.f)[,6][which(is.na(tax_table(psITS.f.wk.t.f)[,6]))]=paste0(tax_table(psITS.f.wk.t.f)[which(is.na(tax_table(psITS.f.wk.t.f)[,6])),5])
tax_table(psITS.f.wk.t.f)[,7][which(is.na(tax_table(psITS.f.wk.t.f)[,7]))]=paste0(tax_table(psITS.f.wk.t.f)[which(is.na(tax_table(psITS.f.wk.t.f)[,7])),6],'_sp')

aggregate(value ~ Var1, FUN = sum,data = df_na)

mis.sum = aggregate(value ~ Var2, FUN = sum,data = df_na)
mis.sum[order(mis.sum$value,decreasing = T),]

### Remove mock
psci =  subset_samples(psITS.f.wk.t.f,Day !='Mock')
psci
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 14 taxa and 54 samples ]
# sample_data() Sample Data:       [ 54 samples by 4 sample variables ]
# tax_table()   Taxonomy Table:    [ 14 taxa by 7 taxonomic ranks ]

### Overall distribution
plot_bar(psci, x = 'Horse', fill="species") +
  facet_wrap(~ Group + Day , ncol = 4, scales = 'free_x') + 
  theme_classic() +
  theme(legend.position = 'bottom',text = element_text(size = 16),
        axis.text.x = element_text(size = 10),
        legend.title = element_blank(),
        strip.background = element_blank()) + 
  guides(fill = guide_legend(nrow = 3))

### Remove d0 and keep horses with 3 occurrences
psci.f = subset_samples(psci,Day !='d0')
Htokeep=names(table(sample_data(psci.f)$Horse)[which(table(sample_data(psci.f)$Horse)==3)])
ci = subset_samples(psci.f,Horse %in% Htokeep)

nb.cols <- 14
mycolors = colorRampPalette(brewer.pal(8, "Set2"))(nb.cols)

Abundance=plot_bar(ci, x = 'Horse', fill="species") +
  scale_fill_manual(values = mycolors)+ 
  facet_wrap(~ Group + Day , ncol = 3, scales = 'free_x') + 
  theme_classic() +
  theme(legend.position = 'bottom',text = element_text(size = 16),
        axis.text.x = element_text(size = 10),
        legend.title = element_blank(),
        strip.background = element_blank()) + 
  theme(axis.line = element_line(size = 1, linetype = "solid"),
        legend.title = element_blank(),legend.position="bottom",
        legend.text = element_text(size=33, family = "Lato",face = "italic"),
        axis.title.x = element_text(size=38, family = "Lato", margin = margin(t = 0.4, unit="cm")), 
        axis.title.y = element_text(size=38, family = "Lato", margin = margin(r = 0.4, unit="cm")),
        axis.text.x = element_text(size=29, face = "bold", family = "Lato"),
        axis.text.y = element_text(size=29, face='bold', family = "Lato"))+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))+
  theme(strip.text.x = element_text(size=29, family = "Lato", face="bold"))+
  guides(fill = guide_legend(nrow = 4))
Abundance

###Differential study on count
psITchico = subset_samples(psITS.f.wk,Group %in% c('Chicory','Control') & Day %in% c("d16","d45"))
psITchico
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 15 taxa and 31 samples ]
# sample_data() Sample Data:       [ 31 samples by 4 sample variables ]
# tax_table()   Taxonomy Table:    [ 15 taxa by 7 taxonomic ranks ]

table(sample_data(psITchico)$Group, sample_data(psITchico)$Day)
# d16 d45
# Control   7   9
# Chicory  10   5

# ### Trajectory from day16 to day45
ct=psITS.f.wk
ct = subset_samples(ct,Day !='d0' & Day!="d31" & Group!='Mo')

dfTraj = reshape2::melt(as.matrix(otu_table(ct)))
colnames(dfTraj) = c('otu','sample.id','count')
dfTraj$Horse = sample_data(ct)$Horse[match(dfTraj$sample.id,sample_data(ct)$sample.id)]
dfTraj$Group = sample_data(ct)$Group[match(dfTraj$sample.id,sample_data(ct)$sample.id)]
idx = match(dfTraj$otu,rownames(tax_table(ct)))
dfTraj$species = paste0(gsub('s__','',tax_table(ct)[idx,7]))
dfTraj$species = factor(dfTraj$species)
dfTraj$Group = factor(dfTraj$Group)
dfTraj$Id = factor(sapply(str_sub(dfTraj$Horse,3), function(x) x[1]))
dfTraj$Day = sapply(str_split(dfTraj$sample.id,"D"), function(x) x[2])
dfTraj$Day = sapply(str_pad(dfTraj$Day, width = 3, side = "left", pad = "D"), function(x) x[1])
dfTraj$Day=as.factor(dfTraj$Day)
mTraj = nlme::lme(sqrt(sqrt(count)) ~ species*Group*Day ,
                  random=~ 1|Id,
                  data = dfTraj)
summary(mTraj)
#                                                               Value Std.Error  DF   t-value p-value
# (Intercept)                                                3.836957  1.970660 395  1.947041  0.0522
# speciesCoronocyclus labiatus                              -0.913588  2.758280 395 -0.331217  0.7407
# speciesCyathostomum catinatum                              3.406367  2.758280 395  1.234960  0.2176
# speciesCyathostomum pateratum                              7.697939  2.758280 395  2.790847  0.0055
# speciesCylicocyclus ashworthi                             14.939365  2.758280 395  5.416188  0.0000
# speciesCylicocyclus insigne                               -1.121834  2.758280 395 -0.406715  0.6844
# speciesCylicocyclus leptostomus                           11.201301  2.758280 395  4.060972  0.0001
# speciesCylicocyclus nassatus                               3.807850  2.758280 395  1.380516  0.1682
# speciesCylicostephanus calicatus                          -2.184587  2.758280 395 -0.792010  0.4288
# speciesCylicostephanus goldi                              -2.997239  2.758280 395 -1.086633  0.2779
# speciesCylicostephanus longibursatus                       2.214812  2.758280 395  0.802969  0.4225
# speciesCylicostephanus minutus                             5.759070  2.758280 395  2.087920  0.0374

ggplot(dfTraj,aes(x = Day, y = sqrt(count), col = Group, group = Id)) +
  geom_point() + geom_line() + facet_wrap(~ species, scales = 'free')

dfTraj=dfTraj[dfTraj$Group=="Chicory",]

ggplot(dfTraj,aes(x = Day, y = sqrt(count), col = Group, group = Id)) +
  geom_point() + geom_line() + facet_wrap(~ species, scales = 'free')

ggplot(dfTraj,aes(x = Day, y = sqrt(sqrt(count)), group = paste0(Day,Group), fill = Group)) +
  #geom_point(size = 3,alpha = .6,position = position_dodge(width = .5))  +
  geom_boxplot(alpha=0.4) +
  theme_classic() +
  theme(legend.position='bottom',strip.background = element_blank(), text = element_text(size = 12),
        strip.text = element_text(size = 12),
        axis.title = element_text(size = 14)) +
  ylab('Transformed counts') +
  scale_fill_manual(values = c('#006d2c')) +
  facet_wrap(~ species, scales='free')

dfTraj45=dfTraj[dfTraj$Day=="D45",]

dfcountd45=aggregate(count ~ species, 
                     FUN = sum, data = dfTraj45)
#                          species   count
# 1         Coronocyclus coronatus   53033
# 2          Coronocyclus labiatus     277
# 3         Cyathostomum catinatum   68538
# 4         Cyathostomum pateratum    1188
# 5         Cylicocyclus ashworthi     394
# 6           Cylicocyclus insigne       0
# 7       Cylicocyclus leptostomus     166
# 8          Cylicocyclus nassatus     541
# 9      Cylicostephanus calicatus      60
# 10         Cylicostephanus goldi      23
# 11 Cylicostephanus longibursatus  513455
# 12       Cylicostephanus minutus 1771732

dfcountd45 <- dfcountd45[order(-dfcountd45$count),]
# species   count
# 12       Cylicostephanus minutus 1771732
# 11 Cylicostephanus longibursatus  513455
# 3         Cyathostomum catinatum   68538
# 1         Coronocyclus coronatus   53033
# 4         Cyathostomum pateratum    1188
# 8          Cylicocyclus nassatus     541
# 5         Cylicocyclus ashworthi     394
# 2          Coronocyclus labiatus     277
# 7       Cylicocyclus leptostomus     166
# 9      Cylicostephanus calicatus      60
# 10         Cylicostephanus goldi      23
# 6           Cylicocyclus insigne       0


### Alpha chicory to be tested within each day
alpha_div <- estimate_richness(ci, split = TRUE, measure = c("Shannon", "Simpson"))
alpha_div$sample.id <- rownames(alpha_div) %>%  as.factor()

alphaplot <- sample_data(ci) %>%
  unclass() %>%
  data.frame() %>%
  left_join(alpha_div, by = "sample.id") %>%
  reshape2::melt(measure.vars = c("Shannon","Simpson"),
                 variable.name = "diversity_measure",
                 value.name = "alpha_diversity")

alphaplot$Group<- factor(alphaplot$Group, levels = c('Control','Chicory'))

Shannon=ggplot(alphaplot) +
  geom_boxplot(aes(x = diversity_measure, y = alpha_diversity, fill = Group), alpha = .4)+
  facet_wrap(~ Day) +
  scale_y_continuous(limits = c(0, 2.2), breaks = seq(0, 3, .2)) +
  labs(x = "Group", y = "Alpha diversity", color = "Pipeline") +
  theme(legend.position = 'bottom', text = element_text(size = 16))+
  theme(axis.line = element_line(size = 1, linetype = "solid"),
        legend.title = element_blank(),legend.position="bottom",
        legend.text = element_text(size=33, family = "Lato"),
        axis.title.x = element_text(size=35, family = "Lato", margin = margin(t = 0.4, unit="cm")), 
        axis.title.y = element_text(size=35, family = "Lato", margin = margin(r = 0.4, unit="cm")),
        axis.text.x = element_text(size=25, face = "bold", family = "Lato"),
        axis.text.y = element_text(size=25, face='bold', family = "Lato"))+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))+
  theme(strip.background = element_rect(fill="white", size=1.5, linetype="solid"),
        strip.text.x = element_text(size=20, family = "Lato", face="bold"))+
  scale_fill_manual(values=c('#737373','#006d2c'))
Shannon

ragg::agg_tiff("~/save/Chicory/Chicory paper/Alpha_plot.tiff", width = 13, height = 7, units = "in", res = 300)
Shannon
dev.off()

Shan=alphaplot[alphaplot$diversity_measure=="Shannon",]

t.test(alpha_diversity ~ Group, 
       data = Shan[Shan$Day=='d16',])
# Welch Two Sample t-test
# 
# data:  alpha_diversity by Group
# t = 1.4112, df = 7.8517, p-value = 0.1966
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.2249860  0.9286215
# sample estimates:
#   mean in group Control mean in group Chicory 
# 1.3012671             0.9494493 

t.test(alpha_diversity ~ Group, 
       data = Shan[Shan$Day=='d31',])
# Welch Two Sample t-test
# 
# data:  alpha_diversity by Group
# t = 3.4912, df = 3.6498, p-value = 0.0291
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   0.1735565 1.8246722
# sample estimates:
#   mean in group Control mean in group Chicory 
# 1.8120491             0.8129347

t.test(alpha_diversity ~ Group, 
       data = Shan[Shan$Day=='d45',])
# Welch Two Sample t-test
# 
# data:  alpha_diversity by Group
# t = 4.1705, df = 3.4397, p-value = 0.01912
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   0.3123737 1.8476638
# sample estimates:
#   mean in group Control mean in group Chicory 
# 1.6307025             0.5506838 

Simps=alphaplot[alphaplot$diversity_measure=="Simpson",]

t.test(alpha_diversity ~ Group, 
       data = Simps[Simps$Day=='d16',])
# Welch Two Sample t-test
# 
# data:  alpha_diversity by Group
# t = 1.259, df = 5.8166, p-value = 0.2562
# alternative hypothesis: true difference in means between group Control and group Chicory is not equal to 0
# 95 percent confidence interval:
#   -0.1420475  0.4384355
# sample estimates:
#   mean in group Control mean in group Chicory 
# 0.6323882             0.4841942 

t.test(alpha_diversity ~ Group, 
       data = Simps[Simps$Day=='d31',])
# Welch Two Sample t-test
# 
# data:  alpha_diversity by Group
# t = 2.5655, df = 3.1745, p-value = 0.07826
# alternative hypothesis: true difference in means between group Control and group Chicory is not equal to 0
# 95 percent confidence interval:
#   -0.07375822  0.80118306
# sample estimates:
#   mean in group Control mean in group Chicory 
# 0.7812203             0.4175078 

t.test(alpha_diversity ~ Group, 
       data = Simps[Simps$Day=='d45',])
# Welch Two Sample t-test
# 
# data:  alpha_diversity by Group
# t = 2.9516, df = 3.1723, p-value = 0.05593
# alternative hypothesis: true difference in means between group Control and group Chicory is not equal to 0
# 95 percent confidence interval:
#   -0.01997705  0.89153890
# sample estimates:
#   mean in group Control mean in group Chicory 
# 0.7501682             0.3143872 

### NMDS chicory
pslog <- transform_sample_counts(ci, function(x) log(1 + x))
psbin  = transform_sample_counts(ci, function(x) ifelse(x>0,1,0))

######------- PCoA Bray
out.bra.log <- ordinate(pslog, method = "PCoA", distance = "bray")
evals <- out.bra.log$values$Eigenvalues

## Method
p = plot_ordination(pslog, out.bra.log, 
                    color = "Group", shape  = "Day") +
  #labs(col = "Pipeline") + 
  geom_point(size = 6, alpha = .4) +
  #coord_fixed(sqrt(evals[2] / evals[1])) + 
  ggtitle(paste0(''))+
  theme(legend.position = 'bottom', text = element_text(size = 16),legend.direction = "vertical")+
  guides(color = guide_legend(reverse=TRUE))+
  theme(axis.line = element_line(size = 1, linetype = "solid"),
        legend.title = element_blank(),
        legend.text = element_text(size=33, family = "Lato"),
        axis.title.x = element_text(size=35, family = "Lato", margin = margin(t = 0.4, unit="cm")), 
        axis.title.y = element_text(size=35, family = "Lato", margin = margin(r = 0.4, unit="cm")),
        axis.text.x = element_text(size=25, face = "bold", family = "Lato"),
        axis.text.y = element_text(size=25, face='bold', family = "Lato"))+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))+
  scale_color_manual(values=c('#737373','#006d2c'))
p


Shannon_Bray<-(Shannon/p)+
  plot_layout(guides = "collect")
# theme(legend.position = "right")
Shannon_Bray

########## In vitro anthelmintic activity evaluation of the SLs extract ##########
#=== IC50 mesurement (LDA) ====
win = read.csv(file='LDA_in_vitro_SLs_XP_Chicory.csv',header=T,sep=';',dec=',',fileEncoding="latin1")
win=win[win$Conc!="2500" & win$Conc!="3500" ,]
win$id=row.names(win)

## Express LD as a percentage of CTL
ctl = mean(win$L3[win$Conc==0 ]/(win$L3[win$Conc==0]+win$L1.L2[win$Conc==0]))
win$pld = (win$L3/(win$L1.L2+win$L3))/ctl

###Outliers
# OUT=win %>%
#   group_by(Conc) %>%
#   identify_outliers(pld)
# dtout=data.frame(Date=OUT$Date,Assay=OUT$Assay, Group=OUT$Group, Conc=OUT$Conc,
#                  Unit=OUT$Unit, L1.L2=OUT$L1.L2, L3=OUT$L3, Dev=OUT$Dev,id=OUT$id, pld=OUT$pld)
# win=anti_join(win,dtout, by = "id")

win_N=win[win$Assay=="Nouzilly",]
win_C=win[win$Assay=="Chamberet",]
## Model
#-- 1st model : different slopes, different ED betw plates 
mod = drm(pld ~ Conc,Group, data = win, fct = LL.2())
summary(mod)
plot(mod,type="all",col=3,lwd=2)
win$Assay=factor(win$Assay)

mod_N = drm(pld ~ Conc, data = win_N, fct = LL.2())
summary(mod_N)
plot(mod_N,type="all",col=3,lwd=2)
win_N$Assay=factor(win_N$Assay)

mod_C= drm(pld ~ Conc, Assay, data = win_C, fct = LL.2())
summary(mod_C)
plot(mod_C,type="all",col=3,lwd=2)
win_C$Assay=factor(win_C$Assay)

# predictions and confidence intervals.
demo.fits_N <- expand.grid(Conc = exp(seq(log(1.00e-02), log(1000000), length=1000)))
demo.fits_N$Group=rep(levels(win_N$Assay))

demo.fits_C <- expand.grid(Conc = exp(seq(log(1.00e-02), log(1000000), length=1000)))
demo.fits_C$Group=rep(levels(win_C$Assay))

# demo.fits=rbind(demo.fits, demo.fits1)

# new data with predictions
#Nouzilly
pm <- predict(mod_N, newdata=demo.fits_N, interval="confidence") 
demo.fits_N$p <- pm[,1]
demo.fits_N$pmin <- pm[,2]
demo.fits_N$pmax <- pm[,3]

win_N$XX = win_N$Conc
win_N$XX[win_N$XX == 0] = 0.01
win_N$pt = factor(paste0(win_N$Group,'_',win_N$Conc))
win_N$XX = factor(win_N$XX)

win_N$Conc[win_N$Conc==0]= 0.01 
win_N$facet <- ifelse(win_N$Conc == min(win_N$Conc), 1, 2)
demo.fits_N$facet <- ifelse(demo.fits_N$Conc == min(demo.fits_N$Conc), 1, 2)
demo.fits_N = subset(demo.fits_N, Conc>0.01) 

df_N = data.frame(win_N %>% 
                    group_by(Conc, pt,Group) %>% 
                    summarize(avg = mean(pld), n=n(), sd =sd(pld), se=sd/sqrt(n)))
df_N$facet <- ifelse(df_N$Conc == min(df_N$Conc), 1, 2)
df_N$Conc[df_N$Conc==min(df_N$Conc)]=0.01 
win_N$Conc=factor(win_N$Conc)

#Chamberet
pm <- predict(mod_C, newdata=demo.fits_C, interval="confidence") 
demo.fits_C$p <- pm[,1]
demo.fits_C$pmin <- pm[,2]
demo.fits_C$pmax <- pm[,3]

win_C$XX = win_C$Conc
win_C$XX[win_C$XX == 0] = 0.01
win_C$pt = factor(paste0(win_C$Group,'_',win_C$Conc))
win_C$XX = factor(win_C$XX)

win_C$Conc[win_C$Conc==0]= 0.01 
win_C$facet <- ifelse(win_C$Conc == min(win_C$Conc), 1, 2)
demo.fits_C$facet <- ifelse(demo.fits_C$Conc == min(demo.fits_C$Conc), 1, 2)
demo.fits_C = subset(demo.fits_C, Conc>0.01) 

df_C = data.frame(win_C %>% 
                    group_by(Conc, pt,Group) %>% 
                    summarize(avg = mean(pld), n=n(), sd =sd(pld), se=sd/sqrt(n)))
df_C$facet <- ifelse(df_C$Conc == min(df_C$Conc), 1, 2)
df_C$Conc[df_C$Conc==min(df_C$Conc)]=0.01 
win_C$Conc=factor(win_C$Conc)
demo.fits_N$Group="INRAE"
df_N$Group="INRAE"

#Graph
Plot_IC50_Chicory=ggplot() +
  geom_ribbon(data = demo.fits_N,aes(x=Conc, y=p, ymin=pmin, ymax=pmax, group=Group, fill = Group), alpha=0.2)+
  geom_line(data = demo.fits_N, aes(x=Conc, y=p, group=Group))+
  geom_ribbon(data = demo.fits_C,aes(x=Conc, y=p, ymin=pmin, ymax=pmax, group=Group, fill = Group), alpha=0.2)+
  geom_line(data = demo.fits_C, aes(x=Conc, y=p, group=Group))+
  annotation_logticks(scaled = TRUE,sides="b")+
  scale_x_log10() +
  geom_point(data = df_N,aes(x = Conc, y = avg, shape=Group, size=Group, fill=Group))+
  geom_errorbar(data = df_N, aes(x = Conc, ymin = df_N$avg-df_N$se,ymax = df_N$avg + df_N$se,width = 0.02))+
  geom_point(data = df_C,aes(x = Conc, y = avg, shape=Group, size=Group, fill=Group))+
  geom_errorbar(data = df_C, aes(x = Conc, ymin = df_C$avg-df_C$se,ymax = df_C$avg + df_C$se,width = 0.02))+
  scale_shape_manual(values=c(20,15 ,20))+
  scale_size_manual(values=c(2,2,2))+
  theme_classic() +
  scale_x_log10(breaks = c(0.01,0.1,1, 10,100, 1000, 10000,100000,1000000),  limits = c(0.01,1000000),
                labels = c("0","0.1","1", "10","100", "1000", "10000","100000","1000000"))+
  labs(title=, y='Development percentage \n to the control', x='Concentrations (µg/mL)')+
  theme(axis.line = element_line(size = 1, linetype = "solid"),
        legend.title = element_blank(),legend.position="none",
        legend.text = element_text(size=33, family = "Lato"),
        axis.title.x = element_text(size=35, family = "Lato", margin = margin(t = 0.4, unit="cm")), 
        axis.title.y = element_text(size=35, family = "Lato", margin = margin(r = 0.4, unit="cm")),
        axis.text.x = element_text(size=25, face = "bold", family = "Lato"),
        axis.text.y = element_text(size=25, face='bold', family = "Lato"))+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))+
  guides(fill = guide_legend(nrow = 1))+
  # scale_fill_discrete(labels = c("Chamberet", "Control", "INRAE"))+ 
  scale_fill_manual(values=c("#a50f15","#969696","#08306b"))+
  scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1))
Plot_IC50_Chicory

## IC50 values
CI50<-ED(mod,50,interva="delta")
# Estimated effective doses
# 
#                Estimate Std. Error    Lower    Upper
# e:Chamberet:50 1.000813   0.183122 0.637117 1.364510
# e:Nouzilly:50  3.429450   0.098921 3.232985 3.625915

## Comparaison between compounds
compParm(mod,"e","-")
#                    Estimate Std. Error t-value   p-value    
# Nouzilly-Chamberet  2434.44     200.06  12.169 < 2.2e-16 ***

### === FECRT against pyrantel ====
FECr = read.csv(file='FECRT_Pyrantel_INRAE_&_Chamberet_XP_Chicory.csv',header=T,sep=';',dec=',',fileEncoding="latin1")

I_fec0.trt = FECr$FEC[FECr$Day =='d0' & FECr$Treatment=='HIGH_TRT']
I_fec15.trt = FECr$FEC[FECr$Day =='d15' & FECr$Treatment=='HIGH_TRT']
I_fec0.ctl = FECr$FEC[FECr$Day =='d0' & FECr$Treatment=='HIGH_CTL']
I_fec15.ctl = FECr$FEC[FECr$Day =='d15' & FECr$Treatment=='HIGH_CTL']
C_fec0_trt= FECr$FEC[FECr$Day =='d0' & FECr$Treatment=='TRT']
C_fec14_trt = FECr$FEC[FECr$Day =='d14' & FECr$Treatment=='TRT']
C_fec0_ctl= FECr$FEC[FECr$Day =='d0' & FECr$Treatment=='CTL']
C_fec14_ctl = FECr$FEC[FECr$Day =='d14' & FECr$Treatment=='CTL']

###----- FECr day15
FECRT_I0 <- eggCounts::fecr_stan(I_fec0.ctl, I_fec0.trt, rawCounts = TRUE, preCF=50, postCF=50,
                                paired = TRUE, indEfficacy = TRUE)
#                         mean        sd      2.5%        50%      97.5%  HPDLow95      mode  HPDHigh95
# FECR                  0.6143    0.1556    0.3284     0.6138     0.9147    0.3330    0.6613     0.9167
# meanEPG.untreated 11361.2914 3182.2345 6050.8741 11114.7635 18259.8758 5660.7082 9675.7541 17727.6197
# meanEPG.treated    4396.1355 2240.6386  848.8773  4077.6448  9410.4345  388.6232 3400.0162  8709.2150

FECRT_I15 <- eggCounts::fecr_stan(I_fec15.ctl, I_fec15.trt, rawCounts = TRUE, preCF=50, postCF=50,
                                 paired = TRUE, indEfficacy = TRUE)
#                         mean        sd      2.5%        50%      97.5%  HPDLow95       mode HPDHigh95
# FECR                  0.9152    0.0684    0.7345     0.9312     0.9956    0.7811     0.9667     1.000
# meanEPG.untreated 12489.6926 2913.7049 7385.6274 12366.9004 18690.4976 6919.9703 12008.1379 18072.498
# meanEPG.treated    1056.2624  896.9010   50.1153   822.4317  3349.2235    0.4151   386.8969  2829.083


FECRT_C0 <-  eggCounts::fecr_stan(C_fec0_ctl, C_fec0_trt, rawCounts = TRUE, preCF=50, postCF=50,
                                  paired = TRUE, indEfficacy = TRUE)
#                         mean        sd      2.5%       50%      97.5%  HPDLow95       mode  HPDHigh95
# FECR                  0.2816    0.0959    0.1461     0.263     0.5163    0.1261     0.2271     0.4729
# meanEPG.untreated 12660.4628 3672.6803 6596.3783 12374.638 20655.2063 6203.2847 11824.9979 20061.9716
# meanEPG.treated    9080.1705 2869.7094 4355.1869  8868.221 15390.9887 3325.7777  8718.0609 14305.3592

FECRT_C14 <- eggCounts::fecr_stan(C_fec14_ctl, C_fec14_trt, rawCounts = TRUE, preCF=50, postCF=50,
                                  paired = TRUE, indEfficacy = TRUE)
#                         mean        sd      2.5%        50%      97.5%  HPDLow95       mode  HPDHigh95
# FECR                  0.9927    0.0073    0.9750     0.9943     0.9995    0.9818     0.9955     1.0000
# meanEPG.untreated 12084.8031 3354.1243 6445.3670 11658.7404 19577.4183 6069.5224 10591.5410 18834.0604
# meanEPG.treated      88.4058   92.5119    5.3286    67.3077   309.5048    0.0579    41.7179   228.2737
####################################################################################################################