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
library("knitr")
library("BiocStyle")
require('dada2')
require('phyloseq')
library(data.table)
library(scales)
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
path<-"~/work/R/Chicory"
setwd(path)

########## Herbage characteristics and dietary choices ##########
# === Body weight data ====
BW=read.csv(file='Weight_XP_Chicory.csv',header=T,sep=';',dec=',')
BW$Group<- factor(BW$Group, levels = c('Control','Chicory'))

BW_summary=data_summary(BW, varname="Weight",groupnames=c("Day",'Group'))
BW_summary
#   Day   Group  Weight       sd
# 1  d0 Control 465.150 16.01055
# 2  d0 Chicory 466.648 16.14718
# 3 d45 Control 454.708 19.43229
# 4 d45 Chicory 470.103 19.65939

# === Statistical analysis ====
bwC=BW[BW$Group=='Chicory',]
bwT=BW[BW$Group=='Control',]

D45chico=bwC[bwC$Day=="d45",]
D45control=bwT[bwT$Day=="d45",]
D45chico=D45chico$Weight
D45control=D45control$Weight

wilcox.test(D45control,D45chico)
# data:  D45control and D45chico
# W = 31, p-value = 0.1655
# alternative hypothesis: true location shift is not equal to 0

########## Chicory effect on FEC ##########
# === Parasite data ====
Chico = read.csv(file='Parasite_data_XP_Chicory.csv',header=T,sep=';',dec=',')
Chico$Group<- factor(Chico$Group, levels = c('Control','Chicory'))

FEC_summary=data_summary(Chico, varname="EPG",groupnames=c("Day",'Group'))
FEC_summary
# Day  Group    EPG       sd
# D0 Control 2169.0 853.5573
# D0 Chicory 2169.0 794.7390
# D16 Control 1084.5 383.0176
# D16 Chicory  285.0 175.7840
# D31 Control 1233.0 654.8206
# D31 Chicory  282.0 307.9610
# D45 Control  864.0 316.9543
# D45 Chicory  153.0 215.3834

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
D0Con=Chico[Chico$Day=='d0' & Chico$Group=="Control",]
D16Con=Chico[Chico$Day=='d16' & Chico$Group=="Control",]
x=D0Con$EPG
y=D16Con$EPG
t.test(x,y)
# data:  x and y
# t = 3.6657, df = 12.483, p-value = 0.003037
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   442.6585 1726.3415
# sample estimates:
#   mean of x mean of y 
# 2169.0    1084.5 

D0Chico=Chico[Chico$Day=='d0' & Chico$Group=="Chicory",]
D16Chico=Chico[Chico$Day=='d16' & Chico$Group=="Chicory",]
D0Chico=D0Chico$EPG
D16Chico=D16Chico$EPG
t.test(D0Chico,D16Chico)
# data:  x and y
# t = 3.6657, df = 12.483, p-value = 0.003037
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   442.6585 1726.3415
# sample estimates:
#   mean of x mean of y 
# 2169.0    1084.5 
D0Chico=Chico[Chico$Day=='d0' & Chico$Group=="Chicory",]
D16Chico=Chico[Chico$Day=='d16' & Chico$Group=="Chicory",]
D0Chico=D0Chico$EPG
D16Chico=D16Chico$EPG
t.test(D0Chico,D16Chico)
# data:  D0Chico and D16Chico
# t = 7.3196, df = 9.8785, p-value = 2.717e-05
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   1309.536 2458.464
# sample estimates:
#   mean of x mean of y 
# 2169       285 

gee.fit.EPG <- geeglm(EPG ~ Group * Day,id = Horses, data = Chico, family = poisson,
                      corstr = "ar1", scale.fix = TRUE, std.err = "san.se")
summary(gee.fit.EPG)

# 
# Coefficients:
#                         Estimate   Std.err    Wald Pr(>|W|)    
# (Intercept)          7.68e+00  1.18e-01 4234.11  < 2e-16 ***
# GroupChicory        -5.77e-18  1.61e-01    0.00  1.00000    
# DayD16              -6.93e-01  1.59e-01   19.09  1.2e-05 ***
# DayD31              -5.65e-01  1.98e-01    8.11  0.00439 ** 
# DayD45              -9.20e-01  1.61e-01   32.52  1.2e-08 ***
# GroupChicory:DayD16 -1.34e+00  2.67e-01   24.98  5.8e-07 ***
# GroupChicory:DayD31 -1.48e+00  3.98e-01   13.71  0.00021 ***
# GroupChicory:DayD45 -1.73e+00  4.65e-01   13.84  0.00020 ***
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#==== Chicory efficacy ====
#data EPG chicory and control for d0 and d45
d=Chico[Chico$Day=="d0"|Chico$Day=="d45",]
d=d[,c('Horses', 'Group','Day','EPG')]
colnames(d)=c('ID','Group','Day','EPG')
d$Block="Chico"

##-- sample control
tm0 = sample(d$EPG[d$Group=='Control' & d$Day=='d0'],size = 10)
tmend = sample(d$EPG[d$Group=='Control' & d$Day=='d45'],size = 10)

##-- sample treated
trt0 = sample(d$EPG[d$Group=='Chicory' & d$Day=='d0'],size = 10)
trtend = sample(d$EPG[d$Group=='Chicory' & d$Day=='d45'],size = 10)

###----- Bayesian hierarchical models
FECRT_Chico <- eggCounts::fecr_stan(trt0, trtend, rawCounts = TRUE,
                                    paired = TRUE, indEfficacy = TRUE)
FECRT_Chico$posterior.summary
#                         mean        sd      2.5%        50%      97.5%  HPDLow95       mode  HPDHigh95
# FECR                  0.9462    0.0341    0.8567     0.9541     0.9879    0.8793     0.9591     0.9923
# meanEPG.untreated 12836.5068 3681.7347 6554.6347 12542.2863 20789.5270 6180.2051 12423.2995 19905.5972
# meanEPG.treated     691.9600  506.9108  134.5365   571.1150  2017.3476   59.7573   439.5252  1660.9121

###----- Bootstrap aproach 
bdatok=d
vecfarm = unique(bdatok$Block)
cicross = NULL
n = 0

for(f in vecfarm){
  ##-- subset farm of interest
  d = bdatok[bdatok$Block==f,]
  ##-- iteration
  n = n + 1
  Eff = array(NA,1000)
  
  for(i in 1:1000){
    ##-- sample control
    tm0 = sample(d$EPG[d$Group=='Control' & d$Day=='d0'],size = 10,replace = T)
    tmend = sample(d$EPG[d$Group=='Control' & d$Day=='d45'],size = 10,replace = T)
    
    ##-- sample treated
    trt0 = sample(d$EPG[d$Group=='Chicory' & d$Day=='d0'],size = 10,replace = T)
    trtend = sample(d$EPG[d$Group=='Chicory' & d$Day=='d45'],size = 10,replace = T)
    
    ##-- FECR
    tm.0=mean(tm0)
    tm.end=mean(tmend)
    trt.0=mean(trt0)
    trt.end=mean(trtend)
    trt = trt.end/trt.0
    tm = tm.end/tm.0
    Eff[i] = ((1-((trt)/(tm)))*100)
    if(Eff[i]<0){Eff[i]=0} 
  }
  #b0 = fecrtCI(d$PEQ0,d$PEQ14,paired=TRUE,R=1000,alpha=.05)
  cicross$Lot[n] = f
  tm=mean(bdatok$EPG[bdatok$Block==f & bdatok$Day=='d45'& bdatok$Group=='Control'])/mean(bdatok$EPG[bdatok$Block==f & bdatok$Day=='d0'& bdatok$Group=='Control'])
  trt=mean(bdatok$EPG[bdatok$Block==f & bdatok$Day=='d45'& bdatok$Group=='Chicory'])/mean(bdatok$EPG[bdatok$Block==f & bdatok$Day=='d0'& bdatok$Group=='Chicory'])
  cicross$Efficacy[n] =((1-((trt)/(tm)))*100) 
  cicross$CIdw[n] = quantile(Eff, 0.025, na.rm=T)
  cicross$CIup[n] = quantile(Eff, 0.975, na.rm=T)
  rm(Eff)
  #  rm(b0)
}
cicross=data.frame(cicross)
cicross
#     Lot Efficacy     CIdw     CIup
# 1 Chico 82.29167 62.63515 96.0607

########## Chicory effect on larval development ##########
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
  labs(title=, y='Developpement percentage (%)', x='Day')+
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

########## Chicory effect on equine cyathostomins larval community structure ##########
path<-"~/work/R/Chicory/Nem"
path2 <- "~/work/R/Chicory/Nem/chicory_dada2_R_ci_mxee25_trunc200_BS16"
setwd(path2)

tr=read.table("track.tsv")
head(tr)
summary(tr$filtered)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 820  402444  458629  465366  548932  703166 

setwd(path)
###----====== 
set.seed(100) # Initialize random number generator for reproducibility

### Training set - once and for all
taxmeth = 'idtaxa'
train <- readDNAStringSet("~/work/R/db/curated_nemaca_ITS2_Strongylidae.fasta") 
tax <- read_tsv("~/work/R/db/idtaxa_03022022.tax") 

trainingSet <- LearnTaxa(train, names(train), tax)

####-------- Create ITS phyloseq object
truncL = 200
mxee = 25
bc = 'ITS'
BS = 16

### Read in seqtab
seqtab_nochim = read.table(file = paste0('~/work/R/Chicory/Nem/chicory_dada2_R_ci_mxee',mxee,'_trunc',truncL,'_BS',BS,'/output.tsv'),
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
#10

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

### Remove d0 and keep horses with 3 occurrences
psci.f = subset_samples(psITS.f.wk,Day !='d0')
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
ct = subset_samples(ct,Day !='d0' & Day!="d31" &Group!='Mo')

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
  theme(strip.background = element_rect(fill="white", size=1.5, linetype="solid"),
        strip.text.x = element_text(size=20, family = "Lato", face="bold"))+
  scale_fill_manual(values=c('#737373','#006d2c'))
Shannon

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
win = read.csv(file='Test_in_vitro_SL.csv',header=T,sep=';',dec=',',fileEncoding="latin1")
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
#                Estimate Std. Error t-value   p-value    
# Nouzilly-Chamberet  2434.44     200.06  12.169 < 2.2e-16 ***

### === FECRT against pyrantel ====
FECr = read.csv(file='FECRT_Nouzilly_Chamberet.csv',header=T,sep=';',dec=',',fileEncoding="latin1")
FECr = FECr[FECr$day %in% c('D0','D14','D15'),]

N_fec0.trt = FECr$FEC[FECr$day =='D0' & FECr$Treatment=='HIGH_TRT']
N_fec15.trt = FECr$FEC[FECr$day =='D15' & FECr$Treatment=='HIGH_TRT']
C_fec0_trt= FECr$FEC[FECr$day =='D0' & FECr$Treatment=='TRT']
C_fec14_trt = FECr$FEC[FECr$day =='D14' & FECr$Treatment=='TRT']

###----- FECr day15
FECRT_N <- eggCounts::fecr_stan(N_fec0.trt, N_fec15.trt, rawCounts = TRUE,
                                paired = TRUE, indEfficacy = TRUE)
fecR_N=FECRT_N$posterior.summary
fecR_N=fecR_N[1,c(1:8)]

FECRT_C <- eggCounts::fecr_stan(C_fec0_trt, C_fec14_trt, rawCounts = TRUE,
                                paired = TRUE, indEfficacy = TRUE)
fecR_C=FECRT_C$posterior.summary
fecR_C=fecR_C[1,c(1:8)]

fecR_N$Group="Nouzilly"
fecR_C$Group="Chamberet"

FECRT=bind_rows(fecR_N,fecR_C)

FECRT$mean_100=FECRT$mean*100
FECRT$HPDLow95_100=FECRT$HPDLow95*100
FECRT$HPDHigh95_100=FECRT$HPDHigh95*100
FECRT
#            mean     sd   2.5%    50%  97.5% HPDLow95   mode HPDHigh95     Group mean_100  HPDLow95_100 HPDHigh95_100
# FECR...1 0.9105 0.0567 0.7662 0.9208 0.9855   0.8105 0.9267         1  Nouzilly    91.05         81.05           100
# FECR...2 0.9927 0.0094 0.9723 0.9950 0.9998   0.9800 0.9965         1 Chamberet    99.27         98.00           100