setwd('C:/Users/cmatt/Documents/sarscov2/fontenele/tsvsanalysis/env3/results')
library("dplyr")
library("tidyr")
library("ggplot2")
library("lattice")
library("data.table")
library("tidyverse")


#read in data

c<- read.csv("all.coverage.csv", header=TRUE) #, colClasses = c("date"="character"))
c <- c %>% filter(!grepl("Undetermined", ID))
v <- read.csv("ww.allvariants.csv")
v <- v %>% filter(!grepl("Undetermined", ID)) %>% subset(select=-c(GFF_FEATURE)) %>% distinct()
d <- read.csv("ww.alldepth.csv")
d <- d %>% filter(!grepl("Undetermined", ID))
c$batch <- as.integer(c$batch)
d$batch <- as.integer(d$batch)
v$POS <- as.integer(v$POS)

#fix dates: 
#turn dates into R type Dates
#some dates missing leading 0s, fix by adding 0 to date strings less than 6 digits
v$date[nchar(v$date)<6] <- paste("0", as.character(v$date[nchar(v$date)<6]), sep="")
v$plot_date <- as.Date(paste("2", substring(as.character(as.Date(v$date, '%m%d%Y')), 2), sep=""))

vc <- left_join(v, c, by=c("ID", "batch", "sequencer"))

#all SNVs & coverage from samtools

#get only SNVs significant by ivar fishers exact 
vsig <- vc %>% filter(sequencer == "NextSeq") %>%  filter(PVAL <= 0.05) 
#turn position col into integers (POS == position)
vsig$POS <- as.integer(vsig$POS)

#plot all snvs by sample ID\
#color by coverage?
ggplot(vsig, aes(x=POS, y=ID, color=coverage))+
  geom_point(size=1)+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  theme_minimal()


#plot all SNVS for high coverage samples
vhighcov <- vsig %>% filter(coverage >= 70)

ggplot(vhighcov, aes(x=POS, y=ID, color=coverage))+
  geom_point(size=1)+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  theme_minimal()

ggplot(vhighcov, aes(x=POS, y=ID, color=meandepth))+
  geom_point(size=1)+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  theme_minimal()



#plot depth and SNVs together

#add missing/depth 0 positions
d_nxsq <- d %>% filter(sequencer == "NextSeq")
max_pos <- d_nxsq %>% select(ID, POS) %>% group_by(ID) %>% summarise(mx=max(POS))
#min_pos <- d %>% select(ID, POS) %>% group_by(ID) %>% summarise(mn=min(POS))
missing <- data.frame()
#miss_tmp <- data.frame()
for (i in levels(as.factor(max_pos$ID))){
  #mid <- max_pos %>% filter(ID == i)
  #m <- max(mid$mx)
  empty<- data.frame(POS=1:29903, depth=0, ID = i)
  di <- d_nxsq %>% filter(ID == i)
  missing <- bind_rows(missing, anti_join(empty, di, by=c("ID", "POS"))) 
}


d_all <- data.frame()
d_all <- full_join(d_nxsq, missing, by=c("ID", "POS", "depth"))
d_all <- d_all %>% group_by(ID)
#sort by position, grouped by ID
d_all <- arrange(d_all, POS, .by_group=TRUE)

#d_all now contains the depth along the sars-cov2 genome at all positions for every sample

#get average of bins by position for d_all 
#use 'each' below to change window size. bin is a sequential digit, can increase rep if needed (not decr). dont mess w length.out
d_all <- d_all %>% group_by(ID) %>% mutate(bin = rep(1:350, each=200, length.out=length(ID)))
dall_avgs <- d_all %>% group_by(ID, bin) %>% summarise(dpavg = mean(depth))
d_all <- left_join(d_all, dall_avgs, by=c("ID", "bin"))
d_all <- d_all %>% group_by(ID, bin) %>% mutate(posavg = mean(POS))


#be careful when rerunning code for d_all. several steps use joins, and if you dont initialize the object you can do werid things and 
#paste in data you didnt mean to

dallplot <- d_all %>% select(ID, bin, dpavg, posavg) %>% distinct()


#add info for extraction method & influent vs sludge (location) to depth df. and dates
#two samples had no SNVs detected and thus were not included in vsig:
distinct(anti_join(select(d, ID), select(v, ID)))

methods <- v %>% select(ID, extr_method) %>% distinct()
methods<- bind_rows(methods, data.frame(ID=c("21_25NW_120320_Si3E1_S4", "CODWW_060121_KFE1_S34"), extr_method=c("Si", "KF")))
dallplot <- left_join(dallplot, methods, by=c("ID"))

wwsource <- v %>% select(ID, location) %>% distinct() 
wwsource <- bind_rows(wwsource, data.frame(ID=c("21_25NW_120320_Si3E1_S4", "CODWW_060121_KFE1_S34"), location=c("21_25NW", "CODWW")))
#change 'codwwtp' to 'codww'
wwsource["location"][wwsource["location"]=="CODWWTP"] <- "CODWW"
dallplot <- left_join(dallplot, wwsource, by=c("ID"))

wwdates <- v %>% select(ID, plot_date) %>% distinct()
add_dates<- data.frame(ID=c("21_25NW_120320_Si3E1_S4", "CODWW_060121_KFE1_S34"), plot_date=c('12032020', '06012020')) 
add_dates$plot_date <- as.Date(add_dates$plot_date, "%m%d%Y")
wwdates <- bind_rows(wwdates, add_dates)
dallplot <- left_join(dallplot, wwdates, by=c("ID"))



#plot depth from d_all and SNVs from vsig, ALL SAMPLES:
ggplot(data=filter(dallplot,dpavg > 0), aes(x=posavg, y=ID))+
  geom_point(aes(alpha=dpavg), stroke=0, shape=15, size=2.5, color= "cadetblue4")+
  scale_alpha_continuous()+
  geom_point(data=vsig, aes(x=POS, y=ID), size=0.5, color="darkorange")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  labs(title="Depth and SNVs for all samples", x="Position along the SARS-CoV2 Wuhan-1 reference genome",
       y="sample name", alpha = "depth")+
  theme_minimal()

#shading as color instead? doesnt work?
ggplot(data=filter(dallplot,dpavg > 0), aes(x=posavg, y=ID))+
  geom_point(aes(fill=dpavg), stroke=0, size=3)+ 
  scale_fill_continuous(high="blue", low="white")+#, color= "lightblue3")+
  geom_point(data=vsig, aes(x=POS, y=ID), size=1, color="darkorange")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  theme_minimal()


#compare extraction methods:
d_extr <- dallplot %>% group_by(location, plot_date) %>% filter(n()>599) %>% filter(extr_method != "KF") %>% ungroup()
vsig_extr <- vsig %>% filter(ID %in% d_extr$ID) 
#orig code 

#plot depth from d_all and SNVs from vsig
#facets = extraction meth
ggplot(data=filter(d_extr,dpavg > 0), aes(x=posavg, y=location))+
  geom_point(aes(alpha=dpavg), stroke=0, shape=15, size=4, color= "cadetblue4")+
  scale_alpha_continuous()+
  geom_point(data=vsig_extr, aes(x=POS, y=location), size=1, color="darkorange")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  facet_wrap(~extr_method, ncol=2)+
  labs(title="Depth and SNVs by extraction method", x="Position along the SARS-CoV2 Wuhan-1 reference genome",
       y="sampling location", alpha = "depth")+
  theme_minimal()

#facets=sample loc
ggplot(data=filter(d_extr,dpavg > 0), aes(x=posavg, y=extr_method))+
  geom_point(aes(alpha=dpavg), stroke=0, shape=15, size=2.5, color= "cadetblue4")+
  scale_alpha_continuous()+
  geom_point(data=vsig_extr, aes(x=POS, y=extr_method), size=0.5, color="darkorange")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  facet_wrap(~location, ncol=2)+
  theme_minimal()


#compare ww source methods:
d_codWW <- dallplot %>% group_by(plot_date) %>% filter(n()==300) %>% filter(grepl("COD", location)) %>% ungroup()
vsig_codWW <- vsig %>% filter(ID %in% d_codWW$ID) 
#orig code 

#plot depth from d_all and SNVs from vsig
#facets = ww vs sludge
ggplot(data=filter(d_codWW,dpavg > 0), aes(x=posavg, y=plot_date))+
  geom_point(aes(alpha=dpavg), stroke=0, shape=15, size=6, color= "cadetblue4")+
  scale_alpha_continuous()+
  geom_point(data=vsig_codWW, aes(x=POS, y=plot_date), size=1.5, color="darkorange")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  facet_wrap(~location, ncol=2)+
  labs(title="Depth and SNVs for City of Davis Waste Water Treatment Plant waste water vs. sludge", x="Position along the SARS-CoV2 Wuhan-1 reference genome",
       y="sampling date", alpha = "depth")+
  theme_minimal()

#facets=sample loc
ggplot(data=filter(d_codWW,dpavg > 0), aes(x=posavg, y=location))+
  geom_point(aes(alpha=dpavg), stroke=0, shape=15, size=2.5, color= "cadetblue4")+
  scale_alpha_continuous()+
  geom_point(data=vsig_codWW, aes(x=POS, y=location), size=0.5, color="darkorange")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  facet_wrap(~plot_date, ncol=2)+
  theme_minimal()







#--------Variants of concern---------------------------

#read in list of csv files created from NextStrain
filelist = list.files(pattern="*.ref.csv$")
vl <- lapply(filelist, read.csv, header=FALSE)
nn <- gsub(filelist, pattern="\\..*", replacement="")

#parse csvs
vl2 <- list()
for (i in 1:length(vl)){
  vl2[[i]] <- t(vl[[i]])
  vl2[[i]] <- as.data.frame(vl2[[i]]) %>% rename(mutation = V1) #%>% mutate(variant = rep(nn[i], length(vl2[[i]])))
  vl2[[i]]$variant <- nn[i]
  rownames(vl2[[i]]) <- NULL
  vl2[[i]]$mutation <- trimws(vl2[[i]]$mutation)
  vl2[[i]]$REF <- substr(vl2[[i]]$mutation, 1, 1) #(nchar(vl2[[i]]$mutation)-(nchar(vl2[[i]]$mutation) - 1)))
  n_last <- 1
  vl2[[i]]$ALT <- substr(vl2[[i]]$mutation, nchar(vl2[[i]]$mutation) - n_last + 1, nchar(vl2[[i]]$mutation))
  vl2[[i]]$POS <- as.integer(parse_number(vl2[[i]]$mutation))
}

#create dataframe with sets of SNVs characteristic of VOCs
vocs <- bind_rows(vl2)


voc_less <- vocs %>% select(variant, POS, REF, ALT)
hits <- data.frame()

for (i in levels(as.factor(voc_less$variant))){
  samples <- v %>% filter(PVAL<= 0.05) %>% filter(sequencer=="NextSeq")
  samples$POS <- as.integer(samples$POS)
  variant_snvs <- voc_less %>% filter(variant == i)
  temp <- inner_join(samples, variant_snvs, by=c("POS", "REF", "ALT"))
  hits <- bind_rows(hits, temp)
}


unique_vocs <- voc_less %>% group_by(POS, REF, ALT) %>% 
  filter(n()==1)%>%
  ungroup()

shared_vocs <- voc_less %>% group_by(POS, REF, ALT) %>% 
  filter(n()>1)%>%
  ungroup()


unique_hits <- data.frame()
for (i in levels(as.factor(voc_less$variant))){
  samples <- v %>% filter(PVAL<= 0.05) %>% filter(sequencer=="NextSeq")
  samples$POS <- as.integer(samples$POS)
  variant_snvs_un <- unique_vocs %>% filter(variant == i)
  temp <- inner_join(samples, variant_snvs_un, by=c("POS", "REF", "ALT"))
  unique_hits <- bind_rows(unique_hits, temp)
}

#extract month
unique_hits$month <- format(unique_hits$plot_date, "%m%y")
unique_hits$day <- format(unique_hits$plot_date, "%m/%d/%y")
unique_hits <- unique_hits[order(as.Date(unique_hits$plot_date, format="%Y-%m-%d")),]


shared_hits <- data.frame()
for (i in levels(as.factor(voc_less$variant))){
  samples <- v %>% filter(PVAL<= 0.05) %>% filter(sequencer=="NextSeq")
  samples$POS <- as.integer(samples$POS)
  variant_snvs_shared <- shared_vocs %>% filter(variant == i)
  temp <- inner_join(samples, variant_snvs_shared, by=c("POS", "REF", "ALT"))
  shared_hits <- bind_rows(shared_hits, temp)
}

#extract month
shared_hits$month <- format(shared_hits$plot_date, "%m")


hits %>% group_by(variant) %>% 
  ggplot(aes(x=POS, y=ID))+
  geom_point(size=3)+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  theme_minimal()




for (i in levels(as.factor(hits$variant))){
  p <- hits %>% filter(variant==i)
  print(ggplot(p, aes(x=POS, y=ID))+
          geom_point(size=3)+
          theme(axis.text.x=element_blank(),
                axis.ticks.x=element_blank())+
          theme_minimal()+
          ggtitle(i))
    
  Sys.sleep(2)
}




for (i in levels(as.factor(hits$variant))){
  s <- shared_hits %>% filter(variant==i)
  u <- unique_hits %>% filter(variant==i)
  varoc <- voc_less %>% filter(variant==i)
  print(ggplot(data=s, aes(x=POS, y=ID), color="blue")+
          geom_point(size=3)+
          geom_point(data=u, aes(x=POS, y=ID), color="green", size=3)+
          #geom_point(data=varoc, aes(x=POS, y=variant), color="red", size=3)+
        theme(axis.text.x=element_blank(),
              axis.ticks.x=element_blank())+
          theme_minimal()+
          ggtitle(i))
  
  Sys.sleep(2)
}


#unique_hits$plot_date <- as.Date(unique_hits$plot_date, "%Y-%m-%d")


#only sample with lots of delta hits
ggplot(data=filter(unique_hits, ID =="N12_066_072621_KFE1_S24"), aes(x=POS, y=ALT_FREQ, color=variant))+
  geom_bar(position="stack", stat = "identity")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  theme_minimal()+
  facet_wrap(~ID, ncol=1)

#by ID, ordered by date:
unique_hits$ID <- reorder(unique_hits$ID, desc(unique_hits$plot_date))
ggplot(data=unique_hits, aes(x=POS, y=ID, color=variant))+
  geom_point(aes(size=ALT_FREQ))+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  geom_jitter(width=0.01)+
  theme_minimal()

#by month, point size
unique_hits$month <- reorder(unique_hits$month, desc(unique_hits$plot_date))
ggplot(data=unique_hits, aes(x=POS, y=month, color=variant))+
  geom_point(aes(size=ALT_FREQ))+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  geom_jitter(width=0.01)+
  theme_minimal()

#function to reorder dates (this is the only way i could figure out to put newer dates at the bottom):
c_trans <- function(a, b, breaks = b$breaks, format = b$format) {
  a <- as.trans(a)
  b <- as.trans(b)
  
  name <- paste(a$name, b$name, sep = "-")
  
  trans <- function(x) a$trans(b$trans(x))
  inv <- function(x) b$inverse(a$inverse(x))
  
  trans_new(name, trans, inverse = inv, breaks = breaks, format=format)
  
}
rev_date <- c_trans("reverse", "date")

#FINAL PLOT, variants
#use scale_y_continuous with rev_date above to flip dates on y
#by day, point size, desc date
#unique_hits$day <- reorder(unique_hits$day, desc(unique_hits$plot_date))
ggplot(data=unique_hits, aes(x=POS, y=plot_date, color=variant))+
  geom_point(aes(size=ALT_FREQ))+
  scale_size(range=c(4,10))+
  geom_jitter(width=0.0)+
  scale_y_continuous(trans=rev_date, n.breaks=30)+
  theme_minimal()+
  labs(title="Alternate nucleotide frequency for SNVs unique to variants of concern over time", x="Position along the SARS-CoV2 Wuhan-1 reference genome",
       y="sampling date", size= "alternate frequency", alpha = "alternate nucleotide frequency")+
  theme(axis.text = element_text(size=20), axis.title = element_text(size=20), title=element_text(size=20))+
  theme(legend.key.size = unit(1, 'cm'),
        legend.key.height = unit(1, 'cm'))+
  theme(legend.title=element_text(size=20))+
  theme(legend.text = element_text(size=20))+
  guides(colour = guide_legend(override.aes = list(size=8)))



#does plotting smaller points last reveal more points? i think no 
small_last <- unique_hits[order(-as.numeric(factor(unique_hits$ALT_FREQ))),]

ggplot(data=small_last, aes(x=POS, y=plot_date, color=variant))+
  geom_point(aes(size=ALT_FREQ))+ #,colour="black", stroke=0.0001)+
  scale_size(range=c(4,10))+
  geom_jitter(width=0.0)+
  scale_y_continuous(trans=rev_date, n.breaks=30)+
  theme_minimal()+
  labs(title="Alternate nucleotide frequency for SNVs unique to variants of concern over time", x="Position along the SARS-CoV2 Wuhan-1 reference genome",
       y="sampling date", alpha = "alternate nucleotide frequency")+
  theme(axis.text = element_text(size=20), axis.title = element_text(size=20), title=element_text(size=20))+
  theme(legend.key.size = unit(1, 'cm'),
        legend.key.height = unit(1, 'cm'))+
  theme(legend.title=element_text(size=20))+
  theme(legend.text = element_text(size=20))+
  guides(colour = guide_legend(override.aes = list(size=8)))


#other versions of the plot above:

#with y axis by week
ggplot(data=unique_hits, aes(x=POS, y=desc(plot_date), color=variant))+
  geom_point(aes(size=ALT_FREQ))+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  geom_jitter(width=0.0)+
  scale_y_date(date_breaks="1 week")+
  theme_minimal()


ggplot(data=unique_hits, aes(x=POS, y=desc(plot_date), color=variant))+
  geom_point(aes(size=ALT_FREQ))+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  geom_jitter(width=0.0)+
  theme_minimal()



#by month, stacked bars:
ggplot(unique_hits, aes(x=POS, y=ALT_FREQ, color=variant))+
  geom_bar(position="stack", stat = "identity", size=1)+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  theme_minimal()+
  facet_wrap(~month, ncol=1)



#adding clinical data!
#NOTE: must filter for nsamp > 0 or ggplot will assign point to 0 values

clin <- read.csv("yolo_county.csv")
clin$week_start <- as.Date(clin$week_start, '%m/%d/%Y')
clin$week_end <- as.Date(clin$week_end, '%m/%d/%Y')
#long format
clin <- clin %>% pivot_longer(., cols=c("alpha", "beta", "delta", "gamma", "epsilon", "iota", "kappa", "mu"),
                              names_to="variant", values_to="nsamp")

#only dates before last ssample date
clinless <- clin %>% filter(week_end <= "2021-08-07")

#clinical VOCs in separate columns:
ggplot(data=filter(clinless, nsamp>0), aes(x=variant, y=week_start, color=variant))+
  geom_point(aes(size=nsamp))+
  scale_y_continuous(trans=rev_date, n.breaks=30)+
  scale_size(range=c(1,10))+
  theme_minimal()+
  labs(title="SARS_CoV2 positive samples collected in Yolo County: Count of Variants by week", x="Variant of Concern",
       y="week start", size = "number of positive samples")+
  theme(axis.text = element_text(size=20), axis.title = element_text(size=20), title=element_text(size=20))+
  theme(legend.key.size = unit(1, 'cm'),
        legend.key.height = unit(1, 'cm'))+
  theme(legend.title=element_text(size=20))+
  theme(legend.text = element_text(size=20))+
  guides(colour = guide_legend(override.aes = list(size=8)))


#with all variants in one column
#something is wrong. there are points being plotted that should not be, for ex the first row has large point for epsilon, but no such point exits. 
#these points arent being jittered for some reason. not present in the plot above
ggplot(data=filter(clinless, nsamp>0), aes(x=TRUE, y=week_start, color=variant))+
  geom_point(aes(size=nsamp))+
  scale_y_continuous(trans=rev_date, n.breaks=30)+
  theme_minimal()+
  geom_jitter(height=0, width=0.2)+
  scale_size(range=c(2,10))+
  labs(title="SARS_CoV2 positive samples collected in Yolo County: Count of Variants by week", x=" ",
       y="week start", size = "number of positive samples")+
  theme(axis.text = element_text(size=20), axis.title = element_text(size=20), title=element_text(size=20))+
  theme(legend.key.size = unit(1, 'cm'),
        legend.key.height = unit(1, 'cm'))+
  theme(legend.title=element_text(size=20))+
  theme(legend.text = element_text(size=20))+
  guides(colour = guide_legend(override.aes = list(size=8)))


#with count on the Y axis instead of represented by point size:
ggplot(data=filter(clinless, nsamp>0), aes(x=nsamp, y=week_start, color=variant))+
  geom_point(size=2)+
  scale_y_continuous(trans=rev_date, n.breaks=30)+
  #scale_size(range=c(2,10))+
  theme_minimal()

#same as above but w lines
ggplot(data=filter(clinless, nsamp>0), aes(x=week_start, y=nsamp, color=variant))+
  geom_point(size=3)+
  geom_line(aes(data=variant))+
  scale_x_continuous(trans=rev_date)+
  #scale_size(range=c(2,10))+
  theme_minimal()+
  coord_flip()



#delta

delta_u <- unique_hits %>% filter(variant=="delta")
delta_s <- shared_hits %>% filter(variant=="delta")

ggplot(data=filter(dallplot,dpavg > 0), aes(x=posavg, y=ID))+
  geom_point(aes(alpha=dpavg), stroke=0, shape=15, size=2.5, color= "cadetblue4")+
  geom_point(data=delta_s, aes(x=POS, y=ID), color="black", size = 1)+
  geom_point(data=delta_u, aes(x=POS, y=ID), color="darkorange", size=1)+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  theme_minimal()+
  ggtitle("delta")+
  labs(x="Position along the SARS-CoV2 Wuhan-1 reference genome",
       y="sample name", alpha = "depth")




ggplot(data=shared_hits, aes(x=POS, y=ID), color="blue")+
  geom_point(size=3)+
  geom_point(data=unique_hits, aes(x=POS, y=ID, color=variant), size=3)+
theme(axis.text.x=element_blank(),
      axis.ticks.x=element_blank())+
  theme_minimal()




