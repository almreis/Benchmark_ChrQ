setwd('/simon/ANDRE/Figures/github/Figure4')

library(data.table)
library(plotROC)
library(ggplot2)
library(MESS)


pdf('plots_figure3_v4.pdf',width=12,height=8,useDingbats = FALSE)

# SV calling analysis for short read libraries
# I can benchmark the SV calling at the breakend level or at the SV level
# At the breakend level all that matters is that the breakend be identified correctly within
# a specific window around the original breakpoint. Breakends window that window are classified
# True Positive (TP), outside the window as False Positive (FP) and missing breakends at expected breakpoints
# are classified as False Negative (FN).

# On the other to bechmark at the SV level, both breakpoints for a given SV need to be correctly identified
# for it to be classfied as True Positive (TP) and missing SVs are classified as False Negatives (FN).

# I'll first process the data at the Breakend level:

# Read the data for SV calling at different coverage levels

# At each coverage level I need to have:
# 1) all the breakend calls that were made. (*_call)
# 2) all the expected breakends. (*_base)
# 3) all the breakend calls that were made and overlap windows (20nt) around expected breakends. (*_overlap)

# I also did this analysis with both Lumpy and Manta and I think I need to keep them separate
# so that I can compare aples to aples.

# First I'll process the data for Lumpy:

lumpy_call <- list()
lumpy_base <- list()
lumpy_overlap <- list()
lumpy_vcf <- list()

all_cov <- list()

h <- 1
for (i in c('NA12878_SEQUIN_HISEQ_X_TEN_PCR_FREE_2','SEQUIN_HISEQ_2500_PCR_BASED_1','SEQUIN_HISEQ_2500_PCR_BASED_2','SEQUIN_NEXTSEQ_500_PCR_BASED_1','SEQUIN_NEXTSEQ_500_PCR_BASED_2','NA12878_SEQUIN_MGISEQ2000')) {
  for (j in as.character(1:20)) {
    dt1 <- fread(sprintf('curves/%s.sv.sorted.%s.lumpy.sorted.overlap.tab',i,j))
    dt1[,FILE:=sprintf('%s',i)]
    dt1[,COVERAGE:=j]
    lumpy_overlap[[h]] <- dt1
    
    dt2 <- fread(sprintf('curves/%s.sv.sorted.%s.lumpy.sorted.bed',i,j))
    dt2[,FILE:=sprintf('%s',i)]
    dt2[,COVERAGE:=j]
    lumpy_call[[h]] <- dt2
    
    dt3 <- fread('curves/sequin_ref_breakends.bed',header=FALSE)
    dt3[,FILE:=sprintf('%s',i)]
    dt3[,COVERAGE:=j] 
    lumpy_base[[h]] <- dt3
    
    dt4 <- fread(sprintf('zcat curves/%s.sv.sorted.%s.lumpy.sorted.vcf.gz',i,j),skip="#CHROM")
    dt4[,FILE:=sprintf('%s',i)]
    dt4[,COVERAGE:=j]
    lumpy_vcf[[h]] <- dt4
    
    dt5 <- fread(sprintf('curves/%s.sv.sorted.%s.summary.tab',i,j))
    dt5[,FILE:=sprintf('%s',i)]
    dt5[,COVERAGE:=j]
    all_cov[[h]] <- dt5
    
    h <- h+1
  }
}
lumpy_call <- rbindlist(lumpy_call)
lumpy_call[,SOFTWARE:="Lumpy"]

lumpy_base <- rbindlist(lumpy_base)
lumpy_base[,SOFTWARE:="Lumpy"]

lumpy_overlap <- rbindlist(lumpy_overlap)
lumpy_overlap[,SOFTWARE:="Lumpy"]

lumpy_vcf <- rbindlist(lumpy_vcf)

all_cov <- rbindlist(all_cov)
cov <- all_cov[,.(median(V2)),by=list(FILE,COVERAGE)]

# Now I'll process the data for Manta:

manta_call <- list()
manta_base <- list()
manta_overlap <- list()
manta_vcf <- list()

h <- 1
for (i in c('NA12878_SEQUIN_HISEQ_X_TEN_PCR_FREE_2','SEQUIN_HISEQ_2500_PCR_BASED_1','SEQUIN_HISEQ_2500_PCR_BASED_2','SEQUIN_NEXTSEQ_500_PCR_BASED_1','SEQUIN_NEXTSEQ_500_PCR_BASED_2','NA12878_SEQUIN_MGISEQ2000')) {
  for (j in as.character(1:20)) {
    dt1 <- fread(sprintf('curves/%s.sv.sorted.%s.manta.sorted.overlap.tab',i,j))
    dt1[,FILE:=sprintf('%s',i)]
    dt1[,COVERAGE:=j]
    manta_overlap[[h]] <- dt1
    
    dt2 <- fread(sprintf('curves/%s.sv.sorted.%s.manta.sorted.bed',i,j))
    dt2[,FILE:=sprintf('%s',i)]
    dt2[,COVERAGE:=j]
    manta_call[[h]] <- dt2
    
    dt3 <- fread('curves/sequin_ref_breakends.bed',header=FALSE)
    dt3[,FILE:=sprintf('%s',i)]
    dt3[,COVERAGE:=j] 
    manta_base[[h]] <- dt3
    
    dt4 <- fread(sprintf('zcat curves/%s.sv.sorted.%s.manta.sorted.vcf.gz',i,j),skip="#CHROM")
    dt4[,FILE:=sprintf('%s',i)]
    dt4[,COVERAGE:=j]
    manta_vcf[[h]] <- dt4
    
    h <- h+1
  }
}
manta_call <- rbindlist(manta_call)
manta_call[,SOFTWARE:="Manta"]

manta_base <- rbindlist(manta_base)
manta_base[,SOFTWARE:="Manta"]

manta_overlap <- rbindlist(manta_overlap)
manta_overlap[,SOFTWARE:="Manta"]

manta_vcf <- rbindlist(manta_vcf)

####################################################################

# Sniffles:

sniffles_call <- list()
sniffles_base <- list()
sniffles_overlap <- list()
sniffles_vcf <- list()

nano_cov <- list()

h <- 1
for (i in c('NA12878_SEQUIN_PROMETHION')) {
  for (j in as.character(1:20)) {
    info <- file.info(sprintf('curves/%s.sv.sorted.%s.sniffles.overlap.tab',i,j))
    if (info$size != 0) {
      dt1 <- fread(sprintf('curves/%s.sv.sorted.%s.sniffles.overlap.tab',i,j))
      dt1[,FILE:=sprintf('%s',i)]
      dt1[,COVERAGE:=j]
      sniffles_overlap[[h]] <- dt1
      
      dt2 <- fread(sprintf('curves/%s.sv.sorted.%s.sniffles.bed',i,j))
      dt2[,FILE:=sprintf('%s',i)]
      dt2[,COVERAGE:=j]
      sniffles_call[[h]] <- dt2
      
      dt3 <- fread('curves/sequin_ref_breakends.bed',header=FALSE)
      dt3[,FILE:=sprintf('%s',i)]
      dt3[,COVERAGE:=j]
      sniffles_base[[h]] <- dt3
      
      dt4 <- fread(sprintf('zcat curves/%s.sv.sorted.%s.sniffles.vcf.gz',i,j),skip="#CHROM")
      dt4[,FILE:=sprintf('%s',i)]
      dt4[,COVERAGE:=j]
      sniffles_vcf[[h]] <- dt4
      
      dt5 <- fread(sprintf('curves/%s.sv.sorted.%s.summary.tab',i,j))
      dt5[,FILE:=sprintf('%s',i)]
      dt5[,COVERAGE:=j]
      nano_cov[[h]] <- dt5
    }
      # dt1 <- data.table(V1=character(),V2=numeric(),V3=numeric(),V4=character(),
      #            V5=character(),V6=numeric(),V7=character(),V8=numeric(),
      #            V9=numeric(),V10=numeric(),V11=character(),V12=character(),
      #            V13=numeric(),FILE=character(),COVERAGE=character())
    h <- h+1
  }
}

sniffles_call <- rbindlist(sniffles_call)
sniffles_call[,SOFTWARE:="Sniffles"]

sniffles_base <- rbindlist(sniffles_base)
sniffles_base[,SOFTWARE:="Sniffles"]

sniffles_overlap <- rbindlist(sniffles_overlap)
sniffles_overlap[,SOFTWARE:="Sniffles"]

sniffles_vcf <- rbindlist(sniffles_vcf)

nano_cov <- rbindlist(nano_cov)
nano_cov_sniffles <- nano_cov[,.(median(V2)),by=list(FILE,COVERAGE)]

# NanoSV:

nanosv_call <- list()
nanosv_base <- list()
nanosv_overlap <- list()
nanosv_vcf <- list()

nano_cov <- list()

h <- 1
for (i in c('NA12878_SEQUIN_PROMETHION')) {
  for (j in as.character(1:20)) {
    info <- file.info(sprintf('curves/%s.sv.sorted.%s.nanosv.overlap.tab',i,j))
    if (info$size != 0) {
      dt1 <- fread(sprintf('curves/%s.sv.sorted.%s.nanosv.overlap.tab',i,j))
      dt1[,FILE:=sprintf('%s',i)]
      dt1[,COVERAGE:=j]
      nanosv_overlap[[h]] <- dt1
      
      dt2 <- fread(sprintf('curves/%s.sv.sorted.%s.nanosv.bed',i,j))
      dt2[,FILE:=sprintf('%s',i)]
      dt2[,COVERAGE:=j]
      nanosv_call[[h]] <- dt2
      
      dt3 <- fread('curves/sequin_ref_breakends.bed',header=FALSE)
      dt3[,FILE:=sprintf('%s',i)]
      dt3[,COVERAGE:=j]
      nanosv_base[[h]] <- dt3
      
      dt4 <- fread(sprintf('zcat curves/%s.sv.sorted.%s.nanosv.vcf.gz',i,j),skip="#CHROM")
      dt4[,FILE:=sprintf('%s',i)]
      dt4[,COVERAGE:=j]
      nanosv_vcf[[h]] <- dt4
      
      dt5 <- fread(sprintf('curves/%s.sv.sorted.%s.summary.tab',i,j))
      dt5[,FILE:=sprintf('%s',i)]
      dt5[,COVERAGE:=j]
      nano_cov[[h]] <- dt5
    }
    # dt1 <- data.table(V1=character(),V2=numeric(),V3=numeric(),V4=character(),
    #            V5=character(),V6=numeric(),V7=character(),V8=numeric(),
    #            V9=numeric(),V10=numeric(),V11=character(),V12=character(),
    #            V13=numeric(),FILE=character(),COVERAGE=character())
    h <- h+1
  }
}

nanosv_call <- rbindlist(nanosv_call)
nanosv_call[,SOFTWARE:="NanoSV"]

nanosv_base <- rbindlist(nanosv_base)
nanosv_base[,SOFTWARE:="NanoSV"]

nanosv_overlap <- rbindlist(nanosv_overlap)
nanosv_overlap[,SOFTWARE:="NanoSV"]

nanosv_vcf <- rbindlist(nanosv_vcf)

nano_cov <- rbindlist(nano_cov)
nano_cov_nanosv <- nano_cov[,.(median(V2)),by=list(FILE,COVERAGE)]

####################################################################
####################################################################

# I'll merge the two datasets so that I can process them together.

call <- rbind(lumpy_call,manta_call,sniffles_call,nanosv_call)
base <- rbind(lumpy_base,manta_base,sniffles_base,nanosv_base)
overlap <- rbind(lumpy_overlap,manta_overlap,sniffles_overlap,nanosv_overlap)

cov <- rbind(cov,nano_cov[,.(V1=median(V2)),by=list(FILE,COVERAGE)])

# Now I need to determine TP, FP and FN at the breakend level
# Based on my previous definition everythin in overlap should be TP
# Everything in call minus overlap will give FP
# Finally everything in base minus overlap will give FN

call[,ID1:=paste(V1,V2,V3,sep = ":")]
base[,ID2:=paste(V1,V2,V3,sep = ":")]
overlap[,ID1:=paste(V7,V8,V9,sep = ":")]
overlap[,ID2:=paste(V1,V2,V3,sep = ":")]

# call[SOFTWARE=="Sniffles" & V5=="BND" & COVERAGE==20]

# I'll create an index to every distintic row so I can easily the detect the ones that are in overlap or not

# I have to fix translocations for Sniffles
# sniffles_bnd <- sniffles_vcf[,.(FILE,COVERAGE,SVTYPE=str_remove(str_extract(INFO,"SVTYPE=([A-Z][A-Z][A-Z])"),"SVTYPE="),
#                                 ALT=str_remove_all(ALT,"[\\[N\\]]"))][SVTYPE=="BND"]

call[,N:=1:dim(call)[1]]
base[,N:=1:dim(base)[1]]

# Testing the overlap based on the IDs
tp1 <- merge(overlap[,.(FILE,COVERAGE,SOFTWARE,ID1)],call[,.(FILE,COVERAGE,SOFTWARE,N,ID1)],by=c("FILE","COVERAGE","SOFTWARE","ID1"),all_x=TRUE)
tp2 <- merge(overlap[,.(FILE,COVERAGE,SOFTWARE,ID2)],base[,.(FILE,COVERAGE,SOFTWARE,N,ID2)],by=c("FILE","COVERAGE","SOFTWARE","ID2"),all_x=TRUE)

fp <- call[!(N %in% tp1$N)]
fn <- base[!(N %in% tp2$N)]

tp <- overlap[,.(CHROM=V7,START=V8,END=V9,WINDOW=ID2,SEQUIN=V4,BREAKEND=V6,SVTYPE=V5,STATUS="TP",FILE,COVERAGE,SOFTWARE,ID=V10)]
fp <- fp[,.(CHROM=V1,START=V2,END=V3,WINDOW=NA,SEQUIN=NA,BREAKEND=NA,SVTYPE=V5,STATUS="FP",FILE,COVERAGE,SOFTWARE,ID=V4)]
fn <- fn[,.(CHROM=NA,START=NA,END=NA,WINDOW=ID2,SEQUIN=V4,BREAKEND=V6,SVTYPE=V5,STATUS="FN",FILE,COVERAGE,SOFTWARE,ID=NA)]

all <- rbind(tp,fp,fn)

# Sanity checks

teste1 <- all[STATUS %in% c("TP")]
teste1 <- teste1[,.N,by=list(FILE,COVERAGE,SOFTWARE,WINDOW)]
teste1 <- teste1[,.(FILE,COVERAGE,SOFTWARE)]
teste2 <- all[STATUS %in% c("FN"),.(FILE,COVERAGE,SOFTWARE)]
teste <- rbind(teste1,teste2)
teste <- teste[,.N,by=list(FILE,COVERAGE,SOFTWARE)]
teste[,COVERAGE:=as.numeric(as.character(COVERAGE))]

# The number of FN+TP in each sample at different coverages should be the same.
# For some reason, the different software will some times report the same variant twice
# in the case of translocations or inversions. 
# To deal with this I'll count a TP for a given breakend only once.

ggplot(teste,aes(COVERAGE,N))+geom_line()+facet_grid(FILE~SOFTWARE)

# I'll add true position of each original breakpoint

ref_pos <- fread("/simon/ANDRE/structural_variants/curves/sequin_structural_hg38_2.6.sample.breakend.bed")
ref_pos <- ref_pos[,.(POS=V2,SEQUIN=V4,BREAKEND=V6)]

all <- merge(all,ref_pos,by=c("SEQUIN","BREAKEND"),all.x=TRUE)

# Now I need to benchmark the True Positive breakends at the SV level.
# In this scenario, I'll consider the SV a true positive if both breakends were identified successfully,
# otherwise I'll classify the SV as a FALSE negative

sv_level <- unique(all[STATUS %in% c("TP","FN"),.(SEQUIN,BREAKEND,STATUS,SVTYPE,FILE,COVERAGE,SOFTWARE)])
sv_level <- data.table(dcast(sv_level,SEQUIN+SVTYPE+FILE+COVERAGE+SOFTWARE~BREAKEND,value.var = "STATUS"))

# Insertions have only one reference breakpoint, so I need to take that into account.

sv_level[,STATUS:="FN"]
sv_level[(`1`=="TP" & `2`=="TP") | (`1`=="TP" & is.na(`2`)),STATUS:="TP"]

sv_level <- merge(sv_level,cov,by=c("FILE","COVERAGE"),all.x=TRUE)

svtype <- data.table(dcast(sv_level[,.N,by=list(FILE,COVERAGE,SOFTWARE,SVTYPE,STATUS,V1)],
                           FILE+COVERAGE+V1+SOFTWARE+SVTYPE~STATUS,value.var = "N"))
svtype[is.na(svtype)] <- 0
svtype[,SN:=TP/(TP+FN)]

svtype$COVERAGE <- factor(svtype$COVERAGE,levels=1:20)

ggplot(svtype,aes(COVERAGE,FILE,fill=SN))+
  geom_tile()+
  facet_grid(SOFTWARE~SVTYPE,scales = "free_y",space = "free_y")+
  scale_fill_gradient2(midpoint = mid, low = "blue", mid = "white",
                       high = "red", space = "Lab",limits=c(0,1),breaks=c(0,0.5,1))+
  theme_bw(8)+
  scale_x_discrete(labels=c("5",rep("",8),"50",rep("",9),"100"))+
  labs(x="Coverage",y=NULL,fill=NULL)

svtype[,READ_LENGTH:="short-read"]
svtype[FILE=="NA12878_SEQUIN_PROMETHION",READ_LENGTH:="long-read"]

label <- svtype[,.(round(mean(SN),2),round(sd(SN),2)),by=list(SVTYPE)]
label[,LABEL:=paste("Mean = ",V1,"; SD = ",V2,sep="")]
label[,X:="Manta"]
label[,Y:=1.1]


ggplot(svtype,aes(SOFTWARE,SN))+
  geom_boxplot(aes(color=READ_LENGTH))+
  facet_wrap(~SVTYPE,nrow=1)+
  scale_y_continuous(limits=c(0,1.1),breaks=c(0,0.5,1))+
  geom_text(data=label,aes(x=X,y=Y,label=LABEL))+
  theme(aspect.ratio = 1)+
  labs(x="SV Type",y="Sensitivity")+
  guides(color=FALSE)

svtype$COV_LABEL <- factor(svtype$COVERAGE,levels=1:20)
levels(svtype$COV_LABEL) <- as.character(seq(5,100,5))
svtype$COV_LABEL <- as.numeric(as.character(svtype$COV_LABEL))
svtype[,QUANTILE:=cut(as.numeric(as.character(COV_LABEL)),breaks=c(0,20,40,60,80,100),labels=c("<20","20-40","40-60","60-80",">80"))]
svtype$PLATFORM <- as.factor(svtype$FILE)
levels(svtype$PLATFORM) <-  c("HiSeqXTen/PCR-free","BGI","ONT","HiSeq2500/PCR-based","HiSeq2500/PCR-based","NextSeq/PCR-based","NextSeq/PCR-based")

wilcox.test(svtype[SVTYPE=="DEL" & QUANTILE == "<20" & PLATFORM =="ONT" & SOFTWARE=="NanoSV",SN],svtype[SVTYPE=="DEL" & QUANTILE == "<20" & PLATFORM =="HiSeqXTen/PCR-free",SN])


ggplot(svtype[SVTYPE=="DEL" & QUANTILE == "<20"],aes(SOFTWARE,SN,color=SOFTWARE))+
  geom_boxplot()+
  facet_wrap(~PLATFORM,nrow=1,drop = TRUE,scales = "free_x")+
  theme(aspect.ratio = 1)+
  scale_y_continuous(limits=c(0,1),breaks=c(0,0.5,1))+
  labs(x=NULL,y="Sensitivity",color="Software")

ggplot(svtype,aes(PLATFORM,SN,color=SOFTWARE))+
  geom_boxplot()+
  facet_wrap(~SVTYPE,nrow=1)+
  theme(aspect.ratio = 1)+
  scale_y_continuous(limits=c(0,1),breaks=c(0,0.5,1))+
  labs(x=NULL,y="Sensitivity",color="Software")


  

svlen <- fread('/simon/ANDRE/structural_variants/sequin_structural_hg38_2.6.svlen.tab')
names(svlen) <- c('SEQUIN','LENGTH') 
svlen <- merge(svlen,sv_level,by="SEQUIN")
svlen[,LENGTH:=abs(LENGTH)]

ggplot(unique(svlen[,.(SEQUIN,SVTYPE,LENGTH)]),aes(SVTYPE,LENGTH))+
  geom_violin(draw_quantiles = c(0.5),width=0.8)+
  geom_jitter(width=0.2,alpha=0.5,size=4)+
  theme(aspect.ratio = 1)+
  labs(x="SV Type",y="Length")

svlen$COVERAGE <- factor(svlen$COVERAGE,levels=1:20)
svlen$PLATFORM <- as.factor(svlen$FILE)
levels(svlen$PLATFORM) <-  c("HiSeqXTen/PCR-free","BGI","ONT","HiSeq2500/PCR-based","HiSeq2500/PCR-based","NextSeq/PCR-based","NextSeq/PCR-based")
svlen$STATUS <- as.factor(svlen$STATUS)
levels(svlen$STATUS) <- c("False-negative","True-positive")

ggplot(svlen[as.numeric(COVERAGE)<=8],aes(COVERAGE,LENGTH,by=FILE,color=PLATFORM))+
  geom_boxplot()+
  facet_grid(SVTYPE~STATUS)+
  scale_x_discrete(labels=seq(5,40,5))+
  scale_color_manual(values=c("HiSeqXTen/PCR-free"="#E185B7","HiSeq2500/PCR-based"="#EB8D82","NextSeq/PCR-based"="#00AEEB","BGI"="#ADAB33","ONT"="#00AB8A"))+
  labs(x="Coverage",y="Length",color=NULL)+
  guides(color=FALSE)

svlen[,READ_LENGTH:="short-read"]
svlen[FILE=="NA12878_SEQUIN_PROMETHION",READ_LENGTH:="long-read"]

ggplot(svlen[as.numeric(COVERAGE) %in% c(2,4,6,8,10,12,14) & SVTYPE == "INS"],aes(COVERAGE,LENGTH,color=PLATFORM))+
  geom_boxplot()+
  #geom_point()+
  #geom_smooth()+
  facet_grid(~STATUS)+
  scale_x_discrete(labels=c(10,20,30,40,50,60,70))+
  theme(aspect.ratio = 0.5)+
  #scale_x_discrete(labels=seq(5,40,5))+
  scale_color_manual(values=c("HiSeqXTen/PCR-free"="#E185B7","HiSeq2500/PCR-based"="#EB8D82","NextSeq/PCR-based"="#00AEEB","BGI"="#ADAB33","ONT"="#00AB8A"))+
  labs(x="Coverage",y="Length",color=NULL)+
  guides(color=FALSE)
  

ggplot(svtype,aes(V1,SN,color=FILE))+geom_line()+facet_grid(SOFTWARE~SVTYPE)

# Sanity check
teste <- sv_level[,.N,by=list(FILE,COVERAGE,SOFTWARE)]
table(teste$N)

# Now I have two datasets to benchmark and make plots. One at the breakpoint level and one at the SV level.

breakend_level <- data.table(all)
breakend_level <- merge(breakend_level,cov,by=c("FILE","COVERAGE"),all.x=TRUE)

stats <- breakend_level[,.N,by=list(FILE,SOFTWARE,STATUS)]
ggplot(stats,aes(FILE,N,fill=STATUS))+
  geom_bar(stat="identity",alpha=0.6)+
  facet_grid(~SOFTWARE,scales = "free_x",space="free_x")+
  labs(x=NULL,y="Frequency",fill=NULL)+
  theme_bw(12)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  scale_y_continuous(breaks=c(0,3000,6000))


# Distance to true breakpoint
dist <- breakend_level[STATUS %in% c("TP")]
dist[,DIST:=START-POS]
per_tech_dist <- data.table(dist)
dist <- dist[,.N,by=list(FILE,SOFTWARE,DIST)]

dist$PLATFORM <- as.factor(dist$FILE)
levels(dist$PLATFORM) <-  c("HiSeqXTen/PCR-free","BGI","ONT","HiSeq2500/PCR-based","HiSeq2500/PCR-based","NextSeq/PCR-based","NextSeq/PCR-based")

# ggplot(dist[DIST>-5 & DIST<5],aes(as.factor(DIST),N,by=FILE,fill=PLATFORM))+
#   geom_bar(color="black",stat="identity",position = position_dodge())+
#   facet_grid(~SOFTWARE)+
#   scale_fill_manual(values=c("HiSeqXTen/PCR-free"="#E185B7","HiSeq2500/PCR-based"="#EB8D82","NextSeq/PCR-based"="#00AEEB","BGI"="#ADAB33","ONT"="#00AB8A"))+
#   theme_bw(12)+
#   theme(#aspect.ratio = 1,
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank())+
#   guides(fill=FALSE)+
#   scale_y_continuous(breaks=c(0,300,600))+
#   labs(x="Distance to true breakpoint",y="Frequency")
total <- dist[,.(TOTAL=sum(N)),by=FILE]

dist <- merge(dist,total,by="FILE")

dist[,N:=N/TOTAL]

dist_rel_freq <- dist[DIST %in% c(-1,0,1),.(sum(N)),by=list(FILE,SOFTWARE)]
dist_rel_freq[,.(Mean=mean(V1),SD=sd(V1)),by=SOFTWARE]

ggplot(dist[DIST>-5 & DIST<5],aes(DIST,N,by=FILE,color=PLATFORM))+
  geom_point(alpha=0.6,size=2)+
  facet_grid(~SOFTWARE)+
  scale_color_manual(values=c("HiSeqXTen/PCR-free"="#E185B7","HiSeq2500/PCR-based"="#EB8D82","NextSeq/PCR-based"="#00AEEB","BGI"="#ADAB33","ONT"="#00AB8A"))+
  theme_bw(12)+
  theme(aspect.ratio = 1,
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())+
  guides(color=FALSE)+
  #scale_y_continuous(breaks=c(0,0.5,1))+
  labs(x="Distance to true breakpoint",y="Frequency")
  
dist <- breakend_level[STATUS %in% c("TP")]
dist[,DIST:=START-POS]
dist <- dist[DIST==0]
dist <- dist[,.N,by=list(FILE,COVERAGE,SOFTWARE,V1)]

dist$PLATFORM <- as.factor(dist$FILE)
levels(dist$PLATFORM) <- c("HiSeqXTen/PCR-free","BGI","ONT","HiSeq2500/PCR-based","HiSeq2500/PCR-based","NextSeq/PCR-based","NextSeq/PCR-based")

ggplot(dist,aes(V1,N,by=FILE,color=PLATFORM))+
  geom_line()+
  facet_grid(~SOFTWARE)+
  scale_color_manual(values=c("HiSeqXTen/PCR-free"="#E185B7","HiSeq2500/PCR-based"="#EB8D82","NextSeq/PCR-based"="#00AEEB","BGI"="#ADAB33","ONT"="#00AB8A"))+
  theme_bw(12)+
  theme(aspect.ratio = 1,
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  labs(x="Coverage",y="Frequency of exact breakpoints")+
  scale_x_continuous(breaks=c(0,50,100))+
  scale_y_continuous(breaks=c(0,20,40))+
  guides(color=FALSE)

##### Per tech dist to original breakpoint
tech_dist <- per_tech_dist[,.(FILE,COVERAGE,BREAKEND,WINDOW,STATUS,SOFTWARE,DIST)]
unique_tech_dist <- unique(tech_dist[STATUS=="TP",.(FILE,COVERAGE,WINDOW,SOFTWARE,DIST)])
unique_tech_dist <- unique_tech_dist[,.(DIST=min(DIST)),by=list(FILE,COVERAGE,WINDOW,SOFTWARE)]
tech_dist <- data.table(dcast(unique_tech_dist,FILE+COVERAGE+WINDOW~SOFTWARE,value.var = "DIST"))
tech_dist[,DIST:=min(Lumpy,Manta,NanoSV,Sniffles,na.rm=TRUE),by=list(FILE,COVERAGE,WINDOW)]
tech_dist <- tech_dist[,.N,by=list(FILE,DIST)]
total <- tech_dist[,.(TOTAL=sum(N)),by=FILE]
tech_dist <- merge(tech_dist,total,by="FILE")
tech_dist[,N:=N/TOTAL]

tech_dist$PLATFORM <- as.factor(tech_dist$FILE)
levels(tech_dist$PLATFORM) <- c("HiSeqXTen/PCR-free","BGI","ONT","HiSeq2500/PCR-based","HiSeq2500/PCR-based","NextSeq/PCR-based","NextSeq/PCR-based")

ggplot(tech_dist,aes(DIST,N,by=FILE,color=PLATFORM))+
  geom_line()+
  geom_point(alpha=0.6)+
  scale_color_manual(values=c("HiSeqXTen/PCR-free"="#E185B7","HiSeq2500/PCR-based"="#EB8D82","NextSeq/PCR-based"="#00AEEB","BGI"="#ADAB33","ONT"="#00AB8A"))+
  theme(aspect.ratio = 1)+
  labs(x="Distance to true breakpoint",y="Relative frequency")+
  guides(color=FALSE)

tech_dist[DIST>=-5 & DIST <= 5,sum(N),by=FILE]
tech_dist[DIST>=-5 & DIST <= 5 & FILE != "NA12878_SEQUIN_PROMETHION",sum(N),by=FILE][,.(mean(V1),sd(V1))]

# First I'll calculate the basic statistics over coverage
# Each true breakpoint has to be considered only once.

all_tp <- breakend_level[STATUS %in% c("TP")]
all_tp <- all_tp[,.N,by=list(FILE,COVERAGE,V1,SOFTWARE,WINDOW)]
all_tp <- all_tp[,.(FILE,COVERAGE,V1,SOFTWARE,WINDOW,STATUS="TP")]
all_fp <- breakend_level[STATUS %in% c("FP")]
all_fp <- all_fp[,.N,by=list(FILE,COVERAGE,V1,SOFTWARE,ID)]
all_fp <- all_fp[,.(FILE,COVERAGE,V1,SOFTWARE,WINDOW=NA,STATUS="FP")]
all_fn<- breakend_level[STATUS %in% c("FN"),.(FILE,COVERAGE,V1,SOFTWARE,WINDOW,STATUS="FN")]
stats <- rbind(all_tp,all_fp,all_fn)

# Performance per technology, regardless of software.
per_tech <- data.table(stats)
per_tech$COV_LABEL <- factor(per_tech$COVERAGE,levels=1:20)
levels(per_tech$COV_LABEL) <- as.character(seq(5,100,5))
per_tech[,QUANTILE:=cut(as.numeric(as.character(COV_LABEL)),breaks=c(0,25,50,75,100),labels=c("<25","25-50","50-75",">75"))]
per_tech[,READ_LENGTH:="short-read"]
per_tech[FILE=="NA12878_SEQUIN_PROMETHION",READ_LENGTH:="long-read"]
per_tech$PLATFORM <- as.factor(per_tech$FILE)
levels(per_tech$PLATFORM) <-  c("HiSeqXTen/PCR-free","BGI","ONT","HiSeq2500/PCR-based","HiSeq2500/PCR-based","NextSeq/PCR-based","NextSeq/PCR-based")

per_software <- per_tech[,.N,by=list(FILE,PLATFORM,COVERAGE,COV_LABEL,QUANTILE,READ_LENGTH,STATUS,SOFTWARE)]
per_software <- data.table(dcast(per_software,FILE+PLATFORM+COVERAGE+COV_LABEL+QUANTILE+READ_LENGTH+SOFTWARE~STATUS,value.var = "N"))
per_software[is.na(per_software)] <- 0

auc <- per_software[,round(MESS::auc(as.numeric(as.character(COV_LABEL))/100,TP/(TP+FP)),2),by=list(FILE,PLATFORM,SOFTWARE)]
auc[,LABEL:=paste("AUC=(",as.character(V1),")",sep="")]
auc[SOFTWARE %in% c("Lumpy","NanoSV"),Y:=0.66]
auc[SOFTWARE %in% c("Manta","Sniffles"),Y:=0.69]
auc[,X:=75]

ggplot(per_software[PLATFORM %in% c("HiSeqXTen/PCR-free","BGI","ONT") | FILE %in% c("SEQUIN_HISEQ_2500_PCR_BASED_2","SEQUIN_NEXTSEQ_500_PCR_BASED_1")],aes(as.numeric(as.character(COV_LABEL)),TP/(TP+FP),by=FILE,color=SOFTWARE))+
  geom_line()+
  facet_wrap(~PLATFORM,nrow=1)+
  scale_x_continuous(limits=c(0,100),breaks=c(0,50,100))+
  scale_y_continuous(limits=c(0.65,1))+
  #geom_vline(xintercept = 55)+
  geom_text(data=auc[PLATFORM %in% c("HiSeqXTen/PCR-free","BGI","ONT") | FILE %in% c("SEQUIN_HISEQ_2500_PCR_BASED_2","SEQUIN_NEXTSEQ_500_PCR_BASED_1")],
            aes(x=X,y=Y,label=LABEL))+
  theme(aspect.ratio = 1)+
  labs(x="Coverage",y="Precision",color="Software")

per_tech <- data.table(dcast(unique(per_tech[,.(FILE,PLATFORM,COVERAGE,COV_LABEL,QUANTILE,READ_LENGTH,WINDOW,SOFTWARE,STATUS)]),FILE+PLATFORM+COVERAGE+COV_LABEL+QUANTILE+READ_LENGTH+WINDOW~SOFTWARE,value.var = "STATUS"))
per_tech[,STATUS:="FN"]
per_tech[(Lumpy=="TP") | (Manta=="TP") | (NanoSV=="TP") | (Sniffles=="TP"),STATUS:="TP"]
per_tech[ STATUS!="TP" & ((Lumpy=="FP") | (Manta=="FP") | (NanoSV=="FP") | (Sniffles=="FP")),STATUS:="FP"]
per_tech <- per_tech[,.N,by=list(FILE,PLATFORM,COVERAGE,COV_LABEL,QUANTILE,READ_LENGTH,STATUS)]
per_tech <- data.table(dcast(per_tech,FILE+PLATFORM+COVERAGE+COV_LABEL+QUANTILE+READ_LENGTH~STATUS,value.var = "N"))
per_tech[is.na(per_tech)] <- 0
per_tech[,SN:=TP/(TP+FN)]
per_tech[,PR:=TP/(TP+FP)]

ggplot(per_tech,aes(as.numeric(as.character(COV_LABEL)),SN,by=FILE,color=PLATFORM))+
  geom_line()+
  #geom_roc()+
  scale_color_manual(values=c("HiSeqXTen/PCR-free"="#E185B7","HiSeq2500/PCR-based"="#EB8D82","NextSeq/PCR-based"="#00AEEB","BGI"="#ADAB33","ONT"="#00AB8A"))+
  theme_bw(12)+
  theme(aspect.ratio = 1,
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  labs(x="Coverage",y="Sensitivity")+
  guides(color=FALSE)



per_tech[,.(mean(SN),sd(SN)),by=list(QUANTILE,READ_LENGTH)]
per_tech[,.(mean(SN),sd(SN)),by=list(QUANTILE,PLATFORM)]

pcr_based <- per_tech[PLATFORM %in% c("HiSeq2500/PCR-based","NextSeq/PCR-based")]
pcr_based$PLATFORM <- as.factor(as.character(pcr_based$PLATFORM))
wilcox.test(pcr_based[QUANTILE=="<25" & PLATFORM == "HiSeq2500/PCR-based",SN],pcr_based[QUANTILE=="<25" & PLATFORM=="NextSeq/PCR-based",SN])
wilcox.test(pcr_based[QUANTILE=="25-50" & PLATFORM == "HiSeq2500/PCR-based",SN],pcr_based[QUANTILE=="25-50" & PLATFORM=="NextSeq/PCR-based",SN])
wilcox.test(pcr_based[QUANTILE=="50-75" & PLATFORM == "HiSeq2500/PCR-based",SN],pcr_based[QUANTILE=="50-75" & PLATFORM=="NextSeq/PCR-based",SN])
wilcox.test(pcr_based[QUANTILE==">75" & PLATFORM == "HiSeq2500/PCR-based",SN],pcr_based[QUANTILE==">75" & PLATFORM=="NextSeq/PCR-based",SN])


format(wilcox.test(per_tech[QUANTILE=="25-50",SN],per_tech[QUANTILE=="25-50",PLATFORM])$p.value,scientific = FALSE)
format(pairwise.wilcox.test(per_tech[QUANTILE=="50-75",SN],per_tech[QUANTILE=="50-75",PLATFORM])$p.value,scientific = FALSE)
format(pairwise.wilcox.test(per_tech[QUANTILE==">75",SN],per_tech[QUANTILE==">75",PLATFORM])$p.value,scientific = FALSE)


ggplot(per_tech,aes(QUANTILE,SN,by=FILE,color=PLATFORM))+
  geom_boxplot()+
  scale_color_manual(values=c("HiSeqXTen/PCR-free"="#E185B7","HiSeq2500/PCR-based"="#EB8D82","NextSeq/PCR-based"="#00AEEB","BGI"="#ADAB33","ONT"="#00AB8A"))+
  theme(aspect.ratio = 1)+
  labs(x="Coverage",y="Sensitivity")+
  guides(color=FALSE)

ggplot(per_tech,aes(QUANTILE,PR,by=FILE,color=PLATFORM))+
  geom_boxplot()+
  facet_wrap(~READ_LENGTH)+
  scale_color_manual(values=c("HiSeqXTen/PCR-free"="#E185B7","HiSeq2500/PCR-based"="#EB8D82","NextSeq/PCR-based"="#00AEEB","BGI"="#ADAB33","ONT"="#00AB8A"))+
  theme(aspect.ratio = 1)+
  labs(x="Coverage",y="Precision")+
  guides(color=FALSE)


auc_per_tech <- per_tech[,auc(as.numeric(as.character(COV_LABEL))/100,SN),by=list(FILE,PLATFORM)]


ggplot(per_tech,aes(as.numeric(as.character(COV_LABEL)),PR,by=FILE,color=PLATFORM))+
  geom_line()+
  #geom_roc()+
  scale_color_manual(values=c("HiSeqXTen/PCR-free"="#E185B7","HiSeq2500/PCR-based"="#EB8D82","NextSeq/PCR-based"="#00AEEB","BGI"="#ADAB33","ONT"="#00AB8A"))+
  theme_bw(12)+
  theme(aspect.ratio = 1,
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  scale_y_continuous(limits = c(0,1),breaks=c(0,0.5,1))+
  labs(x="Coverage",y="Precision")+
  guides(color=FALSE)


stats <- stats[,.N,by=list(FILE,COVERAGE,V1,SOFTWARE,STATUS)]
stats <- data.table(dcast(stats,FILE+COVERAGE+V1+SOFTWARE~STATUS,value.var = "N"))
stats[is.na(stats)] <- 0
stats[,SN:=TP/(TP+FN)]
stats[,PR:=TP/(TP+FP)]

stats$PLATFORM <- as.factor(stats$FILE)
levels(stats$PLATFORM) <-  c("HiSeqXTen/PCR-free","BGI","ONT","HiSeq2500/PCR-based","HiSeq2500/PCR-based","NextSeq/PCR-based","NextSeq/PCR-based")

ggplot(stats,aes(V1,SN,by=FILE,color=PLATFORM))+
  geom_line()+
  #geom_roc()+
  facet_grid(~SOFTWARE)+
  scale_color_manual(values=c("HiSeqXTen/PCR-free"="#E185B7","HiSeq2500/PCR-based"="#EB8D82","NextSeq/PCR-based"="#00AEEB","BGI"="#ADAB33","ONT"="#00AB8A"))+
  theme_bw(12)+
  theme(aspect.ratio = 1,
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  scale_x_continuous(breaks=c(0,50,100))+
  scale_y_continuous(breaks=c(0,0.3,0.6))+
  labs(x="Coverage",y="Sensitivity")+
  guides(color=FALSE)

auc <- stats[,round(MESS::auc(V1/100,SN),2),by=list(FILE,PLATFORM,SOFTWARE)]
auc[,LABEL:=paste("AUC=(",as.character(V1),")",sep="")]
auc[SOFTWARE %in% c("Lumpy","NanoSV"),Y:=0.2]
auc[SOFTWARE %in% c("Manta","Sniffles"),Y:=0.27]
auc[,X:=75]

ggplot(stats[PLATFORM %in% c("HiSeqXTen/PCR-free","BGI","ONT") | FILE %in% c("SEQUIN_HISEQ_2500_PCR_BASED_2","SEQUIN_NEXTSEQ_500_PCR_BASED_1")],aes(V1,SN,by=FILE,color=SOFTWARE))+
  geom_line()+
  facet_wrap(~PLATFORM,nrow=1)+
  theme(aspect.ratio = 1)+
  geom_text(data=auc[PLATFORM %in% c("HiSeqXTen/PCR-free","BGI","ONT") | FILE %in% c("SEQUIN_HISEQ_2500_PCR_BASED_2","SEQUIN_NEXTSEQ_500_PCR_BASED_1")],
            aes(x=X,y=Y,label=LABEL))+
  scale_x_continuous(limits=c(0,100),breaks=c(0,50,100))+
  labs(x="Coverage",y="Sensitivity",color="Software")

auc[SOFTWARE %in% c("Lumpy","NanoSV"),V2:="A1"]
auc[SOFTWARE %in% c("Manta","Sniffles"),V2:="A2"]
auc[SOFTWARE %in% c("Lumpy","NanoSV"),V3:="B1"]
auc[SOFTWARE %in% c("Manta","Sniffles"),V3:="B2"]

auc_sgment <- merge(dcast(auc[,.(FILE,PLATFORM,V2,V1)],FILE+PLATFORM~V2,value.var = "V1"),
      dcast(auc[,.(FILE,PLATFORM,V3,SOFTWARE)],FILE+PLATFORM~V3,value.var = "SOFTWARE"),by=c("FILE","PLATFORM"))

ggplot(auc,aes(SOFTWARE,V1,by=FILE,color=PLATFORM))+
  geom_point(alpha=0.6,size=2)+
  geom_segment(data=auc_sgment,aes(x=B1,y=A1,xend=B2,yend=A2))+
  scale_color_manual(values=c("HiSeqXTen/PCR-free"="#E185B7","HiSeq2500/PCR-based"="#EB8D82","NextSeq/PCR-based"="#00AEEB","BGI"="#ADAB33","ONT"="#00AB8A"))+
  theme(aspect.ratio = 1)+
  scale_y_continuous(limits=c(0.3,0.6),breaks=c(0.3,0.45,0.6))+
  labs(x=NULL,y="AUC")+
  guides(color=FALSE)

  

ggplot(stats,aes(V1,PR,by=FILE,color=PLATFORM))+
  geom_line()+
  facet_grid(~SOFTWARE)+
  scale_color_manual(values=c("HiSeqXTen/PCR-free"="#E185B7","HiSeq2500/PCR-based"="#EB8D82","NextSeq/PCR-based"="#00AEEB","BGI"="#ADAB33","ONT"="#00AB8A"))+
  theme_bw(12)+
  theme(aspect.ratio = 1,
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  scale_x_continuous(breaks=c(0,50,100))+
  scale_y_continuous(limits=c(0,1),breaks=c(0,0.5,1))+
  labs(x="Coverage",y="Precision")+
  guides(color=FALSE)


# I also need to extract other quality information about the call
lumpy <- breakend_level[SOFTWARE=="Lumpy"]
lumpy_qual <- data.table(lumpy_vcf)
info <- strsplit(lumpy_vcf$FORMAT[1],":",fixed = TRUE)[[1]]
lumpy_qual[,c("GT","SU","PE","SR") := tstrsplit(s1,":",fixed = TRUE)]
lumpy_qual[,SU:=as.numeric(SU)]
lumpy_qual[,PE:=as.numeric(PE)]
lumpy_qual[,SR:=as.numeric(SR)]
lumpy_qual <- lumpy_qual[,.(CHROM=`#CHROM`,FILE,COVERAGE,ID,GT,SU,PE,SR)]
lumpy <- merge(lumpy,lumpy_qual,by=c("CHROM","FILE","COVERAGE","ID"),all.x=TRUE)

manta <- breakend_level[SOFTWARE=="Manta"]
manta_qual <- data.table(manta_vcf)
manta_qual[,c("GT","FT","GQ","PL","PR","SR") := tstrsplit(s1,":",fixed = TRUE)]
manta_qual[,GQ:=as.numeric(GQ)]
manta_qual[,c("PL.1","PL.2","PL.3"):=tstrsplit(PL,",",fixed = TRUE)]
manta_qual[,PL.1:=as.numeric(PL.1)]
manta_qual[,PL.2:=as.numeric(PL.2)]
manta_qual[,PL.3:=as.numeric(PL.3)]
manta_qual[,c("PR.REF","PR.ALT"):= tstrsplit(PR,",",fixed = TRUE)]
manta_qual[,PR.REF:=as.numeric(PR.REF)]
manta_qual[,PR.ALT:=as.numeric(PR.ALT)]
manta_qual[,c("SR.REF","SR.ALT"):= tstrsplit(SR,",",fixed = TRUE)]
manta_qual[,SR.REF:=as.numeric(SR.REF)]
manta_qual[,SR.ALT:=as.numeric(SR.ALT)]
manta_qual <- manta_qual[,.(CHROM=`#CHROM`,FILE,COVERAGE,ID,GQ,PL.1,PL.2,PL.3,PR.REF,PR.ALT,SR.REF,SR.ALT)]
manta <- merge(manta,manta_qual,by=c("CHROM","FILE","COVERAGE","ID"),all.x=TRUE)

####################

sniffles <- breakend_level[SOFTWARE=="Sniffles"]
sniffles_qual <- data.table(sniffles_vcf)
sniffles_qual[,c("GT","DR","DV") := tstrsplit(`L.N.176.sv.sorted.2.bam`,":",fixed = TRUE)]
sniffles_qual[,DV:=as.numeric(DV)]
sniffles_qual <- sniffles_qual[,.(CHROM=`#CHROM`,FILE,COVERAGE,ID=as.character(ID),GT,DR,DV)]
sniffles <- merge(sniffles,sniffles_qual,by=c("CHROM","FILE","COVERAGE","ID"),all.x=TRUE)

nanosv <- breakend_level[SOFTWARE=="NanoSV"]
nanosv_qual <- data.table(nanosv_vcf)
nanosv_qual[,c('GT','DR','DV','GQ','HR','PL') := tstrsplit(s1,":",fixed = TRUE)]
nanosv_qual[,c('DV.1','DV.2') := tstrsplit(DV,",",fixed = TRUE)]
nanosv_qual[,DV.1:=as.numeric(DV.1)]
nanosv_qual[,DV.2:=as.numeric(DV.2)]
nanosv_qual <- nanosv_qual[,.(CHROM=`#CHROM`,FILE,COVERAGE,ID=as.character(ID),GT,DR,DV.1,DV.2,GQ,HR,PL)]
nanosv <- merge(nanosv,nanosv_qual,by=c("CHROM","FILE","COVERAGE","ID"),all.x=TRUE)

# Just added the real coverage values for each sample as V1
# lumpy <- merge(lumpy,cov,by=c("FILE","COVERAGE"),all.x=TRUE)
# manta <- merge(manta,cov,by=c("FILE","COVERAGE"),all.x=TRUE)

qual <- rbind(lumpy[STATUS %in% c("TP","FP"),.(FILE,COVERAGE,V1,SOFTWARE,SU,PE,SR,WINDOW,ID,STATUS)],
              manta[STATUS %in% c("TP","FP"),.(FILE,COVERAGE,V1,SOFTWARE,SU=PR.ALT+SR.ALT,PE=PR.ALT,SR=SR.ALT,WINDOW,ID,STATUS)])
              # sniffles[STATUS %in% c("TP","FP"),.(FILE,COVERAGE,V1,SOFTWARE,SU=DV,WINDOW,ID,STATUS)],
              # nanosv[STATUS %in% c("TP","FP"),.(FILE,COVERAGE,V1,SOFTWARE,SU=DV.1,WINDOW,ID,STATUS)])

su <- qual[STATUS=="TP" & SOFTWARE %in% c("Lumpy","Manta"),.(FILE,COVERAGE,WINDOW,SOFTWARE,SR)]
su[,N1:=.N,by=list(FILE,COVERAGE,WINDOW)]
su[,N2:=.N,by=list(FILE,COVERAGE,WINDOW,SOFTWARE)]

su <- data.table(dcast(su[N1==2 & N2==1,],FILE+COVERAGE+WINDOW~SOFTWARE,value.var = "SR"))

lm_eqn <- function(df){
  m <- lm(Manta ~ Lumpy, df);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));
}

ggplot(su,aes(Lumpy,Manta))+
  #geom_smooth()+
  #geom_point(alpha=0.6)+
  geom_point(alpha=0.6)+
  geom_smooth(method="lm")+
  geom_text(x = 20, y = 80, label = lm_eqn(su), parse = TRUE,check_overlap = TRUE)+
  geom_abline(slope=1,color="red")+
  theme(aspect.ratio = 1)

ggplot(qual,aes(STATUS,SU,color=FILE))+
  geom_boxplot()+
  facet_grid(~SOFTWARE)+
  theme_bw(12)+
  theme(aspect.ratio = 1,
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  guides(color=FALSE)+
  labs(x=NULL,y="Number of supporting reads")+
  scale_y_continuous(breaks=c(0,100,200))+
  scale_x_discrete(labels=c("False-positive","True-positive"))

qual <- qual[,.(MEAN=mean(SU,na.rm=TRUE),SD=sd(SU,na.rm=TRUE)),by=list(FILE,COVERAGE,V1,SOFTWARE)]

qual$PLATFORM <- as.factor(qual$FILE)
levels(qual$PLATFORM) <- c("HiSeqXTen/PCR-free","BGI","ONT","HiSeq2500/PCR-based","HiSeq2500/PCR-based","NextSeq/PCR-based","NextSeq/PCR-based")
  
ggplot(qual,aes(V1,MEAN,by=FILE,color=PLATFORM))+
  geom_line()+
  facet_grid(~SOFTWARE)+
  scale_color_manual(values=c("HiSeqXTen/PCR-free"="#E185B7","HiSeq2500/PCR-based"="#EB8D82","NextSeq/PCR-based"="#00AEEB","BGI"="#ADAB33","ONT"="#00AB8A"))+
  theme_bw(12)+
  theme(aspect.ratio = 1,
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  scale_x_continuous(breaks=c(0,50,100))+
  scale_y_continuous(breaks=c(0,30,60))+
  labs(x="Coverage",y="Average number of supporting reads")+
  guides(color=FALSE)
  
dev.off()