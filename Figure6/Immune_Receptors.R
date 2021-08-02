#setwd('/simon/ANDRE/Figures/Figure6/cdr3')
setwd('/simon/ANDRE/Figures/github/Figure6/Immune_Receptors')
library(Metrics)
pdf('Immune_Receptors.pdf',height=8,width=12,useDingbats = FALSE)

all.files <- list.files(pattern = ".split_immune.no_edges.sorted.sub.full.all_metrics.tab")
l <- lapply(all.files, fread, header=TRUE)
all.files <- str_remove(all.files,".split_immune.no_edges.sorted.sub.full.all_metrics.tab")
l <- lapply(seq_along(all.files), function(x,y,i) x[[i]][,FILE:=y[i]],x=l,y=all.files)
dt <- rbindlist(l)
dt[,LABEL:="Rearranged"]
dt[endsWith(CHROM,"_R"),LABEL:="Background"]

consensus <- dt[,.(NAME,CHROM,POS,READS_ALL,REF=REF.x,A,C,G,T,MATCHES,MISMATCHES,INSERTIONS,DELETIONS,FILE)]
consensus[DELETIONS>MATCHES,CONSENSUS:="DEL"]
consensus[INSERTIONS>MATCHES,CONSENSUS:="INS"]
consensus[MISMATCHES>MATCHES & is.na(CONSENSUS),CONSENSUS:="MIS"]

consensus[is.na(CONSENSUS),]


# What is the variation in coverage at CDR3?
cov <- dt[,.(NAME,CHROM,POS,READS_ALL,FILE,LABEL)]
cov <- cov[,.(MEAN=median(READS_ALL),SD=sd(READS_ALL)),by=list(FILE,NAME,CHROM,LABEL)]
ggplot(cov[LABEL=="Rearranged"],aes(FILE,MEAN))+
  geom_boxplot()+
  scale_y_continuous(limits=c(0,100))

# What is the variation in mapping quality at CDR3?
mapq <- dt[,.(NAME,CHROM,POS,READS_MAPQ0,FILE,LABEL)]
mapq <- mapq[,.(MEAN=mean(READS_MAPQ0),SD=sd(READS_MAPQ0)),by=list(FILE,NAME,CHROM,LABEL)]
ggplot(mapq[LABEL=="Rearranged"],aes(FILE,MEAN))+
  geom_boxplot()+
  scale_y_continuous(limits=c(0,100))

# What is the variation in base quality at CDR3?
baseq <- dt[,.(NAME,CHROM,POS,RMS_BASEQ,RMS_BASEQ_MATCHES,RMS_BASEQ_MISMATCHES,FILE)]
baseq <- baseq[,.(RMS_BASEQ=mean(RMS_BASEQ),RMS_BASEQ_MATCHES=mean(RMS_BASEQ_MATCHES),RMS_BASEQ_MISMATCHES=mean(RMS_BASEQ_MISMATCHES)),by=list(FILE,NAME,CHROM)]

ggplot(baseq,aes(FILE,RMS_BASEQ))+
  geom_boxplot()+
  scale_y_continuous(limits=c(0,100))

ggplot(baseq,aes(FILE,RMS_BASEQ_MATCHES))+
  geom_boxplot()

ggplot(baseq,aes(FILE,RMS_BASEQ_MISMATCHES))+
  geom_boxplot()

# What is the error rate?
error <- dt[,.(NAME,CHROM,POS,READS_ALL,MISMATCHES,INSERTIONS,DELETIONS,FILE)]
error[,LENGTH:=nchar(NAME)]
error <- error[,ERROR:=(MISMATCHES+INSERTIONS+DELETIONS)]
error <- error[,.(ERROR=sum(ERROR),READS_ALL=median(READS_ALL)),by=list(FILE,NAME,CHROM,LENGTH)]
error[,ERROR:=ERROR/(READS_ALL*LENGTH)]
error[,LABEL:=tstrsplit(CHROM,"_",fixed=TRUE)[4]]

ggplot(error,aes(FILE,log(ERROR+1)))+
  geom_boxplot()
error[ERROR==0,ERROR:=NA]
ref <- fread('immune_sequins_ref.tab')
ref <- ref[,.(V1,V3,V5)]
names(ref) <- c("SEQUIN","V3","AF")
ref[grepl("TRD",V3),RECEPTOR:="TRD"]
ref[grepl("TRA",V3),RECEPTOR:="TRA"]
ref[grepl("TRB",V3),RECEPTOR:="TRB"]
ref[grepl("TRG",V3),RECEPTOR:="TRG"]
ref[grepl("IGK",V3),RECEPTOR:="IGK"]
ref[grepl("IGL",V3),RECEPTOR:="IGL"]


error <- merge(error,ref[,.(CHROM=SEQUIN,RECEPTOR)],by="CHROM")
error$PLATFORM <- as.factor(as.character(error$FILE))
levels(error$PLATFORM) <- c("HiSeqXTen/PCR-free","BGI","ONT","HiSeq2500/PCR-based","HiSeq2500/PCR-based","NextSeq/PCR-based","NextSeq/PCR-based")

ggplot(error,aes(RECEPTOR,FILE,by=CHROM,fill=log(ERROR)))+
  geom_tile()+
  facet_wrap(~LABEL,drop = TRUE,scales = "free_x",nrow=1)+
  #scale_x_discrete(labels=error$RECEPTOR)+
  #scale_fill_gradient2(low = "blue",mid="white",high = "red", na.value="gray", name = "" )+
  scale_fill_gradient2(midpoint = -4, low = "blue", mid = "white",
                       high = "red", space = "Lab",name=NULL,breaks=c(-2,-4,-6))+
  theme(axis.text.x = element_text(angle = 45,hjust = 1),aspect.ratio = 1,legend.position ="bottom")+
  labs(x=NULL,y=NULL)

error[,.(mean(ERROR*100,na.rm=TRUE),sd(ERROR*100,na.rm=TRUE)),by=PLATFORM]
  

# What is the count ratio between rearranged receptors?
af <- fread('immune_sequins.af.tab',header=FALSE)
names(af) <- c('CHROM','AF')

cov <- merge(cov,af,by='CHROM')

cov$PLATFORM <- as.factor(cov$FILE)
#levels(cov$PLATFORM) <-  c("HiSeqXTen/PCR-free","HiSeq2500/PCR-based","HiSeq2500/PCR-based","NextSeq/PCR-based","NextSeq/PCR-based","BGI","ONT")
levels(cov$PLATFORM) <- c("HiSeqXTen/PCR-free","BGI","ONT","HiSeq2500/PCR-based","HiSeq2500/PCR-based","NextSeq/PCR-based","NextSeq/PCR-based")
cov[,MEAN:=MEAN/100]
rsq <- function (x, y) cor(x, y) ^ 2
rmse <- cov[MEAN<=1,.(V1=round(rmse(MEAN,AF),2),V2=signif(rsq(MEAN,AF),2)),by=list(PLATFORM,FILE)]
rmse[,LABEL:=paste("RMSE=",round(V1,3),"; R2=",V2,sep="")]
#[FILE %in% c('NA12878_SEQUIN_HISEQ_X_TEN_PCR_FREE_2','SEQUIN_HISEQ_2500_PCR_BASED_1','NA12878_SEQUIN_MGISEQ2000','NA12878_SEQUIN_PROMETHION')]
ggplot(cov[FILE %in% c('NA12878_SEQUIN_HISEQ_X_TEN_PCR_FREE_2','SEQUIN_HISEQ_2500_PCR_BASED_1','NA12878_SEQUIN_MGISEQ2000','NA12878_SEQUIN_PROMETHION','SEQUIN_NEXTSEQ_500_PCR_BASED_2')],aes(AF,MEAN,color=PLATFORM))+
  geom_point(alpha=0.6,size=2)+
  #geom_smooth()+
  facet_wrap(~PLATFORM,nrow=1)+
  geom_abline(slope=1,color="red",linetype="dashed")+
  geom_text(data=rmse[FILE %in% c('NA12878_SEQUIN_HISEQ_X_TEN_PCR_FREE_2','SEQUIN_HISEQ_2500_PCR_BASED_1','NA12878_SEQUIN_MGISEQ2000','NA12878_SEQUIN_PROMETHION','SEQUIN_NEXTSEQ_500_PCR_BASED_2')],aes(x=0.35,y=0.75,label=LABEL))+
  scale_color_manual(values=c("HiSeqXTen/PCR-free"="#E185B7","HiSeq2500/PCR-based"="#EB8D82","NextSeq/PCR-based"="#00AEEB","BGI"="#ADAB33","ONT"="#00AB8A"))+
  scale_x_continuous(breaks=c(0,0.25,0.5))+
  scale_y_continuous(limits=c(0,1),breaks=c(0,0.5,1))+
  theme(aspect.ratio = 1)+
  guides(color=FALSE)+
  labs(x="Relative frequency",y="Relative coverage")

ggplot(cov[FILE %in% c('NA12878_SEQUIN_HISEQ_X_TEN_PCR_FREE_2','SEQUIN_HISEQ_2500_PCR_BASED_1','NA12878_SEQUIN_MGISEQ2000','NA12878_SEQUIN_PROMETHION','SEQUIN_NEXTSEQ_500_PCR_BASED_2')],aes(AF,MEAN))+
  geom_point()
# Relative Frequency of MAPQ

all.files <- list.files(pattern = ".split_immune.no_edges.sorted.sub.mapq.tab")
l <- lapply(all.files, fread, header=FALSE)
all.files <- str_remove(all.files,".split_immune.no_edges.sorted.sub.mapq.tab")
l <- lapply(seq_along(all.files), function(x,y,i) x[[i]][,FILE:=y[i]],x=l,y=all.files)
mapq <- rbindlist(l)
mapq[,LABEL:="Rearranged"]
mapq[endsWith(V2,"_R"),LABEL:="Background"]
mapq <- mapq[V2!="*"]

mapq_n <- mapq[,.N,by=list(FILE,V4,V2,LABEL)]
mapq_n <- data.table(dcast(mapq_n[V4 %in% c(0,60)],FILE+V2+LABEL~V4,value.var = "N"))
mapq_n <- mapq_n[,.(FILE,V2,LABEL,MAPQ0=`0`,MAPQ60=`60`)]
mapq_n[,RATIO:=MAPQ60/MAPQ0]
mapq_n[,MAPQ60:=MAPQ60/(MAPQ0+MAPQ60)]
ggplot(mapq_n[LABEL=="Rearranged"],aes(FILE,log(RATIO)))+geom_boxplot()

ggplot(mapq_n[LABEL=="Rearranged"],aes(FILE,MAPQ60))+geom_boxplot()

############

#  Reference file

ref <- fread('immune_sequins_ref.tab')
ref <- ref[,.(V1,V3,V5)]
names(ref) <- c("SEQUIN","V3","AF")
ref[grepl("TRD",V3),RECEPTOR:="TRD"]
ref[grepl("TRA",V3),RECEPTOR:="TRA"]
ref[grepl("TRB",V3),RECEPTOR:="TRB"]
ref[grepl("TRG",V3),RECEPTOR:="TRG"]
ref[grepl("IGK",V3),RECEPTOR:="IGK"]
ref[grepl("IGL",V3),RECEPTOR:="IGL"]
vdj <- ref[,tstrsplit(V3,"/")]
vdj[is.na(V4),V4:=V3]
vdj[grepl("J[0-9A-Z]+\\*?[0-9]?",V2),V3:=V2]
vdj[grepl("J[0-9A-Z]+\\*?[0-9]?",V2),V2:=NA]
names(vdj) <- c("V","D","J","CDR3")
ref <- cbind(ref,vdj)

ref[,V:=gsub("\\*[0-9]+","",V)]
ref[,D:=gsub("\\*[0-9]+","",D)]
ref[,J:=gsub("\\*[0-9]+","",J)]

ref[,VDJ:=paste(V,D,J,sep="/")]
ref[,CDR3:=toupper(CDR3)]

all_diagnostic <- list()
index <- 1
#for (file in c('NA12878_SEQUIN_HISEQ_X_TEN_PCR_FREE_2','SEQUIN_HISEQ_2500_PCR_BASED_1','SEQUIN_HISEQ_2500_PCR_BASED_2','SEQUIN_NEXTSEQ_500_PCR_BASED_1','SEQUIN_NEXTSEQ_500_PCR_BASED_2','NA12878_SEQUIN_MGISEQ2000')) {
for (file in c('NA12878_SEQUIN_HISEQ_X_TEN_PCR_FREE_2','SEQUIN_HISEQ_2500_PCR_BASED_1','SEQUIN_HISEQ_2500_PCR_BASED_2','SEQUIN_NEXTSEQ_500_PCR_BASED_1','SEQUIN_NEXTSEQ_500_PCR_BASED_2','NA12878_SEQUIN_MGISEQ2000')) {
  #for (file in c('NA12878_SEQUIN_HISEQ_X_TEN_PCR_FREE_2')) { 
  #file <- 'NA12878_SEQUIN_HISEQ_X_TEN_PCR_FREE_2'
  dt <- fread(sprintf('%s.split_immune.no_edges.sorted.sub.repertoire.clonotypes.summary.txt',file))
  # dt <- cSplit(dt,"allVHitsWithScore",sep=",",direction = "long")
  # dt[!startsWith(allDHitsWithScore,"TR"),allDHitsWithScore:=NA]
  # dt <- cSplit(dt,"allDHitsWithScore",sep=",",direction = "long")
  # dt <- cSplit(dt,"allJHitsWithScore",sep=",",direction = "long")
  dt_vdj <- dt[,.(cloneId,cloneCount,cloneFraction,allVHitsWithScore,allDHitsWithScore,allJHitsWithScore,targetSequences,nSeqCDR3,FILE=file)]
  names(dt_vdj) <- c("ID","COUNT","FRACTION","V","D","J","targetSequences","CDR3","FILE")
  dt_vdj[,CDR3:=substring(CDR3,4,nchar(CDR3)-3)]
  dt_vdj <-merge(dt_vdj,ref[,.(SEQUIN,VDJ,CDR3,AF,RECEPTOR)],by="CDR3",all=TRUE)
  dt_vdj[,FILE:=file]
  all_diagnostic[[index]] <- dt_vdj
  index <- index+1
}
all_diagnostic <- rbindlist(all_diagnostic)
all_diagnostic[!is.na(SEQUIN) & !is.na(ID),STATUS:="TP"]
all_diagnostic[is.na(SEQUIN),STATUS:="FP"]
all_diagnostic[is.na(ID),STATUS:="FN"]

cdr3 <- unique(all_diagnostic[STATUS %in% c("TP","FN"),.(CDR3,FILE,SEQUIN,STATUS,RECEPTOR)])
cdr3[,LABEL:=tstrsplit(SEQUIN,"_",fixed=TRUE)[4]]
cdr3$PLATFORM <- as.factor(cdr3$FILE)
levels(cdr3$PLATFORM) <-  c("HiSeqXTen/PCR-free","HiSeq2500/PCR-based","HiSeq2500/PCR-based","NextSeq/PCR-based","NextSeq/PCR-based","BGI")
cdr3[STATUS=="TP",STATUS:=PLATFORM]

ggplot(cdr3,aes(RECEPTOR,FILE,by=SEQUIN,fill=STATUS))+
  geom_tile(color="lightgray")+
  facet_wrap(~LABEL,scales = "free_x",nrow=1)+
  scale_fill_manual(values=c("HiSeqXTen/PCR-free"="#E185B7","HiSeq2500/PCR-based"="#EB8D82","NextSeq/PCR-based"="#00AEEB","BGI"="#ADAB33","FN"="gray"))+
  theme(axis.text.x = element_text(angle = 45,hjust = 1),
        aspect.ratio = 1,
        axis.text.y=element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank())+
  guides(fill=FALSE)+
  labs(x=NULL,y=NULL)

all_diagnostic[,c("V_Ref","D_Ref","J_Ref"):=tstrsplit(VDJ,"/")]
all_diagnostic[,N:=seq(1,dim(all_diagnostic)[1])]
all_diagnostic[,V_STATUS:=ifelse(grepl(V_Ref,V),TRUE,FALSE),by=N]
all_diagnostic[,D_STATUS:=ifelse(grepl(D_Ref,D),TRUE,FALSE),by=N]
all_diagnostic[,J_STATUS:=ifelse(grepl(J_Ref,J),TRUE,FALSE),by=N]
all_diagnostic[D_Ref=="NA",D_STATUS:=TRUE]
all_diagnostic[V_STATUS==TRUE &  D_STATUS==TRUE & J_STATUS==TRUE,VDJ_STATUS:="TP"]
all_diagnostic[is.na(VDJ_STATUS) & STATUS=="TP",VDJ_STATUS:="FN"]

vdj <- all_diagnostic[STATUS %in% c("TP","FN"),.(VDJ,FILE,SEQUIN,VDJ_STATUS,RECEPTOR)]
vdj <- vdj[,.(VDJ_STATUS=ifelse("TP" %in% VDJ_STATUS,"TP","FN")),by=list(VDJ,FILE,SEQUIN,RECEPTOR)]
vdj[,LABEL:=tstrsplit(SEQUIN,"_",fixed=TRUE)[4]]
vdj$PLATFORM <- as.factor(vdj$FILE)
levels(vdj$PLATFORM) <-  c("HiSeqXTen/PCR-free","HiSeq2500/PCR-based","HiSeq2500/PCR-based","NextSeq/PCR-based","NextSeq/PCR-based","BGI")
vdj[VDJ_STATUS=="TP",VDJ_STATUS:=PLATFORM]

ggplot(vdj,aes(RECEPTOR,FILE,by=SEQUIN,fill=VDJ_STATUS))+
  geom_tile(color="lightgray")+
  facet_wrap(~LABEL,scales = "free_x",nrow=1)+
  scale_fill_manual(values=c("HiSeqXTen/PCR-free"="#E185B7","HiSeq2500/PCR-based"="#EB8D82","NextSeq/PCR-based"="#00AEEB","BGI"="#ADAB33","FN"="gray"))+
  theme(axis.text.x = element_text(angle = 45,hjust = 1),
        aspect.ratio = 1,
        axis.text.y=element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank())+
  guides(fill=FALSE)+
  labs(x=NULL,y=NULL)

precision <- rbind(unique(all_diagnostic[STATUS %in% c("TP"),.(CDR3,FILE,SEQUIN,STATUS,RECEPTOR)]),
                   all_diagnostic[STATUS %in% c("FP"),.(CDR3,FILE,SEQUIN,STATUS,RECEPTOR)])
precision <- data.table(dcast(precision[,.N,by=list(FILE,STATUS)],FILE~STATUS))
precision[,PRECISION:=TP/(TP+FP)]

fraction <- rbind(unique(all_diagnostic[STATUS %in% c("TP"),.(CDR3,FILE,SEQUIN,STATUS,RECEPTOR,FRACTION)]),
      all_diagnostic[STATUS %in% c("FP"),.(CDR3,FILE,SEQUIN,STATUS,RECEPTOR,FRACTION)])
fraction$PLATFORM <- as.factor(fraction$FILE)
levels(fraction$PLATFORM) <-  c("HiSeqXTen/PCR-free","HiSeq2500/PCR-based","HiSeq2500/PCR-based","NextSeq/PCR-based","NextSeq/PCR-based","BGI")

ggplot(fraction,aes(FILE,FRACTION,color=PLATFORM))+
  geom_point(size=2)+
  scale_color_manual(values=c("HiSeqXTen/PCR-free"="#E185B7","HiSeq2500/PCR-based"="#EB8D82","NextSeq/PCR-based"="#00AEEB","BGI"="#FFCC04","FN"="gray"))+
  facet_wrap(~STATUS)+
  guides(color=FALSE)+
  theme(aspect.ratio = 1)+
  labs(x=NULL,y="Observed frequency")

#all_diagnostic[,CDR3:=substring(CDR3,4,nchar(CDR3)-3)]

#all_diagnostic <- merge(all_diagnostic,ref[,.(SEQUIN,VDJ,CDR3,AF)],by="CDR3",all=TRUE)
dev.off()

