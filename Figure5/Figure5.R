setwd('/simon/ANDRE/Figures/github/Figure5')

pdf('plots_figure4_v4.pdf',height=8,width=13,useDingbats = FALSE)

library(stringr)
library(gridExtra)
library(MESS)

#####################

ref <- fread('zcat sequin_smallvariants_hg38_2.7.hap.vcf.gz',skip="#CHROM")
ref[,c("AF","AL","GT"):=tstrsplit(INFO,";",fixed=TRUE)]
ref[,AF:=str_remove(AF,"AF=")]
ref[,AL:=str_remove(AL,"AL=")]
ref[,GT:=str_remove(GT,"GT=")]
ref <- ref[,.(CHROM=`#CHROM`,POS,ID,REF,ALT,AF,AL,GT)]

dt <- list()
j <- 1

for (i in c('NA12878_SEQUIN_HISEQ_X_TEN_PCR_FREE_2','SEQUIN_HISEQ_2500_PCR_BASED_1','SEQUIN_HISEQ_2500_PCR_BASED_2','SEQUIN_NEXTSEQ_500_PCR_BASED_1','SEQUIN_NEXTSEQ_500_PCR_BASED_2','NA12878_SEQUIN_MGISEQ2000','NA12878_SEQUIN_PROMETHION')) {
  file <- sprintf('%s.sequin.sorted.phased',i)
  tp <- fread(sprintf('zcat %s/tp.vcf.gz',file),skip="#CHROM")
  if (nrow(tp)!=0) {
    tp[,c("GT_CALL"):=tstrsplit(s1,":",fixed=TRUE)[1]]
    tp[,HAP_GROUP:=NA]
    if (i != "NA12878_SEQUIN_PROMETHION") {
      tp[,c("HAP_GROUP"):=tstrsplit(s1,":",fixed=TRUE)[6]]
    } else {
      tp[,c("HAP_GROUP"):=tstrsplit(s1,":",fixed=TRUE)[5]]
    }
    tp <- tp[,.(CHROM=`#CHROM`,POS,REF,ALT,GT_CALL,HAP_GROUP,STATUS="TP")]
  } else {
    tp <- tp[,.(CHROM=character(),POS=numeric(),REF=character(),ALT=character(),GT_CALL=character(),HAP_GROUP=numeric(),STATUS=character())]
  }
  fp <- fread(sprintf('zcat %s/fp.vcf.gz',file),skip="#CHROM")
  if (nrow(fp)!=0) {
    fp[,c("GT_CALL"):=tstrsplit(s1,":",fixed=TRUE)[1]]
    fp[,HAP_GROUP:=NA]
    #fp[,c("HAP_GROUP"):=tstrsplit(s1,":",fixed=TRUE)[6]]
    fp <- fp[,.(CHROM=`#CHROM`,POS,REF,ALT,GT_CALL,HAP_GROUP,STATUS="FP")]
  } else {
    fp <- fp[,.(CHROM=character(),POS=numeric(),REF=character(),ALT=character(),GT_CALL=character(),HAP_GROUP=numeric(),STATUS=character())]
  }
  
  all <- rbind(tp,fp)
  all[,ALT:=NULL]
  
  all <- merge(ref,all,by=c("CHROM","POS","REF"),all.x = TRUE,all.y = TRUE)
  all[,ALLELE1:=0]
  all[,ALLELE2:=0]
  all[AL=="1" | AL=="1,2",ALLELE1:=1]
  all[AL=="2" | AL=="1,2",ALLELE2:=1]
  all[,GT1:=paste(ALLELE1,ALLELE2,sep="|")]
  all[,GT2:=paste(ALLELE2,ALLELE1,sep="|")]
  all[GT_CALL==GT1,ST1:=TRUE]
  all[GT_CALL==GT2,ST2:=TRUE]
  all[,PHASED:="NO"]
  all[!is.na(HAP_GROUP),PHASED:="YES"]
  all[,CORRECTLY:="NO"]
  a2 <- all[!is.na(HAP_GROUP),.(ifelse(all(ST1) | all(ST2),TRUE,FALSE)),by=HAP_GROUP]
  all[HAP_GROUP %in% a2[V1==TRUE,HAP_GROUP],CORRECTLY:="YES"]
  all[,c("V1","V2","V3"):=tstrsplit(ID,"_")[1:3]]
  all[,SEQUIN:=paste(V1,V2,V3,sep="_")]
  all[,c("V1","V2","V3"):=NULL]
  all <- all[order(SEQUIN,CHROM,POS)]
  all[,DISTANCE:=as.numeric()]
  all[,DISTANCE:=c(diff(POS),-1),by=SEQUIN]
  all[is.na(ID)| DISTANCE<0,DISTANCE:=NA]
  all[,FILE:=file]
  all[is.na(STATUS),STATUS:="FN"]
  dt[[j]] <- all
  j <- j+1
}

dt <- rbindlist(dt)
dt[,FILE:=str_remove(FILE,".sequin.sorted.phased")]
dt$FILE <- as.factor(dt$FILE)

x <- table(dt[,FILE],dt[,STATUS])
stats <- dcast(data.table(x),V1~V2,value.var = "N")
stats[,FP_Proportion:=FP/(FN+FP+TP)]
stats[,READ_LENGTH:="short-read"]
stats[V1=="NA12878_SEQUIN_PROMETHION",READ_LENGTH:="long-read"]
stats[,.(mean(FP_Proportion),sd(FP_Proportion)),by=READ_LENGTH]

x <- melt(x,id.var="FILE",variable.name="S4TATUS",value.name = "VALUE")

ggplot(x,aes(Var1,VALUE,fill=Var2))+
  geom_bar(stat="identity")+
  theme(aspect.ratio = 1)+
  labs(y="Frequency")

####################################
# Allele - BED files

for (i in c('NA12878_SEQUIN_HISEQ_X_TEN_PCR_FREE_2','SEQUIN_HISEQ_2500_PCR_BASED_1','SEQUIN_HISEQ_2500_PCR_BASED_2','SEQUIN_NEXTSEQ_500_PCR_BASED_1','SEQUIN_NEXTSEQ_500_PCR_BASED_2','NA12878_SEQUIN_MGISEQ2000','NA12878_SEQUIN_PROMETHION')) {

  file <- sprintf('%s',i)
  tmp <- dt[FILE==file]
  hg <- tmp[!is.na(HAP_GROUP),.(tSTART=min(POS)-1,tEND=max(POS)-1),by=list(FILE,HAP_GROUP,CHROM)]
  hg <- hg[,.(CHROM,tSTART,tEND,HAP_GROUP)]
  #tmp <- merge(tmp,hg,by=c('FILE','HAP_GROUP'),all.x=TRUE)
  tmp[,tmpALLELE1:=ALLELE1]
  tmp[,tmpALLELE2:=ALLELE2]
  tmp[PHASED=="YES" & CORRECTLY=="NO",ALLELE1:=tmpALLELE2]
  tmp[PHASED=="YES" & CORRECTLY=="NO",ALLELE2:=tmpALLELE1]
  al1 <- tmp[ALLELE1==0,COLOR1:='255,0,0']
  al1 <- tmp[ALLELE1!=0,COLOR1:='0,255,0']
  al2 <- tmp[ALLELE2==0,COLOR2:='255,0,0']
  al2 <- tmp[ALLELE2!=0,COLOR2:='0,255,0']
  al1 <- al1[,.(CHROM,POS-1,END=POS,V1='.',V2=0,V3="+",V4=POS-1,V5=POS,V6=COLOR1)]
  al2 <- al2[,.(CHROM,POS-1,END=POS,V1='.',V2=0,V3="+",V4=POS-1,V5=POS,V6=COLOR2)]
  write.table(al1,file=sprintf('alleles/%s.allele1.bed',file),col.names = FALSE,quote=FALSE,row.names = FALSE,sep="\t")
  write.table(al2,file=sprintf('alleles/%s.allele2.bed',file),col.names = FALSE,quote=FALSE,row.names = FALSE,sep="\t")
  write.table(hg,file=sprintf('alleles/%s.haplotypes.bed',file),col.names = FALSE,quote=FALSE,row.names = FALSE,sep="\t")
}

####################################

example <- 'S0333_HP_009'
#example <- 'S0335_HP_010'
#example <- 'S0322_HP_003'

hap_ref <- data.table(ref)
hap_ref[,ID:=str_sub(ID, end=-4)]
hap_ref[,ALLELE1:=0]
hap_ref[,ALLELE2:=0]
hap_ref[AL=="1" | AL=="1,2",ALLELE1:=1]
hap_ref[AL=="2" | AL=="1,2",ALLELE2:=1]
hap_ref <- hap_ref[,.(POS,FILE="REF",STATUS="TP",ALLELE1,ALLELE2,SEQUIN=ID,HAP_GROUP=NA,PHASED="YES")]

#######################################
#######################################

# How many HET variants per sequin

het <- hap_ref[!(ALLELE1==1 & ALLELE2==1)]
het_n <- het[,.(N=.N),by=SEQUIN]

het[,DISTANCE:=as.numeric()]
het <- het[order(SEQUIN,POS)]
het[,DISTANCE:=c(diff(POS),-1),by=SEQUIN]
het[DISTANCE<0,DISTANCE:=NA]
het <- het[order(DISTANCE)]
het[,CUMSUM:=cumsum(DISTANCE)]
het[,PAIR:=c(POS[-1],NA),by=SEQUIN]
het[,ID:=paste(SEQUIN,POS,PAIR,sep=":")]

ggplot(het,aes(FILE,DISTANCE))+
  #geom_violin(width=0.4,draw_quantiles = c(0.5))+
  geom_boxplot(width=0.2,outlier.shape = NA)+
  geom_jitter(width=0.1,alpha=0.4)+
  theme(aspect.ratio = 1,axis.text.x = element_blank(),axis.ticks.x = element_blank())+
  geom_hline(yintercept = 217,linetype="dashed")+
  labs(x=NULL,y="Distance")

# ggplot(het,aes(DISTANCE,CUMSUM/max(het$CUMSUM,na.rm = TRUE)))+
#   geom_line()+
#   theme(aspect.ratio = 1)+
#   labs(y="Cumulative distribution",x="Distance")+
#   geom_vline(xintercept = 217,linetype="dashed")

#######################################



hap <- dt[SEQUIN==example]
hap[,tmpALLELE1:=ALLELE1]
hap[,tmpALLELE2:=ALLELE2]
hap[PHASED=="YES" & CORRECTLY=="NO",ALLELE1:=tmpALLELE2]
hap[PHASED=="YES" & CORRECTLY=="NO",ALLELE2:=tmpALLELE1]
hap_ref <- hap_ref[SEQUIN==example]
hap <- hap[,.(POS,FILE,STATUS,ALLELE1,ALLELE2,SEQUIN,HAP_GROUP,PHASED)]
hap <- rbind(hap,hap_ref)
hap <- melt(hap,id.vars = c("POS","FILE","STATUS","SEQUIN",'HAP_GROUP','PHASED'),variable.name = "ALLELE",value.name = "VALUE")
hap[,LABEL:=paste(FILE,ALLELE,sep=".")]

pair <- dt[FILE!="REF",.(POS,FILE,SEQUIN,PHASED,HAP_GROUP)]
pair <- unique(pair[order(FILE,SEQUIN,POS)])
pair <- pair[PHASED=="YES"]
pair[,PAIR:=c(POS[-1],NA),by=list(FILE,SEQUIN,HAP_GROUP)]
pair <- pair[!is.na(PAIR)]
pair[,DISTANCE:=PAIR-POS]
pair$PLATFORM <- as.factor(as.character(pair$FILE))
levels(pair$PLATFORM) <- c("HiSeqXTen/PCR-free","BGI","ONT","HiSeq2500/PCR-based","HiSeq2500/PCR-based","NextSeq/PCR-based","NextSeq/PCR-based")
labels <- pair[,.(Q1=quantile(DISTANCE,0.25),Q2=quantile(DISTANCE,0.5),Q3=quantile(DISTANCE,0.75)),by=PLATFORM]
labels[,LABEL:=paste("Q1=",format(Q1,2),"; Q2=",format(Q2,2),"; Q3=",format(Q3,2),sep="")]
labels[,X:=2500]
labels[,Y:=seq(0.6,0.76,0.04)]

ggplot(pair,aes(PLATFORM,DISTANCE,color=PLATFORM))+
  geom_boxplot()+
  labs(x=NULL,y="Pairwise distance")+
  scale_color_manual(values=c("HiSeqXTen/PCR-free"="#E185B7","HiSeq2500/PCR-based"="#EB8D82","NextSeq/PCR-based"="#00AEEB","BGI"="#ADAB33","ONT"="#00AB8A"))+
  theme(aspect.ratio = 1)

ggplot(pair,aes(DISTANCE,color=PLATFORM))+
  geom_density(aes(y=..scaled..))+
  labs(x=NULL,y="Pairwise distance")+
  scale_color_manual(values=c("HiSeqXTen/PCR-free"="#E185B7","HiSeq2500/PCR-based"="#EB8D82","NextSeq/PCR-based"="#00AEEB","BGI"="#ADAB33","ONT"="#00AB8A"))+
  theme(aspect.ratio = 1)+
  scale_y_continuous(limits=c(0,1),breaks=c(0,0.5,1))+
  geom_text(data=labels,aes(x=X,y=Y,label=LABEL))+
  geom_vline(xintercept=c(300,1000),linetype="dashed")+
  guides(color=FALSE)


hg <- hap[!is.na(HAP_GROUP),.(tSTART=min(POS)-1,tEND=max(POS)-1),by=list(FILE,HAP_GROUP,SEQUIN)]
hg_ref <- hap[FILE=="REF",.(tSTART=min(POS)-1,tEND=max(POS)-1),by=list(FILE,HAP_GROUP,SEQUIN)]
hg <- rbind(hg,hg_ref)
hg <- rbind(hg[,.(LABEL=paste(FILE,"ALLELE1",sep="."),FILE,SEQUIN,tSTART,tEND,HAP_GROUP)],hg[,.(LABEL=paste(FILE,"ALLELE2",sep="."),FILE,SEQUIN,tSTART,tEND,HAP_GROUP)])

hap$FILE <- factor(hap$FILE,levels=c("REF","L.V.158","NA12878_SEQUIN_HISEQ_X_TEN_PCR_FREE_2","SEQUIN_HISEQ_2500_PCR_BASED_1","SEQUIN_HISEQ_2500_PCR_BASED_2","SEQUIN_NEXTSEQ_500_PCR_BASED_1","SEQUIN_NEXTSEQ_500_PCR_BASED_2","NA12878_SEQUIN_MGISEQ2000","NA12878_SEQUIN_PROMETHION"))

sub <- c('REF','NA12878_SEQUIN_HISEQ_X_TEN_PCR_FREE_2','SEQUIN_HISEQ_2500_PCR_BASED_1','SEQUIN_NEXTSEQ_500_PCR_BASED_1','NA12878_SEQUIN_MGISEQ2000','NA12878_SEQUIN_PROMETHION')

ggplot(hap[STATUS!="FN" & FILE %in%  sub],aes(POS,LABEL))+
  geom_segment(data=hg[FILE %in%  sub],aes(x=tSTART,y=LABEL,xend=tEND,yend=LABEL),color="gray",size=8,alpha=0.4)+
  geom_point(aes(color=as.factor(VALUE)),pch="|",size=8,alpha=0.3)+
  geom_point(data=hap[STATUS!="FN" & PHASED=="YES" & FILE %in%  sub],aes(color=as.factor(VALUE)),pch="|",size=8)+
  facet_grid(FILE~.,scale="free_y",space="free_y",switch = "y")+
  labs(x="Position (bp)",y=NULL,color=NULL)+
  theme(axis.text.y=element_blank(),axis.ticks.y=element_blank())+
  guides(color=FALSE)

####################################
# Insert size

short.is <- fread('all_insert_size_metrics.tab')
short.is[,V4:=seq(1,.N)]
short.is <- short.is[,.(FILE=V1,SIZE=rep(V2,V3)),by=V4]
short.is[,V4:=NULL]
long.rl <- fread('NA12878_SEQUIN_PROMETHION.sequin.sorted.read.lenght.txt')
insert.size <- rbind(short.is,long.rl[,.(FILE='NA12878_SEQUIN_PROMETHION',SIZE=V1)])
insert.size[,PLATFORM:=FILE]
insert.size$PLATFORM <- as.factor(insert.size$PLATFORM )
levels(insert.size$PLATFORM) <- c("HiSeqXTen/PCR-free","HiSeqXTen/PCR-free","HiSeq2500/PCR-based","HiSeq2500/PCR-based","NextSeq/PCR-based","NextSeq/PCR-based","BGI","ONT")

labels <- insert.size[FILE!="NA12878_SEQUIN_PROMETHION",.(round(mean(SIZE),2),round(sd(SIZE),2),V3=round(median(SIZE),2)),by=PLATFORM]
labels[,LABEL:=paste("Mean=",V1,"; SD=",V2,sep="")]
labels[,X:=600]
labels[,Y:=seq(0.6,0.72,0.04)]

ggplot(insert.size[FILE!="NA12878_SEQUIN_PROMETHION"],aes(SIZE,color=PLATFORM))+
  #geom_bar(stat="identity",color="black")+
  geom_density(aes(y=..scaled..))+
  #geom_density()+
  scale_x_continuous(limits=c(0,800),breaks=c(0,400,800))+
  scale_y_continuous(breaks=c(0,0.5,1))+
  scale_color_manual(values=c("HiSeqXTen/PCR-free"="#E185B7","HiSeq2500/PCR-based"="#EB8D82","NextSeq/PCR-based"="#00AEEB","BGI"="#ADAB33","ONT"="#00AB8A"))+
  theme(aspect.ratio = 1)+
  geom_text(data=labels,aes(x=X,y=Y,label=LABEL))+
  geom_vline(xintercept = 289.63,linetype="dashed")+
  labs(x="Fragment length",y="Density")+
  guides(color=FALSE)
  #labs(x=NULL,y="Average fragment length")

test <- insert.size[,.N,by=list(PLATFORM,SIZE)]
test <- merge(test,test[,max(N),by=PLATFORM],by="PLATFORM")
test[,N:=N/V1]
auc <- test[,round(MESS::auc(SIZE,N),2),by=PLATFORM]
subset <- test[SIZE>289.63 & SIZE<800,round(MESS::auc(SIZE,N),2),by=PLATFORM]

# fraction/subset different beetween hiseq x ten and bgi 
(subset[PLATFORM=="HiSeqXTen/PCR-free",V1]-subset[PLATFORM=="BGI",V1])/auc[PLATFORM=="HiSeqXTen/PCR-free",V1]

ggplot(test[PLATFORM!="ONT"],aes(SIZE,N,color=PLATFORM))+
  geom_line()+
  geom_vline(xintercept = 289.63,linetype="dashed")+
  scale_x_continuous(limits=c(0,800),breaks=c(0,400,800))

localMaxima <- function(x) {
  # Use -Inf instead if x is numeric (non-integer)
  y <- diff(c(-.Machine$integer.max, x)) > 0L
  rle(y)$lengths
  y <- cumsum(rle(y)$lengths)
  y <- y[seq.int(1L, length(y), 2L)]
  if (x[[1]] == x[[2]]) {
    y <- y[-1]
  }
  y
}

ggplot(insert.size[FILE=="NA12878_SEQUIN_PROMETHION"],aes(SIZE,by=FILE,color=PLATFORM))+
  #geom_bar(stat="identity",color="black")+
  geom_density(aes(y=..scaled..))+
  scale_x_continuous(limits=c(0,4000))+
  scale_y_continuous(breaks=c(0,0.5,1))+
  geom_vline(xintercept = 510,linetype="dashed")+
  geom_vline(xintercept = 1780,linetype="dashed")+
  scale_color_manual(values=c("HiSeqXTen/PCR-free"="#E185B7","HiSeq2500/PCR-based"="#EB8D82","NextSeq/PCR-based"="#00AEEB","BGI"="#ADAB33","ONT"="#00AB8A"))+
  theme(aspect.ratio = 1)+
  labs(x="Fragment length",y="Density")+
  guides(color=FALSE)
#labs(x=NULL,y="Average fragment length")

####################################

dt[,PLATFORM:=FILE]
dt1 <- dt[GT=="HET"]
#levels(dt1$PLATFORM) <- c("HiSeqXTen/PCR-free","HiSeqXTen/PCR-free","HiSeq2500/PCR-based","HiSeq2500/PCR-based","NextSeq/PCR-based","NextSeq/PCR-based","BGI","ONT")
levels(dt1$PLATFORM) <- c("HiSeqXTen/PCR-free","BGI","ONT","HiSeq2500/PCR-based","HiSeq2500/PCR-based","NextSeq/PCR-based","NextSeq/PCR-based")

####################################
remove <- c('NA_NA_NA','S0362_HP_024')
hp1 <- dt1[,.(SIZE=max(POS)-min(POS),.N),by=list(FILE,PLATFORM,SEQUIN)]
#dcast(hp1[!(SEQUIN %in% remove),.(FILE,SEQUIN,N)],SEQUIN~FILE,value.var = "N")
hp2 <- dt1[!is.na(HAP_GROUP) & PHASED=="YES",.(BLOCK_SIZE=max(POS)-min(POS),PHASED=.N),by=list(FILE,PLATFORM,SEQUIN,HAP_GROUP)]
#hp2 <- dt[!is.na(HAP_GROUP) & PHASED=="YES",.(BLOCK_SIZE=max(POS)-min(POS),PHASED=.N),by=list(FILE,PLATFORM,SEQUIN)]
hp1 <- merge(hp1,hp2,by=c("FILE","PLATFORM","SEQUIN"))
hp2 <- dt1[!is.na(HAP_GROUP) & CORRECTLY=="YES",.(CORRECTLY=.N),by=list(FILE,PLATFORM,SEQUIN,HAP_GROUP)]
hp1 <- merge(hp1,hp2,by=c("FILE","PLATFORM","SEQUIN","HAP_GROUP"))
investigate <- hp1[BLOCK_SIZE==0]
hp1 <- hp1[!BLOCK_SIZE==0]

ggplot(hp1,aes(FILE,log10(BLOCK_SIZE),fill=PLATFORM))+
  geom_violin(draw_quantiles = c(0.5))+
  labs(y="Haplotype block size (Log10)",x=NULL,fill=NULL)+
  #scale_y_continuous(limits=c(0,4),breaks=c(0,2,4))+
  scale_fill_manual(values=c("HiSeqXTen/PCR-free"="#E185B7","HiSeq2500/PCR-based"="#EB8D82","NextSeq/PCR-based"="#00AEEB","BGI"="#ADAB33","ONT"="#00AB8A"))+
  theme(aspect.ratio=1,legend.position = "bottom",axis.text.x=element_blank(),axis.ticks.x = element_blank())+
  guides(fill=FALSE)

#wilcox.test(hp1[PLATFORM =="ONT",BLOCK_SIZE],hp1[PLATFORM =="HiSeqXTen/PCR-free",BLOCK_SIZE])
pvalue <- pairwise.wilcox.test(hp1[,log10(BLOCK_SIZE)],hp1[,PLATFORM],p.adjust.method="BH")$p.value
pvalue <- data.table(melt(pvalue))
pvalue[,Significance:=cut(value,breaks=c(-Inf,0.0001,0.001,0.01,0.05,1),labels=c("****","***","**","*","n.s."))]

ggplot(pvalue,aes(Var1,Var2,fill=value))+
  geom_tile()+
  geom_text(aes(label=Significance),color="white")+
  theme(aspect.ratio = 1)+
  labs(x=NULL,y=NULL,fill="P-value")+
  scale_fill_gradient(limits=c(0,1),breaks=c(0,0.5,1))

hp <- data.table(hp1)
hp[,PHASED:=PHASED/N*100]
hp[,CORRECTLY:=CORRECTLY/N*100]
hp[,BLOCK_SIZE:=BLOCK_SIZE/SIZE*100]

ggplot(hp,aes(FILE,BLOCK_SIZE))+geom_point(alpha=0.6)
ggplot(hp,aes(FILE,PHASED))+geom_point(alpha=0.6)

hp_stats <- hp[,.(BLOCK_SIZE=mean(BLOCK_SIZE),PHASED=mean(PHASED)),by=list(FILE,PLATFORM)]

ggplot(hp,aes(PHASED,BLOCK_SIZE,by=FILE,color=PLATFORM))+
  geom_point(alpha=0.6,size=3)+
  scale_x_continuous(limits=c(0,100),breaks=c(0,50,100))+
  scale_y_continuous(limits=c(0,100),breaks=c(0,50,100))+
  #facet_wrap(~FILE,nrow=1)+
  labs(x="Phased variants (%)",y="Haplotype block size (%)")+
  scale_color_manual(values=c("HiSeqXTen/PCR-free"="#E185B7","HiSeq2500/PCR-based"="#EB8D82","NextSeq/PCR-based"="#00AEEB","BGI"="#ADAB33","ONT"="#00AB8A"))+
  theme(aspect.ratio = 1)+
  guides(color=FALSE)

####################################

# 1) Fraction of variants correctly identified?

stats <- data.table(dcast(dt[,.N,by=list(FILE,STATUS)], FILE~STATUS,value.var = "N"))
stats[is.na(stats)] <- 0
stats[,TOTAL:=(FN+FP+TP)]
stats[,SN:=TP/TOTAL]
stats[,SP:=TP/(TP+FP)]
stats
ss <- tableGrob(stats)
#grid.arrange(ss,top="Fraction of variants correctly identified?")

# 2) Fraction of variants phased? Correctly phased?

phased <- data.table(dcast(dt[,.N,by=list(FILE,PHASED)], FILE~PHASED,value.var = "N"))
phased[,TOTAL:=(YES+NO)]
phased[,SN:=YES/TOTAL]

correctly <- data.table(dcast(dt[,.N,by=list(FILE,CORRECTLY)], FILE~CORRECTLY,value.var = "N"))
correctly <- correctly[,.(FILE,CORRECTLY=YES)]
phased <- merge(phased,correctly,by="FILE")
phased[,SP:=CORRECTLY/YES]
phased
x <- stats[,.(Sample=FILE,Sensitivity=round(SN,2),Specificity=round(SP,2))]
y <- phased[,.(Phased=round(SN,2),Precision=round(SP,2))]
phased <- cbind(x,y)
ss <- tableGrob(phased)
grid.arrange(ss,top="Fraction of variants phased? Correctly phased?")

# 3) size of haplotype blocks
# hp <- dt[!is.na(HAP_GROUP)]
# hp1 <- hp[,.(SIZE=max(POS)-min(POS),.N),by=list(FILE,PLATFORM,HAP_GROUP)]
# hp2hp[PHASED=="YES",.(PHASED=.N),by=list(FILE,PLATFORM,HAP_GROUP)]
# 
# ggplot(hp,aes(FILE,log10(SIZE),fill=PLATFORM))+
#   geom_boxplot()+
#   labs(y="Haplotype block size (bp)",x=NULL,fill=NULL)+
#   scale_y_continuous(limits=c(0,4),breaks=c(0,2,4))+
#   scale_fill_manual(values=c("HiSeqXTen/PCR-free"="#E185B7","HiSeq2500/PCR-based"="#EB8D82","NextSeq/PCR-based"="#00AEEB","BGI"="#ADAB33","ONT"="#00AB8A"))+
#   theme(aspect.ratio=1,legend.position = "bottom",axis.text.x=element_blank(),axis.ticks.x = element_blank())

# 4) Probability of phasing relative to distance between
dist <- dt[GT=="HET"]
#dist[,DISTANCE:=as.numeric()]
dist <- dist[order(FILE,PLATFORM,SEQUIN,POS)]
dist[,DISTANCE:=c(diff(POS),-1),by=list(FILE,PLATFORM,SEQUIN)]
dist[is.na(ID) | DISTANCE<0,DISTANCE:=NA]
dist <- dist[!is.na(DISTANCE)]
dist[,.(PHASED,DISTANCE,FILE,PLATFORM)]
dist[PHASED=="YES",PHASED:="1"]
dist[PHASED=="NO",PHASED:="0"]
dist[,PHASED:=as.numeric(PHASED)]
dist <- dist[order(FILE,DISTANCE,decreasing = TRUE)]
dist <- dist[!is.na(DISTANCE),.(DISTANCE,PHASED=cumsum(PHASED)),by=list(FILE,PLATFORM)]
dist <- merge(dist,dist[,.(TOTAL=max(PHASED)),by=FILE],by="FILE")
dist[,FRACTION:=PHASED/TOTAL*100]
#levels(dist$PLATFORM) <- c("HiSeqXTen/PCR-free","HiSeqXTen/PCR-free","HiSeq2500/PCR-based","HiSeq2500/PCR-based","NextSeq/PCR-based","NextSeq/PCR-based","BGI","ONT")
levels(dist$PLATFORM) <- c("HiSeqXTen/PCR-free","BGI","ONT","HiSeq2500/PCR-based","HiSeq2500/PCR-based","NextSeq/PCR-based","NextSeq/PCR-based")

ggplot(dist,aes(DISTANCE,FRACTION,by=FILE,color=PLATFORM))+
  geom_line()+
  theme(aspect.ratio = 1,legend.position = "bottom")+
  scale_x_continuous(limits=c(0,1000),breaks=c(0,500,1000))+
  scale_y_continuous(limits=c(0,100),breaks=c(0,50,100))+
  scale_color_manual(values=c("HiSeqXTen/PCR-free"="#E185B7","HiSeq2500/PCR-based"="#EB8D82","NextSeq/PCR-based"="#00AEEB","BGI"="#ADAB33","ONT"="#00AB8A"))+
  labs(x="Pairwise SNV distance",y="SNVs phased (%)",color=NULL)+
  geom_vline(xintercept = 217,linetype="dashed")+
  #geom_vline(xintercept = 700,linetype="dashed")+
  guides(color=FALSE)

auc1 <- dist[,auc(DISTANCE,FRACTION),by=list(FILE,PLATFORM)]
auc1[,V1:=V1/(100*1000)]


dev.off()