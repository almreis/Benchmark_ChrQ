setwd('/simon/ANDRE/Figures/github/Figure2')

# R packages

library(cowplot)
library(ggpubr)
library(quantmod)
library(changepoint)
library(data.table)
library(stringr)
library(ggplot2)
library(reshape2)

theme_set(theme_cowplot(font_size=20))

pdf('Figure2.pdf')

# Upload pysamstats statistics

all.files <- list.files(pattern = ".all_metrics.tab")
l <- lapply(all.files, fread, header=TRUE)
all.files <- str_remove(all.files,".sequin.sorted.full.all_metrics.tab")
l <- lapply(seq_along(all.files), function(x,y,i) x[[i]][,FILE:=y[i]],x=l,y=all.files)
dt <- rbindlist(l)
dt[,PLATFORM:=as.factor(FILE)]
levels(dt$PLATFORM) <- c("HiSeqXTen/PCR-free","HiSeqXTen/PCR-free","BGI","ONT","HiSeq2500/PCR-based","HiSeq2500/PCR-based","NextSeq/PCR-based","NextSeq/PCR-based")


# File containing 100bp sliding window to calculate GC content

gc_context <- fread('sequin_smallvariants_hg38_2.6.gc_window.tab',header=FALSE)
names(gc_context) <- c('CHROM','POS','NAME','SEQ')
gc_context[,A:= str_count(SEQ, "A")]
gc_context[,C:= str_count(SEQ, "C")]
gc_context[,G:= str_count(SEQ, "G")]
gc_context[,T:= str_count(SEQ, "T")]
gc_context[,SUM:=A+C+G+T]
gc_context <- gc_context[,.(CHROM,POS=as.numeric(POS),GC_CONTEXT=(C+G)/SUM*100)]

# Merge per base summary statistics with GC contennt information

dt <- merge(dt,gc_context,by=c("CHROM","POS"))

dt[,ID:=paste(CHROM,POS,sep=':')]

# Calculate a normalization factor based on the total sum of read counts for each library

norm_factor <- dt[,sum(READS_ALL),by=FILE][,.(NORM_FACTOR=V1/10^6),by=FILE]

dt <- merge(dt,norm_factor,by="FILE")

# Distribution of per base normalized coverage per library 

p1 <- ggplot(dt,aes(FILE,READS_ALL/NORM_FACTOR,fill=as.factor(PLATFORM)))+
  geom_violin(trim=FALSE, alpha = 0.6,draw_quantiles = 0.5,adjust=2,width=0.5)+
  coord_flip()+
  scale_y_continuous(limits=c(0,4),breaks=c(0,2,4),position="right")+
  labs(y=NULL,x=NULL)+
  guides(fill=FALSE)+
  theme(aspect.ratio = 1,
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

dt[,.(mean(READS_ALL/NORM_FACTOR),sd(READS_ALL/NORM_FACTOR)),by=PLATFORM]

vcf <- fread('zcat sequin_smallvariants_hg38_2.7.vcf.gz',skip="#CHROM")
vcf[,POS:=POS-1]
vcf[,N:=nchar(REF)]
vcf[N==1,END:=POS]
vcf[N!=1,END:=POS+(N-1)]
vcf <- vcf[,c(1,2,10,3)]
names(vcf) <- c('CHROM','START','END','NAME')
vcf <- vcf[,.(CHROM,POS=seq(START,END)),by=NAME]
vcf[,ID:=paste(CHROM,POS,sep=':')]

error <- dt[!(ID %in% vcf$ID),]
error <- error[,.(CHROM,NAME,POS,ID,FILE,PLATFORM,READS_ALL,MISMATCHES,DELETIONS,INSERTIONS)]
window <- 1000
windows <- unique(error[,.(CHROM,POS,NAME,ID)])
windows <- windows[order(CHROM,POS,NAME)]
windows[,WINDOW:=paste(NAME,rep(1:round(.N/window+0.5),each=window),sep="_"),by=NAME]
windows <- windows[,.(ID,WINDOW)]
error <- merge(error,windows,by="ID",allow.cartesian = TRUE)
error_long <- melt(error,measure.vars = c("MISMATCHES","DELETIONS","INSERTIONS"),variable.name = "ERROR_TYPE",value.name="VALUE")
error_summary <- error_long[,.(VALUE=sum(VALUE),READS_ALL=sum(READS_ALL)),by=list(FILE,PLATFORM,ERROR_TYPE,WINDOW)]
error_summary[,FRACTION:=VALUE/READS_ALL]

error_summary[ERROR_TYPE=="MISMATCHES",.(mean(FRACTION),sd(FRACTION)),by=list(PLATFORM)]

# Distribution of frequency of mismatches per library 

p2 <- ggplot(error_summary[ERROR_TYPE=="MISMATCHES"],aes(FILE,log2(FRACTION),fill=PLATFORM))+
  geom_violin(trim=FALSE,draw_quantiles = 0.5,scale="area",alpha=0.6,width=0.5)+
  labs(y=NULL,x=NULL)+
  guides(fill=FALSE)+
  coord_flip()+
  scale_y_continuous(breaks=c(-5,-10),position="right")+
  theme(aspect.ratio = 1,
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

errors <- data.table(melt(dt[,.(FILE,PLATFORM,ID,REF.x,A,C,T,G)],id.vars = c("FILE","PLATFORM","ID","REF.x"),variable.name = "BASE", value.name = "COUNT"))
errors <- errors[COUNT>0 & REF.x != BASE]
errors[,MISMATCH:=paste(REF.x,BASE,sep=">")]
errors[MISMATCH=="G>T",MISMATCH:="C>A"]
errors[MISMATCH=="G>C",MISMATCH:="C>G"]
errors[MISMATCH=="G>A",MISMATCH:="C>T"]
errors[MISMATCH=="A>T",MISMATCH:="T>A"]
errors[MISMATCH=="A>G",MISMATCH:="T>C"]
errors[MISMATCH=="A>C",MISMATCH:="T>G"]
errors[,TYPE:="Transversion"]
errors[MISMATCH %in% c("C>T","T>C"),TYPE:="Transition"]
expand <- errors[rep(1:.N,COUNT)]
expand <- expand[,.N,by=list(FILE,PLATFORM,TYPE,MISMATCH)]
total <- expand[,.(TOTAL=sum(N)),by=list(FILE,PLATFORM)]
expand <- merge(expand,total,by=c("PLATFORM","FILE"))
expand[,N:=N/TOTAL*100]

ggplot(expand,aes(MISMATCH,N,fill=PLATFORM,by=FILE))+
  geom_bar(stat="identity",position=position_dodge(),alpha=0.6,color="black")+
  #geom_point(alpha=0.6,size=3)+
  labs(y="Relative frequency (%)",x=NULL)+
  guides(fill=FALSE)+
  theme(aspect.ratio = 1)

by_type <- expand[,.(N=sum(N)),by=list(PLATFORM,FILE,TYPE)]
by_type_plat <- by_type[,.(N=mean(N)),by=list(PLATFORM,TYPE)]
#mean(by_type[!(FILE %in% c("L.V.268","L.V.269","L.V.NANO")) & TYPE=="Transition",N])
#sd(by_type[!(FILE %in% c("L.V.268","L.V.269","L.V.NANO")) & TYPE=="Transition",N])

# Relative frequency of transitions

p3 <- ggplot(by_type[TYPE=="Transition"],aes(FILE,N,fill=PLATFORM))+
  geom_bar(stat="identity",position=position_dodge(),alpha=0.6,color="black",width=0.5)+
  #geom_point(alpha=0.6,size=3)+
  labs(y=NULL,x=NULL)+
  guides(fill=FALSE)+
  coord_flip()+
  scale_y_continuous(limits=c(0,100),breaks=c(0,50,100),position="right")+
  theme(aspect.ratio = 1,
        axis.text.y = element_blank(),
        axis.ticks.y=element_blank())

# Relative frequency of transversions

p4 <- ggplot(by_type[TYPE=="Transversion"],aes(FILE,N,fill=PLATFORM))+
  geom_bar(stat="identity",position=position_dodge(),alpha=0.6,color="black",width=0.5)+
  #geom_point(alpha=0.6,size=3)+
  labs(y=NULL,x=NULL)+
  guides(fill=FALSE)+
  coord_flip()+
  scale_y_continuous(limits=c(0,100),breaks=c(0,50,100),position="right")+
  theme(aspect.ratio = 1,
        axis.text.y = element_blank(),
        axis.ticks.y=element_blank())
grid.arrange(p1,p2,p3,p4,nrow=1)

# Calculate relationship between GC content and normalized coverage

gc <- unique(dt[,.(CHROM,NAME,POS,ID,REF.x,GC_CONTEXT)])
gc[,index:=seq(1,.N),by=NAME]

i <- gc[,.(index=c(1,cpts(cpt.mean(GC_CONTEXT, test.stat="Normal", method = "BinSeg", penalty = "Manual", pen.value = "log(2*log(n))",Q=4,
                                   minseglen=100)))),by=list(NAME)]
i[,BREAK:=1]
gc <- merge(gc,i,by=c("NAME",'index'),all.x = TRUE)
gc[is.na(BREAK),BREAK:=0]
gc[,BREAK:=replace(BREAK,seq_len(.N)==.N,0),by=NAME]
gc[,LABEL:=paste("gc",cumsum(BREAK),sep=""),by=NAME]

size <- gc[,.N,by=list(NAME,LABEL)]
content <- gc[,.(CONTENT=.N),by=list(NAME,LABEL,REF.x)][REF.x %in% c('C','G')]
content <- content[,.(CONTENT=sum(CONTENT)),by=list(NAME,LABEL)]
content <- merge(size,content,by=c('NAME','LABEL'))
content[,GC:=CONTENT/N*100]

gc <- merge(gc,content[,.(NAME,LABEL,GC,N)],by=c('NAME','LABEL'))
fragments <- unique(gc[,.(NAME,LABEL,GC,N)])
fragments <- fragments[,.(LABEL,FOLD=c(1,GC[-1*length(GC)]/c(GC[-1]))),by=NAME]

gc <- merge(gc,fragments,by=c("NAME","LABEL"))

gc <- merge(gc,dt[,.(FILE,PLATFORM,NAME,ID,READS_ALL,MISMATCHES,NORM_FACTOR)],by=c("NAME","ID"))

limits <- c(min(gc$GC)-1,25,50,75,max(gc$GC)+1)
gc[,GC_LABEL:=cut(GC,breaks=limits,labels=c(1,2,3,4))]

ggplot(gc[,.SD[sample(.N, min(1000000,.N))]],aes(GC,READS_ALL/NORM_FACTOR,color=PLATFORM,by=FILE))+
  geom_smooth()+
  labs(y="Normalized coverage",x="GC content (%)")+
  #scale_y_continuous(limits=c(0.5,2.5),breaks=c(0.5,1.5,2.5))+
  scale_y_continuous(limits=c(0.5,2.5),breaks=seq(0.5,2.5,0.2))+
  guides(color=FALSE)+
  geom_vline(xintercept = c(30,65),linetype="dashed")+
  theme(aspect.ratio = 1)


ggplot(gc[GC<=30,.SD[sample(.N, min(1000000,.N))]],aes(GC,READS_ALL/NORM_FACTOR,color=PLATFORM,by=FILE))+
  geom_smooth(method="lm")+
  labs(y="Normalized coverage",x="GC content (%)")+
  #scale_y_continuous(limits=c(0.5,2.5),breaks=c(0.5,1.5,2.5))+
  scale_y_continuous(limits=c(0.5,2.5),breaks=seq(0.5,2.5,0.2))+
  guides(color=FALSE)+
  geom_vline(xintercept = c(30,65),linetype="dashed")+
  theme(aspect.ratio = 1)

test <- gc[GC<=30,.SD[sample(.N, min(1000000,.N))]]
test[,.(coef(lm(READS_ALL/NORM_FACTOR~GC, .SD))),by=PLATFORM]


lm_eqn <- function(df){
  m <- lm(READS_ALL/NORM_FACTOR ~ GC, df);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));
}



gc[,GC:=GC/100]


gc_normal <- gc[GC>0.30 & GC<0.65,.SD[sample(.N, min(1000000,.N))]]
#gc_normal <- gc[GC>0.30 & GC<0.65,]
eq1 <- gc_normal[,.(a=unname(coef(lm(READS_ALL/NORM_FACTOR ~ GC))[1]),
                  b=unname(coef(lm(READS_ALL/NORM_FACTOR ~ GC))[2])),by=PLATFORM]
eq1[,SIGN:=ifelse(as.numeric(b)<0,-1,1),by=PLATFORM]
eq1[,b:=as.character(abs(as.numeric(b)))]
eq1 <- eq1[order(PLATFORM)]
eq1[,SEGMENT:=2]

gc_poor <- gc[GC<=0.30,.SD[sample(.N, min(1000000,.N))]]
#gc_poor <- gc[GC<=0.30,]
eq2 <- gc_poor[,.(a=unname(coef(lm(READS_ALL/NORM_FACTOR ~ GC))[1]),
                  b=unname(coef(lm(READS_ALL/NORM_FACTOR ~ GC))[2])),by=PLATFORM]
eq2[,SIGN:=ifelse(as.numeric(b)<0,-1,1),by=PLATFORM]
eq2[,b:=as.character(abs(as.numeric(b)))]
eq2 <- eq2[order(PLATFORM)]
eq2[,SEGMENT:=1]

gc_rich <- gc[GC>=0.65,.SD[sample(.N, min(1000000,.N))]]
#gc_rich <- gc[GC>=0.65,]
eq3 <- gc_rich[,.(a=unname(coef(lm(READS_ALL/NORM_FACTOR ~ GC))[1]),
                  b=unname(coef(lm(READS_ALL/NORM_FACTOR ~ GC))[2])),by=PLATFORM]
eq3[,SIGN:=ifelse(as.numeric(b)<0,-1,1),by=PLATFORM]
eq3[,b:=as.character(abs(as.numeric(b)))]
eq3 <- eq3[order(PLATFORM)]
eq3[,SEGMENT:=3]

eq <- rbind(eq1,eq2,eq3)
eq[SEGMENT==2,XLABEL:=0.50]
eq[SEGMENT==2,YLABEL:=seq(1.9,2.3,0.1)]
eq[SEGMENT==2,X:=0.3]
eq[SEGMENT==2,XEND:=0.65]

eq[SEGMENT==1,XLABEL:=0.15]
eq[SEGMENT==1,YLABEL:=seq(1.9,2.3,0.1)]
eq[SEGMENT==1,X:=0]
eq[SEGMENT==1,XEND:=0.30]

eq[SEGMENT==3,XLABEL:=0.80]
eq[SEGMENT==3,YLABEL:=seq(1.9,2.3,0.1)]
eq[SEGMENT==3,X:=0.65]
eq[SEGMENT==3,XEND:=1]

eq[,A:=as.numeric(as.character(a))]
eq[,B:=as.numeric(as.character(b))]

eq[,Y:=(A + (B*SIGN)*X)]
eq[,YEND:=(A + (B*SIGN)*XEND)]

eq[,a:=signif(as.numeric(as.character(a)),digits=2)]
eq[,b:=signif(as.numeric(as.character(b)),digits=2)]

eq[,EQ:=ifelse(as.numeric(SIGN)<0,
                as.character(as.expression(substitute(italic(y) == a - b %.% italic(x)))),
                as.character(as.expression(substitute(italic(y) == a + b %.% italic(x))))),by=list(SEGMENT,PLATFORM)]


ggplot()+
  geom_segment(data = eq,aes(x=X,xend=XEND,y=Y,yend=YEND,color=PLATFORM))+
  geom_text(data=eq,aes(x = XLABEL, y = YLABEL, label = EQ,color=PLATFORM),parse=TRUE)+
  labs(y="Normalized coverage",x="GC content")+
  scale_color_manual(values=c("HiSeqXTen/PCR-free"="#E185B7","HiSeq2500/PCR-based"="#EB8D82","NextSeq/PCR-based"="#00AEEB","BGI"="#ADAB33","ONT"="#00AB8A"))+
  scale_x_continuous(limits=c(0,1),breaks=c(0,0.5,1))+
  scale_y_continuous(limits=c(0.5,2.5),breaks=c(0.5,1.5,2.5))+
  #geom_abline(slope=2.5,intercept = 0.7,linetype="dashed")+
  guides(color=FALSE)+
  geom_vline(xintercept = c(0.30,0.65),linetype="dashed")+
  theme(aspect.ratio = 1)

eq[,RANGE:=((Y-YEND)/YEND)*100]
eq[SEGMENT==3,RANGE:=((YEND-Y)/Y)*100]
eq$SEGMENT <- as.factor(eq$SEGMENT)
levels(eq$SEGMENT) <- c("GC-poor","GC-normal","GC-rich")

ggplot(eq,aes(PLATFORM,RANGE,fill=PLATFORM))+
  geom_bar(stat="identity")+
  scale_fill_manual(values=c("HiSeqXTen/PCR-free"="#E185B7","HiSeq2500/PCR-based"="#EB8D82","NextSeq/PCR-based"="#00AEEB","BGI"="#FFCC04","ONT"="#00AB8A"))+
  facet_wrap(~SEGMENT)+
  labs(x=NULL,y="Normalized coverage range (%)")+
  guides(fill=FALSE)+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),aspect.ratio = 3)


test <- gc[,.(MEAN=mean(READS_ALL/NORM_FACTOR),SD=sd(READS_ALL/NORM_FACTOR),N=.N),by=list(PLATFORM,GC)]
ggplot(test,aes(GC,MEAN,color=PLATFORM))+
  geom_line()

ggplot(gc[GC<=30],aes(PLATFORM,READS_ALL/NORM_FACTOR))+
  geom_boxplot()+
  scale_y_continuous(limits=c(0.5,2.5),breaks=c(0.5,1.5,2.5))+
  theme(aspect.ratio = 1)

pvalue <- pairwise.wilcox.test(gc[GC<=30,READS_ALL/NORM_FACTOR],gc[GC<=30,PLATFORM],p.adjust.method="bonferroni")$p.value
summary(aov(READS_ALL/NORM_FACTOR ~ PLATFORM,data=gc[GC<=30]))
set.seed(245)
subset1 <- gc[GC<=30,.SD[sample(.N, min(100000,.N))]]
subset1[,TEST:="GC_poor"]
subset2 <- gc[,.SD[sample(.N, min(100000,.N))]]
subset2[,TEST:="Global"]
subset4 <- gc[GC>=55,.SD[sample(.N, min(100000,.N))]]
subset4[,TEST:="GC_rich"]
subset <- rbind(subset1,subset2,subset4)
#pairwise.wilcox.test(subset[,READS_ALL/NORM_FACTOR], subset[,PLATFORM],p.adjust.method = "bonferroni")
#pairwise.t.test(subset[,READS_ALL/NORM_FACTOR], subset[,PLATFORM],p.adjust.method = "bonferroni")

ggplot(subset,aes(PLATFORM,READS_ALL/NORM_FACTOR,color=TEST))+
  #geom_violin(trim=FALSE, alpha = 0.6,draw_quantiles = 0.5,adjust=2,width=0.5)+
  geom_boxplot()+
  #facet_wrap(~PLATFORM)+
  scale_y_continuous(limits=c(0.5,2.5),breaks=c(0.5,1.5,2.5))+
  theme(aspect.ratio = 1)+
  theme(aspect.ratio = 1)

stats <- subset[,.(mean(READS_ALL/NORM_FACTOR),sd(READS_ALL/NORM_FACTOR)),by=list(PLATFORM,TEST)]
stats <- dcast(data = stats[,.(PLATFORM,TEST,V1)],PLATFORM~TEST)
stats[,GC_poor:=(1-(GC_poor/Global))*100]
stats[,GC_rich:=(1-(GC_rich/Global))*100]
stats

subset <- gc[GC>=65,]
set.seed(123)
subset <- subset[,.SD[sample(.N, min(5000,.N))],by = PLATFORM]
pairwise.wilcox.test(subset[,READS_ALL/NORM_FACTOR], subset[,PLATFORM],p.adjust.method = "bonferroni")
pairwise.t.test(subset[,READS_ALL/NORM_FACTOR], subset[,PLATFORM],p.adjust.method = "bonferroni")

ggplot(subset,aes(PLATFORM,READS_ALL/NORM_FACTOR))+
  geom_boxplot()+
  scale_y_continuous(limits=c(0.5,2.5),breaks=c(0.5,1.5,2.5))+
  theme(aspect.ratio = 1)
# Normalized coverage around simple repeats
sample(.N,3)
w <- 400

# annotation of simple repeats and mobile elements within sequin regions
annotation <- fread('sequin_features_hg38_per_base_2.8.tab')
sr <- annotation[!is.na(Repeat_ID),.(CHROM=Chromosome,POS=Position,NAME=Name,R.ID=Repeat_ID,R.SIZE=Repeat_Overall_Size,R.NAME=NA,R.SEQ=Repeat,R.TYPE="SR")]
me <- annotation[!is.na(Mobile_Element_ID),.(CHROM=Chromosome,POS=Position,NAME=Name,R.ID=Mobile_Element_ID,R.SIZE=Mobile_Element_Size,R.NAME=Mobile_Element,R.SEQ=NA,R.TYPE="ME")]
repeats <- rbind(sr,me)
# get the start postion of each repeat
rpos <- repeats[,.(R.POS=min(POS)),by=R.ID]
repeats <- merge(repeats,rpos,by="R.ID")
repeats <- unique(repeats[,.(R.ID,CHROM,START=R.POS-w,END=R.POS,R.TYPE)])
repeats <- repeats[,.(CHROM,POS=seq(START,END),R.POS=END,R.TYPE),by=R.ID]
repeats <- merge(dt[,.(CHROM,POS,NAME,REF.x,FILE,PLATFORM,READS_ALL,NORM_FACTOR)],repeats,by=c("CHROM","POS"))
repeats <- repeats[order(FILE,CHROM,POS)]
repeats[,DISTANCE:=round(POS-R.POS)]


before <- repeats[R.TYPE=="SR"]
before[,PANEL:="Before"]

# ggplot(sr,aes(DISTANCE,READS_ALL/NORM_FACTOR,by=FILE,color=PLATFORM))+
#   geom_smooth(alpha=0.5)+
#   #scale_x_continuous(limits=c(-30,30),breaks=seq(-30,30,10))+
#   labs(x="Distance to simple repeat",y="Normalized coverage")+
#   guides(color=FALSE)+
#   theme(aspect.ratio = 1)

###### After Repeat

annotation <- fread('sequin_features_hg38_per_base_2.8.tab')
#annotation <- annotation[!is.na(Repeat_ID) | !is.na(Mobile_Element_ID),.(CHROM=Chromosome,POS=Position,NAME=Name,Repeat_ID,Repeat,Repeat_Size,Repeat_Overall_Size,Mobile_Element,Mobile_Element_ID,Mobile_Element_Size)]
sr <- annotation[!is.na(Repeat_ID),.(CHROM=Chromosome,POS=Position,NAME=Name,R.ID=Repeat_ID,R.SIZE=Repeat_Overall_Size,R.NAME=NA,R.SEQ=Repeat,R.TYPE="SR")]
me <- annotation[!is.na(Mobile_Element_ID),.(CHROM=Chromosome,POS=Position,NAME=Name,R.ID=Mobile_Element_ID,R.SIZE=Mobile_Element_Size,R.NAME=Mobile_Element,R.SEQ=NA,R.TYPE="ME")]
repeats <- rbind(sr,me)
rpos <- repeats[,.(R.POS=max(POS)),by=R.ID]
repeats <- merge(repeats,rpos,by="R.ID")
repeats <- unique(repeats[,.(R.ID,CHROM,START=R.POS,END=R.POS+w,R.TYPE)])
repeats <- repeats[,.(CHROM,POS=seq(START,END),R.POS=START,R.TYPE),by=R.ID]
repeats <- merge(dt[,.(CHROM,POS,NAME,REF.x,FILE,PLATFORM,READS_ALL,NORM_FACTOR)],repeats,by=c("CHROM","POS"))
repeats <- repeats[order(FILE,CHROM,POS)]
repeats[,DISTANCE:=round(POS-R.POS)]

after <- repeats[R.TYPE=="SR"]
after[,PANEL:="After"]

sr <- rbind(before,after)


ggplot(sr,aes(DISTANCE,READS_ALL/NORM_FACTOR,by=FILE,color=PLATFORM))+
  geom_smooth(alpha=0.5)+
  #scale_y_continuous(limits=c(1.2,1.5),breaks=c(1.2,1.5))+
  labs(x="Distance to simple repeat",y="Normalized coverage")+
  guides(color=FALSE)+
  theme(aspect.ratio = 1)


#####

all.files <- list.files(pattern = ".all_metrics.tab")
l <- lapply(all.files, fread, header=TRUE)
all.files <- str_remove(all.files,".sequin.sorted.full.all_metrics.tab")
l <- lapply(seq_along(all.files), function(x,y,i) x[[i]][,FILE:=y[i]],x=l,y=all.files)
dt <- rbindlist(l)
dt[,PLATFORM:=as.factor(FILE)]
levels(dt$PLATFORM) <- c("HiSeqXTen/PCR-free","HiSeqXTen/PCR-free","BGI","ONT","HiSeq2500/PCR-based","HiSeq2500/PCR-based","NextSeq/PCR-based","NextSeq/PCR-based")

gc <- data.table(FILE=character(),CHROM=character(),POS=numeric(),REF=character(),ALT=character(),AD=numeric(),FREQ=numeric(),PVAL=numeric(),GQ=numeric(),LABEL=character(),PLATFORM=character())


for (i in c('NA12878_SEQUIN_HISEQ_X_TEN_PCR_FREE_2','SEQUIN_HISEQ_2500_PCR_BASED_1','SEQUIN_HISEQ_2500_PCR_BASED_2','SEQUIN_NEXTSEQ_500_PCR_BASED_1','SEQUIN_NEXTSEQ_500_PCR_BASED_2','NA12878_SEQUIN_MGISEQ2000','NA12878_SEQUIN_PROMETHION')) {
  tp <- fread(sprintf('zcat %s.sequin.sorted.raw_varscan_variants.no_filter_gc/tp.vcf.gz',i),skip="#CHROM")
  tp[, AD := tstrsplit(Sample1, ":", fixed=TRUE)[6]]
  tp[,AD := as.numeric(AD)]
  tp[, FREQ := tstrsplit(Sample1, ":", fixed=TRUE)[7]]
  tp[,FREQ := as.numeric(str_remove(FREQ,"%"))/100]
  tp[, PVAL := tstrsplit(Sample1, ":", fixed=TRUE)[8]]
  tp[,PVAL:=as.numeric(PVAL)]
  tp[, GQ := tstrsplit(Sample1, ":", fixed=TRUE)[2]]
  tp[,GQ:=as.numeric(GQ)]
  tp[,LABEL:="TP"]
  tp <- tp[,c(1,2,4,5,11,12,13,14,15)]
  tp[,FILE:=sprintf("%s",i)]
  
  fp <- fread(sprintf('zcat %s.sequin.sorted.raw_varscan_variants.no_filter_gc/fp.vcf.gz',i),skip="#CHROM")
  fp[, AD := tstrsplit(Sample1, ":", fixed=TRUE)[6]]
  fp[,AD := as.numeric(AD)]
  fp[, FREQ := tstrsplit(Sample1, ":", fixed=TRUE)[7]]
  fp[,FREQ := as.numeric(str_remove(FREQ,"%"))/100]
  fp[, PVAL := tstrsplit(Sample1, ":", fixed=TRUE)[8]]
  fp[,PVAL:=as.numeric(PVAL)]
  fp[, GQ := tstrsplit(Sample1, ":", fixed=TRUE)[2]]
  fp[,GQ:=as.numeric(GQ)]
  fp[,LABEL:="FP"]
  fp <- fp[,c(1,2,4,5,11,12,13,14,15)]
  fp[,FILE:=sprintf("%s",i)]
  
  all <- rbind(tp,fp)
  all <- merge(all,unique(dt[,.(FILE,PLATFORM)]),by="FILE")
  names(all) <- c("FILE","CHROM","POS","REF","ALT","AD","FREQ","PVAL","GQ","LABEL","PLATFORM")
  gc <- rbind(gc,all)
}

gc <- gc[order(FILE,-FREQ)]
gc[LABEL=="FP",COUNT:=1]
gc[LABEL=="TP",COUNT:=0]
gc[,CUMSUM:=cumsum(COUNT),by=list(FILE)]
gc <- gc[order(FILE,FREQ)]

# ggplot(gc,aes(FREQ,CUMSUM,by=FILE,color=PLATFORM))+
#   geom_line()+
#   theme(aspect.ratio = 1)+
#   labs(x="Allele frequency",y="False positive rate",title="Extreme GC : No filtering")

gc[,CUMSUM:=CUMSUM/max(gc$CUMSUM)]
auc <- gc[,.(MESS::auc(FREQ,CUMSUM)),by=PLATFORM]
auc[,LABEL:=paste("AUC=",signif(V1,2),sep="")]
auc[,X:=0.25]
auc[,Y:=seq(0.5,0.9,0.1)]

ggplot(gc)+
  geom_line(aes(FREQ,CUMSUM,by=FILE,color=PLATFORM))+
  theme(aspect.ratio = 1)+
  scale_x_continuous(limits=c(0,0.4),breaks=c(0,0.1,0.2,0.3))+
  scale_y_continuous(limits=c(0,1),breaks=c(0,0.5,1))+
  geom_text(data=auc,aes(x=X,y=Y,label=LABEL,color=PLATFORM))+
  scale_color_manual(values=c("HiSeqXTen/PCR-free"="#E185B7","HiSeq2500/PCR-based"="#EB8D82","NextSeq/PCR-based"="#00AEEB","BGI"="#FFCC04","ONT"="#00AB8A"))+
  labs(x="Allele frequency",y="False positive rate")+
  guides(color=FALSE)



dev.off()
