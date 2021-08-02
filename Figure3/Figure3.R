setwd('/simon/ANDRE/Figures/github/Figure3')

library(cowplot)
library(ggpubr)
library(data.table)
library(stringr)

theme_set(theme_cowplot(font_size=20))

dodge <- position_dodge(width = 0.8)

whichmedian <- function(x) which.min(abs(x - median(x)))

pdf("Plots_Figure3.pdf",useDingbats=FALSE)

# Upload pysamstats statistics

all.files <- list.files(pattern = ".all_metrics.tab")
l <- lapply(all.files, fread, header=TRUE)
all.files <- str_remove(all.files,".sequin.sorted.full.all_metrics.tab")
l <- lapply(seq_along(all.files), function(x,y,i) x[[i]][,FILE:=y[i]],x=l,y=all.files)
dt <- rbindlist(l)

# Calculate a normalization factor based on the total sum of read counts for each library

norm_factor <- dt[,sum(READS_ALL),by=FILE][,.(NORM_FACTOR=V1/10^6),by=FILE]

dt <- merge(dt,norm_factor,by="FILE")

# Upload all positions containing variants within sequins

# both small variants
vcf1 <- fread('sequin_smallvariants_hg38_2.7.vcf',skip="#CHROM")
vcf1[,POS:=POS-1]
vcf1[,N:=nchar(REF)]
vcf1[N==1,END:=POS]
vcf1[N!=1,END:=POS+(N-1)]
vcf1 <- vcf1[,c(1,2,10,3)]
names(vcf1) <- c('CHROM','START','END','NAME')
vcf1 <- vcf1[,.(CHROM,POS=seq(START,END)),by=NAME]
vcf1[,ID:=paste(CHROM,POS,sep=':')]

# and structural variants
vcf2 <- fread('sequin_structural_hg38_2.6.vcf',skip="#CHROM")
vcf2 <- vcf2[ALT %in% c("<INS>","<DEL>")]
vcf2[,END:=str_extract(INFO,"END=([0-9-Z]+)")]
vcf2[,END:=as.numeric(str_remove(END,"END="))]
vcf2 <- vcf2[,c(1,2,9,3)]
names(vcf2) <- c('CHROM','START','END','NAME')
vcf2 <- vcf2[,.(CHROM,POS=seq(START,END)),by=NAME]
vcf2[,ID:=paste(CHROM,POS,sep=':')]

vcf <- rbind(vcf1,vcf2)

# Remove any position overlapping variants so that they are not wrongly accounted as errors
dt[,ID:=paste(CHROM,POS,sep=':')]
hp <- dt[!(ID %in% vcf$ID),]

# Upload annotation of homopolymer of varying sizes within sequin regions
annotation <- fread('sequin_regions_hg38_2.6.homopolymers_no_edges.tab',header=FALSE)
annotation <- annotation[!duplicated(annotation)]
annotation <- annotation[,.(V1,seq(V2-1,V3-2),V4,V5,V6),by=V7]
annotation[,V2:=V2+1]
names(annotation) <- c('HP_ID','CHROM','POS','NAME','HP_BASE','HP_LENGTH')

# Create a homopolymer length distribution (small <= 5bp, medium >5bp & <=15bp, large >15bp)
limits <- c(0,5,15,100)
hist <- unique(annotation[,.(HP_ID,HP_LENGTH)])
hist <- hist[,.N,by=HP_LENGTH]
hist[,HP_LABEL:=cut(HP_LENGTH,breaks=limits,labels=c("small","medium","large"))]

ggplot(hist,aes(HP_LENGTH,log(N)))+
  geom_line()+
  theme(aspect.ratio = 1)+
  geom_vline(xintercept = 5,linetype="dashed")+
  geom_vline(xintercept = 15,linetype="dashed")+
  labs(x="Homopolymer length (nt)", y="Frequency (Log2)")


# Total number of homopolymers of each length category
n <- hist[,.(sum(N)),by=HP_LABEL]
 
# restrict the pysamstats to positions overlapping homopolymers 
hp <- merge(hp,annotation,by=c("CHROM","POS","NAME"))

# Calculate frequency of different error types (mismatches, insertions and deletions) within homopolymers
hp <- hp[,.(FILE,PLATFORM,NAME,CHROM,POS,REF.x,READS_ALL,MISMATCHES,INSERTIONS,DELETIONS,HP_ID,HP_BASE,HP_LENGTH)]
hp <- hp[order(FILE,CHROM,POS)]
hp[,c("MISMATCHES_F","INSERTIONS_F","DELETIONS_F"):=list(MISMATCHES/READS_ALL,INSERTIONS/READS_ALL,DELETIONS/READS_ALL)]

limits <- c(min(hp$HP_LENGTH)-1,5,15,max(hp$HP_LENGTH)+1)
hp[,HP_LABEL:=cut(HP_LENGTH,breaks=limits,labels=c("small","medium","large"))]

# Calculate some stats for each platform
hp[,.(DEL=format(mean(DELETIONS/READS_ALL,na.rm = TRUE),scientific=FALSE),INS=format(mean(INSERTIONS/READS_ALL,na.rm = TRUE),scientific=FALSE)),by=list(PLATFORM)]

hp[,PLATFORM:=as.factor(PLATFORM)]

# Calculate the average mismatche, insertion and deletion frequency for each homopolymer length across the different libraries
hp_sum <- hp[,.(.N,M=sum(MISMATCHES_F,na.rm=TRUE),D=sum(DELETIONS_F,na.rm=TRUE),I=sum(INSERTIONS_F,na.rm=TRUE)),by=list(FILE,PLATFORM,HP_LENGTH)]
hp_sum[,HP_LABEL:=cut(HP_LENGTH,breaks=limits,labels=c("small","medium","large"))]

# Distribution average mismatch frequency
ggplot(hp_sum,aes(HP_LABEL,log(M/N),fill=PLATFORM,by=FILE))+
  geom_boxplot(alpha=0.6)+
  labs(x="Homopolymer Length",y="Mismatch fraction (Log)",color=NULL)+
  scale_y_continuous(limits=c(-7,-3),breaks=c(-7,-5,-3))+
  theme(aspect.ratio = 1)+
  guides(fill=FALSE)

# Distribution average deletion frequency
ggplot(hp_sum,aes(HP_LABEL,log(D/N),fill=PLATFORM,by=FILE))+
  geom_boxplot(alpha=0.6)+
  labs(x="Homopolymer Length",y="Deletion fraction (Log)",color=NULL)+
  scale_y_continuous(limits=c(-11,0),breaks=c(-10,-5,0))+
  theme(aspect.ratio = 1)+
  guides(fill=FALSE)

ggplot(hp_sum,aes(HP_LENGTH,log(D/N),color=PLATFORM,fill=PLATFORM,by=FILE))+
  #geom_boxplot(alpha=0.6)+
  geom_smooth()+
  labs(x="Homopolymer Length",y="Deletion fraction (Log)",color=NULL)+
  scale_y_continuous(limits=c(-11,0),breaks=c(-10,-5,0))+
  geom_vline(xintercept = c(5,15),linetype="dashed")+
  theme(aspect.ratio = 1)+
  guides(fill=FALSE,color=FALSE)

correlation_del <- hp_sum[HP_LENGTH<=15,.(cor(HP_LENGTH,D/N)),by=list(FILE)]
mean(correlation_del[,V1])
sd(correlation_del[,V1])

# Distribution average insertion frequency
ggplot(hp_sum,aes(HP_LABEL,log(I/N),fill=PLATFORM,by=FILE))+
  geom_boxplot(alpha=0.6)+
  labs(x="Homopolymer Length",y="Insertion fraction (Log)",color=NULL)+
  scale_y_continuous(limits=c(-13,-4),breaks=c(-12,-8,-4))+
  theme(aspect.ratio = 1)+
  guides(fill=FALSE)

hp_sum[,I:=I/N]
hp_sum[,D:=D/N]
hp_sum[,PLATFORM:=as.factor(PLATFORM)]
hp_sum[,FILE:=as.factor(FILE)]

kruskal.test(I ~ FILE, data = hp_sum)
kruskal.test(D ~ PLATFORM, data = hp_sum)

format(pairwise.wilcox.test(hp_sum[,D],hp_sum[,PLATFORM],p.adjust.method = "bonferroni")$p.value,scientific=FALSE)
adjust <- "BH"
a1 <- data.table(melt(format(pairwise.wilcox.test(hp_sum[HP_LABEL=="small" & FILE!= "NA12878_SEQUIN_PROMETHION",D],hp_sum[HP_LABEL=="small" & FILE!= "NA12878_SEQUIN_PROMETHION",PLATFORM],p.adjust.method = adjust)$p.value,scientific=FALSE)))
a1$Label <- "Deletion: small"
a2 <- data.table(melt(format(pairwise.wilcox.test(hp_sum[HP_LABEL=="medium" & FILE!= "NA12878_SEQUIN_PROMETHION",D],hp_sum[HP_LABEL=="medium" & FILE!= "NA12878_SEQUIN_PROMETHION",PLATFORM],p.adjust.method = adjust)$p.value,scientific = FALSE)))
a2$Label <- "Deletion: medium"
a3 <- data.table(melt(format(pairwise.wilcox.test(hp_sum[HP_LABEL=="large" & FILE!= "NA12878_SEQUIN_PROMETHION",D],hp_sum[HP_LABEL=="large" & FILE!= "NA12878_SEQUIN_PROMETHION",PLATFORM],p.adjust.method = adjust)$p.value,scientific = FALSE)))
a3$Label <- "Deletion: large"
a4 <- data.table(melt(format(pairwise.wilcox.test(hp_sum[HP_LABEL=="small" & FILE!= "NA12878_SEQUIN_PROMETHION",I],hp_sum[HP_LABEL=="small" & FILE!= "NA12878_SEQUIN_PROMETHION",PLATFORM],p.adjust.method = adjust)$p.value,scientific=FALSE)))
a4$Label <- "Insertion: small"
p.value <- data.table(rbind(a1,a2,a3,a4))
p.value[,value:=as.numeric(as.character(format(value,scientific=FALSE)))]
p.value <- p.value[!is.na(value)]
p.value[,Significance:=cut(value,breaks=c(-Inf,0.0001,0.001,0.01,0.05,1),labels=c("****","***","**","*","n.s."))]
levels(p.value$Var1) <- c("MGISEQ2000","NextSeq500","HiSeqXTen")
levels(p.value$Var2) <- c("HiSeq2500","MGISEQ2000","NextSeq500")
ggplot(p.value,aes(Var1,Var2,fill=value))+
  geom_tile()+
  facet_wrap(~Label,nrow=2,ncol=2)+
  geom_text(aes(label=Significance),color="white")+
  theme(aspect.ratio = 1,axis.text.x=element_text(angle=30,hjust = 1))+
  labs(x=NULL,y=NULL,fill="P-value")+
  scale_fill_gradient(limits=c(0,1),breaks=c(0,0.5,1))

axis.text.x=theme_text(angle=-90)

ggplot(hp_sum,aes(HP_LENGTH,log(I),color=PLATFORM,fill=PLATFORM,by=FILE))+
  #geom_boxplot(alpha=0.6)+
  geom_smooth()+
  labs(x="Homopolymer Length",y="Insertion fraction (Log)",color=NULL)+
  scale_y_continuous(limits=c(-13,-4),breaks=c(-12,-8,-4))+
  geom_vline(xintercept = c(5,15),linetype="dashed")+
  theme(aspect.ratio = 1)+
  guides(fill=FALSE,color=FALSE)

correlation_ins <- hp_sum[HP_LENGTH<=15,.(cor(HP_LENGTH,I/N)),by=list(FILE)]
mean(correlation_ins[,V1])
sd(correlation_ins[,V1])

hp_sum_stats <- hp_sum[,.(Avg.Deletion=mean(D/N),SD.Deletion=sd(D/N),Avg.Insertion=mean(I/N),SD.Insertion=sd(I/N)),by=list(FILE,HP_LABEL)]

hp_sum_stats <- merge(hp_sum_stats,hp_sum_stats[FILE=="NA12878_SEQUIN_PROMETHION",.(HP_LABEL,NANO.DEL=Avg.Deletion,NANO.INS=Avg.Insertion)],by="HP_LABEL")
hp_sum_stats[,NANO.FOLD.DEL:=log(NANO.DEL/Avg.Deletion)]
hp_sum_stats[,NANO.FOLD.INS:=log(NANO.INS/Avg.Insertion)]

hp_sum_stats[FILE!="NA12878_SEQUIN_PROMETHION",.(mean(NANO.FOLD.DEL),sd(NANO.FOLD.DEL)),by=list(HP_LABEL)]

hp_sum_stats[FILE!="NA12878_SEQUIN_PROMETHION",.(mean(NANO.FOLD.INS),sd(NANO.FOLD.INS)),by=list(HP_LABEL)]

hp_sum_stats <- hp_sum[,.(Avg.Deletion=mean(D/N),SD.Deletion=sd(D/N),Avg.Insertion=mean(I/N),SD.Insertion=sd(I/N)),by=list(FILE,HP_LABEL)]

pcr.free <- hp_sum_stats[FILE %in% c("NA12878_SEQUIN_HISEQ_X_TEN_PCR_FREE_1","NA12878_SEQUIN_HISEQ_X_TEN_PCR_FREE_1"),.(HP_LABEL,NANO.DEL=Avg.Deletion,NANO.INS=Avg.Insertion)][,.(FREE.DEL=mean(NANO.DEL),FREE.INS=mean(NANO.INS)),by=HP_LABEL]
hp_sum_stats <- merge(hp_sum_stats,pcr.free,by="HP_LABEL")
hp_sum_stats[,FREE.FOLD.DEL:=FREE.DEL/Avg.Deletion]
hp_sum_stats[,FREE.FOLD.INS:=FREE.INS/Avg.Insertion]

hp_sum_stats[!(FILE %in% c("NA12878_SEQUIN_HISEQ_X_TEN_PCR_FREE_1","NA12878_SEQUIN_HISEQ_X_TEN_PCR_FREE_2")),.(mean(FREE.FOLD.DEL),sd(FREE.FOLD.DEL)),by=list(HP_LABEL)]

hp_sum_stats[!(FILE %in% c("NA12878_SEQUIN_HISEQ_X_TEN_PCR_FREE_1","NA12878_SEQUIN_HISEQ_X_TEN_PCR_FREE_2")),.(mean(FREE.FOLD.INS),sd(FREE.FOLD.INS)),by=list(HP_LABEL)]

#########################################################

b_names <- c('BAT-25','BAT-26','D2S123','D5S346','D17S250')

all.files <- list.files(pattern = ".sequin.bethesda.sorted.tab")
l <- lapply(all.files, fread, header=FALSE)
all.files <- str_remove(all.files,".sequin.bethesda.sorted.tab")
l <- lapply(seq_along(all.files), function(x,y,i) x[[i]][,FILE:=y[i]],x=l,y=all.files)
bethesda <- rbindlist(l)
names(bethesda) <- c('NAME','REPEAT','SIZE','READ','SEQ','LENGTH','FILE')
bethesda[NAME=='D2S123',SIZE:=42]
bethesda[NAME=='D5S346',SIZE:=40]

#bethesda <- bethesda[!(NAME %in% b_names)]
bethesda[,DISTANCE:=SIZE-LENGTH]

bethesda <- merge(bethesda,unique(dt[,.(FILE,PLATFORM)]),by="FILE")

bethesda_matrix <- merge(bethesda[,.(TOTAL=.N),by=list(FILE,NAME)],bethesda[DISTANCE==0,.(MATCH=.N),by=list(FILE,PLATFORM,NAME)],by=c("FILE","NAME"))
bethesda_matrix[,MATCH:=MATCH/TOTAL]
test <- data.table(dcast(bethesda_matrix[,.(FILE,NAME,MATCH)],FILE~NAME,value.var = "MATCH"))
test[is.na(test)] <- 0
bethesda_matrix <- melt(test,id.vars = "FILE",variable.name = "NAME",value.name = "MATCH")
bethesda_matrix[,TYPE:="Stable"]
bethesda_matrix[NAME %in% b_names,TYPE:="Unstable"]
mid <- 0.5
ggplot(bethesda_matrix,aes(FILE,NAME,fill=MATCH))+
  geom_tile()+
  scale_fill_gradient2(midpoint = mid, low = "blue", mid = "white",
                        high = "red", space = "Lab",limits=c(0,1),breaks=c(0,0.5,1))+
  facet_grid(TYPE~.,scales="free_y",space = "free_y")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  labs(x=NULL,y=NULL,fill=NULL)

matrix_values <- merge(bethesda_matrix,unique(dt[,.(FILE,PLATFORM)]),by="FILE")
matrix_values[,.(mean(MATCH),sd(MATCH)),by=list(PLATFORM,TYPE)]

unstable_comparison <- data.table(dcast(bethesda_matrix,FILE+NAME~TYPE,value.var = "MATCH"))
unstable_comparison <- unstable_comparison[,.(Stable=mean(Stable,na.rm=TRUE),Unstable=mean(Unstable,na.rm=TRUE)),by=FILE]
unstable_comparison[,Ratio:=Unstable/Stable]
mean(unstable_comparison[FILE!="NA12878_SEQUIN_PROMETHION",1-Ratio])
sd(unstable_comparison[FILE!="NA12878_SEQUIN_PROMETHION",1-Ratio])
##############################################

all.files <- list.files(pattern = ".sequin.homopolymers.sorted.tab")
l <- lapply(all.files, fread, header=FALSE)
all.files <- str_remove(all.files,".sequin.homopolymers.sorted.tab")
l <- lapply(seq_along(all.files), function(x,y,i) x[[i]][,FILE:=y[i]],x=l,y=all.files)
homo_reads <- rbindlist(l)
names(homo_reads) <- c('NAME','REPEAT','SIZE','READ','SEQ','LENGTH','FILE')
homo_reads[,SIZE:=NULL]

homo_sizes <- fread('sequin_regions_hg38_2.6.homopolymers_no_edges.tab',header=FALSE)
homo_sizes <- homo_sizes[,.(NAME=V7,SIZE=V6)]

homo_reads <- merge(homo_reads,homo_sizes,by="NAME")

na <- homo_reads[is.na(SEQ)]

homo_reads <- merge(homo_reads,unique(dt[,.(FILE,PLATFORM)]),by="FILE")
homo_reads[,DISTANCE:=SIZE-LENGTH]

homo_reads <- homo_reads[!is.na(SEQ)]
homo_reads <- homo_reads[SIZE > 10]
homo_reads_sum <- homo_reads[!is.na(SEQ),.N,by=list(FILE,PLATFORM,DISTANCE)]
homo_reads_total <- homo_reads_sum[,.(TOTAL=sum(N)),by=list(FILE)]
homo_reads_sum <- merge(homo_reads_sum,homo_reads_total,by='FILE')

ggplot(homo_reads_sum,aes(DISTANCE,N/TOTAL,fill=PLATFORM,by=FILE))+
  geom_bar(stat="identity",position=position_dodge(),alpha=0.6,color="black")+
  scale_x_continuous(limits=c(-3,3),breaks = seq(-3,3))+
  #scale_y_continuous(limits=c(0,0.9),breaks=c(0,0.4,0.8))+
  guides(fill=FALSE)+
  theme(aspect.ratio = 1)+
  labs(x="Distance to true length",y="Relative frequency")

exact <- homo_reads_sum[DISTANCE==0,.(FRACTION=N/TOTAL,PLATFORM)][,(FRACTION=mean(FRACTION)),by=PLATFORM]
exact[!(PLATFORM %in% c("xten/pcr_free","nanopore")),.(mean(V1),sd(V1))]

# ##############################################
# 
# av <- dt[!(ID %in% vcf$ID),]
# av <- av[,.(FILE,PLATFORM,CHROM,POS,NAME,MISMATCHES_F=MISMATCHES/READS_ALL,INSERTIONS_F=INSERTIONS/READS_ALL,DELETIONS_F=DELETIONS/READS_ALL)]
# av <- av[,.(MISMATCHES=sum(MISMATCHES_F,na.rm=TRUE)/.N,DELETIONS=sum(DELETIONS_F,na.rm=TRUE)/.N,INSERTIONS=sum(INSERTIONS_F,na.rm=TRUE)/.N),by=list(FILE,PLATFORM)]
# av <- melt(av,id.vars = c("FILE","PLATFORM"),variable.name = "ERROR_TYPE",value.name = "VALUE")
# av[,FILE:=as.factor(FILE)]
# levels(av$FILE) <- 1:length(levels(av$FILE))
# av$FILE <- as.numeric(av$FILE)
# 
# bethesda <- fread('bethesda_microsatelites.bed',header=FALSE)
# bethesda[, c("V4","V5") := tstrsplit(V4, ":", fixed=TRUE)]
# names(bethesda) <- c('CHROM','START','END','NAME.B','REPEAT.B')
# bethesda <- bethesda[,.(CHROM,POS=seq(START,END-1),REPEAT.B),by=NAME.B]
# 
# bethesda <- merge(bethesda,dt,by=c('CHROM','POS'))
# bethesda[,ID:=paste(CHROM,POS,sep=':')]
# 
# bethesda <- bethesda[!(ID %in% vcf$ID),]
# 
# bethesda <- bethesda[,.(FILE,PLATFORM,CHROM,POS,NAME.B,MISMATCHES_F=MISMATCHES/READS_ALL,INSERTIONS_F=INSERTIONS/READS_ALL,DELETIONS_F=DELETIONS/READS_ALL)]
# bethesda  <- bethesda[,.(MISMATCHES=sum(MISMATCHES_F,na.rm=TRUE)/.N,DELETIONS=sum(DELETIONS_F,na.rm=TRUE)/.N,INSERTIONS=sum(INSERTIONS_F,na.rm=TRUE)/.N),by=list(FILE,PLATFORM,NAME.B)]
# bethesda  <- melt(bethesda ,id.vars = c("FILE","PLATFORM","NAME.B"),variable.name = "ERROR_TYPE",value.name = "VALUE")
# 
# bethesda[,FILE:=as.factor(FILE)]
# levels(bethesda$FILE) <- 1:length(levels(bethesda$FILE))
# bethesda$FILE <- as.numeric(bethesda$FILE)

# ggplot(bethesda[ERROR_TYPE=="DELETIONS"],aes(FILE,log(VALUE),color=PLATFORM,by=as.factor(FILE)))+
#   geom_point(pch=21,size=3)+
#   geom_segment(data=av[ERROR_TYPE=="DELETIONS"],aes(x=FILE-0.25,xend=FILE+0.25,y=log(VALUE),yend=log(VALUE)),color="red",size=1)+
#   theme(aspect.ratio = 1,axis.ticks.x=element_blank(),axis.text.x=element_blank())+
#   guides(color=FALSE)+
#   labs(x=NULL,y="Average frequency (Log)",title="Deletions")
# 
# ggplot(bethesda[ERROR_TYPE=="INSERTIONS"],aes(FILE,log(VALUE),color=PLATFORM,by=as.factor(FILE)))+
#   geom_point(pch=21,size=3)+
#   geom_segment(data=av[ERROR_TYPE=="INSERTIONS"],aes(x=FILE-0.25,xend=FILE+0.25,y=log(VALUE),yend=log(VALUE)),color="red",size=1)+
#   theme(aspect.ratio = 1,axis.ticks.x=element_blank(),axis.text.x=element_blank())+
#   guides(color=FALSE)+
#   labs(x=NULL,y="Average frequency (Log)",title="Insertions")


dev.off()