#setwd('/simon/ANDRE/Figures/Figure5/dilution_series')
setwd('/simon/ANDRE/Figures/github/Figure6/HLA')
pdf('plots_figure5_v4.pdf',width=15,paper="a4r",useDingbats = FALSE)


all.files <- list.files(pattern = ".HLA.calibrated.ref.sorted.variation.txt")
l <- lapply(all.files, fread, header=TRUE)
all.files <- str_remove(all.files,".HLA.calibrated.ref.sorted.variation.txt")
l <- lapply(seq_along(all.files), function(x,y,i) x[[i]][,FILE:=y[i]],x=l,y=all.files)
dt <- rbindlist(l)

dt$PLATFORM <- as.factor(dt$FILE)
#levels(dt$PLATFORM) <-  c("HiSeqXTen/PCR-free","HiSeq2500/PCR-based","HiSeq2500/PCR-based","NextSeq/PCR-based","NextSeq/PCR-based","BGI","ONT")
levels(dt$PLATFORM) <- c("HiSeqXTen/PCR-free","BGI","ONT","HiSeq2500/PCR-based","HiSeq2500/PCR-based","NextSeq/PCR-based","NextSeq/PCR-based")

mean_cov <- dt[,.(MEAN=mean(READS_ALL)),by=list(FILE)]

ggplot(dt,aes(FILE,READS_ALL,color=PLATFORM))+
  geom_boxplot()+
  scale_color_manual(values=c("HiSeqXTen/PCR-free"="#E185B7","HiSeq2500/PCR-based"="#EB8D82","NextSeq/PCR-based"="#00AEEB","BGI"="#ADAB33","ONT"="#00AB8A"))+
  labs(x=NULL,y="Coverage")+
  guides(color=FALSE)+
  theme(aspect.ratio = 1)+
  geom_hline(yintercept = 100,linetype="dashed")

exons <- fread('all_exons.bed')
names(exons) <- c('CHROM','START','END','EXON')
exons <- unique(exons[EXON %in% c('EX2','EX3') & endsWith(CHROM,"_A")])
exons <- exons[,.(POS=seq(START,END)),by=list(CHROM,EXON)]

dt <- merge(dt,exons,by=c("CHROM","POS"))
dt[,MAX:=max(A,C,G,T),by=list(FILE,CHROM,POS)]
dt[MAX==A,CONSENSUS:="A"]
dt[MAX==C,CONSENSUS:="C"]
dt[MAX==`T`,CONSENSUS:="T"]
dt[MAX==G,CONSENSUS:="G"]

variants <- fread('hla_variants.tab')
names(variants) <- c('CHROM','POS','A1','A2','STATUS')

dt <- merge(dt,variants,by=c("CHROM","POS"))

dt$EXON <- as.factor(dt$EXON)
levels(dt$EXON) <- c("Exon 2","Exon 3")

ggplot(dt,aes(EXON,READS_ALL,by=FILE,color=PLATFORM))+
  geom_boxplot()+
  scale_color_manual(values=c("HiSeqXTen/PCR-free"="#E185B7","HiSeq2500/PCR-based"="#EB8D82","NextSeq/PCR-based"="#00AEEB","BGI"="#ADAB33","ONT"="#00AB8A"))+
  labs(x=NULL,y="Coverage")+
  guides(color=FALSE)+
  theme(aspect.ratio = 1)+
  geom_hline(yintercept = 100,linetype="dashed")
  
dt[,MISMATCHES:=MISMATCHES/READS_ALL]
dt[,DELETIONS:=DELETIONS/READS_ALL]
dt[,INSERTIONS:=INSERTIONS/READS_ALL]

error <- dt[STATUS!="Mimatch"]

error <- error[,.(MISMATCHES=sum(MISMATCHES)/length(MISMATCHES),
         DELETIONS=sum(DELETIONS)/length(DELETIONS),
         INSERTIONS=sum(INSERTIONS)/length(INSERTIONS)),by=list(FILE,PLATFORM,EXON)]

ggplot(error,aes(EXON,MISMATCHES,by=FILE,fill=PLATFORM))+
  geom_bar(stat="identity",position = position_dodge(),color="black")+
  scale_fill_manual(values=c("HiSeqXTen/PCR-free"="#E185B7","HiSeq2500/PCR-based"="#EB8D82","NextSeq/PCR-based"="#00AEEB","BGI"="#ADAB33","ONT"="#00AB8A"))+
  labs(x=NULL,y="Average mismatch")+
  guides(fill=FALSE)+
  theme(aspect.ratio = 1)

ggplot(error,aes(EXON,DELETIONS,by=FILE,fill=PLATFORM))+
  geom_bar(stat="identity",position = position_dodge(),color="black")+
  scale_fill_manual(values=c("HiSeqXTen/PCR-free"="#E185B7","HiSeq2500/PCR-based"="#EB8D82","NextSeq/PCR-based"="#00AEEB","BGI"="#ADAB33","ONT"="#00AB8A"))+
  labs(x=NULL,y="Average deletion")+
  guides(fill=FALSE)+
  theme(aspect.ratio = 1)

ggplot(error,aes(EXON,INSERTIONS,by=FILE,fill=PLATFORM))+
  geom_bar(stat="identity",position = position_dodge(),color="black")+
  scale_fill_manual(values=c("HiSeqXTen/PCR-free"="#E185B7","HiSeq2500/PCR-based"="#EB8D82","NextSeq/PCR-based"="#00AEEB","BGI"="#ADAB33","ONT"="#00AB8A"))+
  labs(x=NULL,y="Average insertion")+
  guides(fill=FALSE)+
  theme(aspect.ratio = 1)

consensus <- dt
start <- consensus[,.(START=min(POS),END=max(POS)),by=list(FILE,CHROM)]
consensus <- dt[STATUS=="Mimatch",CONSENSUS:="N"]
consensus_seq <- consensus[,.(CHROM,POS,REF,EXON,FILE,PLATFORM,CONSENSUS)]
consensus_seq <- consensus_seq[order(FILE,PLATFORM,CHROM,POS,EXON)]
consensus_seq <- consensus_seq[,.(REF=paste(REF,collapse = ""),CONSENSUS=paste(CONSENSUS,collapse = "")),by=list(FILE,PLATFORM,CHROM,EXON)]
consensus_seq <- consensus_seq[order(CHROM,EXON,FILE)] 
levels(consensus_seq$EXON) <- c("Exon2","Exon3")
write.table(consensus_seq,file="consensus_seq.tab",quote=FALSE,sep="\t",col.names=FALSE,row.names=FALSE)

consensus <- dt[STATUS!="Mimatch"]
consensus <- merge(consensus,start,by=c("FILE","CHROM"))
consensus[,RELATIVE_POS.START:=POS-START]
consensus[,RELATIVE_POS.END:=POS-END]
consensus[,EDIT:=0]
consensus[REF!=CONSENSUS,EDIT:=1]
consensus <- consensus[,.(EDIT=sum(EDIT)),by=list(FILE,PLATFORM,CHROM)]
consensus$CHROM <- as.factor(consensus$CHROM)
levels(consensus$CHROM) <- c("HLA-A","HLA-B","HLA-C","HLA-DQB1")

ggplot(consensus,aes(CHROM,EDIT,by=FILE,color=PLATFORM))+
  geom_point(alpha=0.6)+
  scale_color_manual(values=c("HiSeqXTen/PCR-free"="#E185B7","HiSeq2500/PCR-based"="#EB8D82","NextSeq/PCR-based"="#00AEEB","BGI"="#ADAB33","ONT"="#00AB8A"))+
  labs(x=NULL,y="Number of edits")+
  guides(color=FALSE)+
  theme(aspect.ratio = 1)

#################################
#################################


# Read reference file with the sequin's HLA 
ref <- fread('hla_sequins.tab',header=FALSE)
ref <- ref[,.(V1,V3,V5)]
names(ref) <- c('Sequin','Ref_Allele','AF')
ref[,c("V1","V2") := tstrsplit(Ref_Allele,"*",fixed=TRUE)]
ref[,V1:=str_remove(V1,"HLA-")]
ref <- dcast(ref[,.(N=seq(1,.N),Ref_Allele),by=V1],V1~N,value.var = "Ref_Allele")
names(ref) <- c("Locus","Ref_Allele1","Ref_Allele2")

dt <- list()
h <- 1
for (i in c('NA12878_SEQUIN_HISEQ_X_TEN_PCR_FREE_2','SEQUIN_HISEQ_2500_PCR_BASED_1','SEQUIN_HISEQ_2500_PCR_BASED_2','SEQUIN_NEXTSEQ_500_PCR_BASED_1','SEQUIN_NEXTSEQ_500_PCR_BASED_2','NA12878_SEQUIN_MGISEQ2000','NA12878_SEQUIN_PROMETHION')) {
  for (j in c(as.character(seq(100,5,-5)))) {
    # if (i=='NANO' & j=='5') {
    #   j <- '5'
    # }
    call <- fread(sprintf('dilution_series/%sHLA%sx/hla/R1_bestguess_G.txt',i,j))
    call[,FILE:=sprintf('%s',i)]
    call[,Expected_Coverage:=j]
    dt[[h]] <- call
    h <- h+1
  }
}
dt <- rbindlist(dt)
# dt[Expected_Coverage=="b",Expected_Coverage:="4"]
# dt[Expected_Coverage=="c",Expected_Coverage:="3"]
# dt[Expected_Coverage=="d",Expected_Coverage:="2"]
# dt[Expected_Coverage=="e",Expected_Coverage:="1"]
# dt[Expected_Coverage=="f",Expected_Coverage:="0.5"]
# dt[Expected_Coverage=="g",Expected_Coverage:="0.1"]
dt[,Expected_Coverage:=as.numeric(Expected_Coverage)]

dt$PLATFORM <- as.factor(dt$FILE)
#levels(dt$PLATFORM) <-  c("HiSeqXTen/PCR-free","HiSeq2500/PCR-based","HiSeq2500/PCR-based","NextSeq/PCR-based","NextSeq/PCR-based","BGI","ONT")
levels(dt$PLATFORM) <- c("HiSeqXTen/PCR-free","BGI","ONT","HiSeq2500/PCR-based","HiSeq2500/PCR-based","NextSeq/PCR-based","NextSeq/PCR-based")
dt <- dt[Locus %in% c("A","B","C","DQB1")]

call <- dt[,.(Locus,Chromosome,Allele,Expected_Coverage,FILE)]
call <- data.table(dcast(call,Locus+FILE+Expected_Coverage~Chromosome,value.var = "Allele"))
names(call) <- c("Locus","FILE","Expected_Coverage","Allele1","Allele2")
call <- merge(call,ref,by=c("Locus"),all.x = TRUE)
call[,Allele1:=str_remove(Allele1,"G")]
call[,Allele2:=str_remove(Allele2,"G")]

call_long <- melt(call,id.vars = c("Locus","FILE","Expected_Coverage"),variable.name = "Label",value.name = "Allele")
call_long[,"V1":= tstrsplit(Allele,"*",fixed=TRUE)[2]]
call_long[,c("V1","V2","V3","V4") := tstrsplit(V1,":",fixed=TRUE)]
call_long[Locus=="C",V3:="01"]
call_long[,Low:=V1]
call_long[,Intermediate:=paste(V1,V2,sep=":")]
call_long[,High:=paste(V1,V2,V3,sep=":")]

low <- data.table(dcast(call_long,Locus+FILE+Expected_Coverage~Label,value.var = "Low"))
low[,STATUS:=ifelse(setequal(c(Allele1,Allele2), c(Ref_Allele1,Ref_Allele2)),1,0),by=list(Locus,FILE,Expected_Coverage)]
ggplot(low[Expected_Coverage>=1],aes(as.factor(Expected_Coverage),FILE,fill=as.factor(STATUS)))+
  geom_tile(color="white")+
  scale_fill_manual(values=c("0"="red","1"="blue"))+
  scale_x_discrete(labels=as.character(c(5,rep("",3),25,rep("",4),50,rep("",4),75,rep("",4),100)))+
  facet_wrap(~Locus,nrow=1)+
  theme(aspect.ratio = 1)+
  guides(fill=FALSE)+
  labs(x="Expected coverage",y=NULL)

intermediate <- data.table(dcast(call_long,Locus+FILE+Expected_Coverage~Label,value.var = "Intermediate"))
intermediate[,STATUS:=ifelse(setequal(c(Allele1,Allele2), c(Ref_Allele1,Ref_Allele2)),1,0),by=list(Locus,FILE,Expected_Coverage)]
ggplot(intermediate[Expected_Coverage>=1],aes(as.factor(Expected_Coverage),FILE,fill=as.factor(STATUS)))+
  geom_tile(color="white")+
  scale_fill_manual(values=c("0"="red","1"="blue"))+
  scale_x_discrete(labels=as.character(c(5,rep("",3),25,rep("",4),50,rep("",4),75,rep("",4),100)))+
  facet_wrap(~Locus,nrow=1)+
  theme(aspect.ratio = 1)+
  guides(fill=FALSE)+
  labs(x="Expected coverage",y=NULL)

dev.off()  

# high <- data.table(dcast(call_long,Locus+FILE+Expected_Coverage~Label,value.var = "High"))
# high[,STATUS:=ifelse(setequal(c(Allele1,Allele2), c(Ref_Allele1,Ref_Allele2)),1,0),by=list(Locus,FILE,Expected_Coverage)]
# ggplot(high[Expected_Coverage>=1],aes(as.factor(Expected_Coverage),FILE,fill=as.factor(STATUS)))+
#   geom_tile(color="white")+
#   scale_fill_manual(values=c("0"="red","1"="blue"))+
#   scale_x_discrete(labels=as.character(c(1,"",3,"",5,rep("",3),25,rep("",4),50,rep("",4),75,rep("",4),100)))+
#   facet_wrap(~Locus,nrow=1)+
#   theme(aspect.ratio = 1)+
#   guides(fill=FALSE)+
#   labs(x="Expected coverage",y=NULL)