## Total R DGRP script no markdown

version


list.of.packages <- c("ggplot2", "dplyr","ggplot2","tidyr","qqman","lme4","outliers")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
library(data.table)
library(tidyr)
library(dplyr)
library(ggplot2)
library(lme4)
library(qqman)
library(outliers)



setwd("~")
#unzip(zipfile="denecke.et.al.2017.raw.data.zip",exdir="denecke.et.al.2017.raw.data")
setwd("./denecke.et.al.2017.raw.data")

source("./functions.list.R")
source("./themes.R")
doses <- c("25ppm","100ppm")

##########Section 3.1: Analysis of the Distribution of Imidacloprid_Resistance##########
dgrp.wiggle.raw <- fread("DGRP_raw_wiggle_data.csv",sep=",",header=T)
dir.create("./Section_3.1_Analysis_of_the_Distribution_of_Imidacloprid_Resistance")
setwd("./Section_3.1_Analysis_of_the_Distribution_of_Imidacloprid_Resistance")



#Summarize your raw wiggle data
sum.wiggle.60 <- dgrp.wiggle.raw %>% wiggle.parse(write=F) %>% 
      wiggle.rmr(write=F) %>% filter(time==60) %>%  group_by(genotype,dose) %>% 
      summarize(rmr.ci=conf.diff(rmr),rmr=mean(rmr),raw=mean(raw)) %>% 
      arrange (desc(dose),rmr)

# Add some modifications that will help for plotting Figure 3.1
sum.wiggle.60$rank <- factor(as.character(rep(1:(nrow(sum.wiggle.60)/2),2)),
                             levels=unique(rep(1:(nrow(sum.wiggle.60)/2),2)))
sum.wiggle.60$dose <- factor(sum.wiggle.60$dose,levels=c("25ppm","100ppm"))
sum.wiggle.60$colour <- rep(colorRampPalette(c("steelblue","orangered"))(nrow(sum.wiggle.60)/2),2)
annot <- data.frame(x = rep(-Inf, length(unique(sum.wiggle.60$dose))), 
                    y = rep(Inf, length(unique(sum.wiggle.60$dose))), 
                    dose=as.character(unique(sum.wiggle.60$dose)),
                    labs=LETTERS[1:length(unique(sum.wiggle.60$dose))])

##### Plot Figure 3.1: Phenotypic Spread of the DGRP
pdf("Figure_1_Phenotypic_Spread_of_the_DGRP.pdf",width=12) 
      gp <- ggplot(data=sum.wiggle.60,aes(x=rank,y=rmr,fill=rank),stat="identity")
      gp <- gp+facet_grid(~dose)
      gp <- gp+geom_bar(position="dodge",stat="identity",colour="black",size=.1)
      gp <- gp+geom_errorbar(aes(ymin=rmr+rmr.ci,ymax=rmr-rmr.ci),width=.8,size=.3)
      gp <- gp+scale_fill_manual(values=sum.wiggle.60$colour)
      gp <- gp+scale_y_continuous(expand = c(0, 0),limits=c(0,1.3),
                                  breaks=c(seq(0,1.2,.2)),labels=c(seq(0,1.2,.2)))
      gp <- gp+labs(y="Mean RMR Value\n",x="\nOrdered Genotype")
      gp <- gp+geom_text(size=7,data=annot,aes(x, y, label=labs,fontface="bold",family="serif"), inherit.aes=FALSE,vjust=1.2,hjust=-.2)
      gp <- gp+theme_bw(base_size=14,base_family="serif")
      gp <- gp+spread_theme
      print(gp)
dev.off()


##### Plot Figure S1: Correlation of 25 and 100ppm Imidacloprid Responses
sum.wiggle.60 <- arrange(sum.wiggle.60,genotype)
dose.compare <- data.frame(genotype=unique(sum.wiggle.60$genotype),
                           high.dose=filter(sum.wiggle.60,dose=="100ppm")$rmr,
                           low.dose=filter(sum.wiggle.60,dose=="25ppm")$rmr)
pdf("Figure_S1_Correlation_of_25_and_100ppm_Imidacloprid_Responses.pdf")
gp <- ggplot(dose.compare,aes(x=high.dose,y=low.dose))
gp <- gp+geom_point(size=2)
gp <- gp+geom_smooth(method="lm")
gp <- gp+labs(y="100ppm RMR Value\n",x="\n25ppm RMR Value")
gp <- gp+theme_bw(base_size=14,base_family="serif")
gp <- gp+corr_theme
print(gp) 
dev.off()


##### Plot Figure S2: Correlation of Motility in the Absence of Imidacloprid with RMR Value

sum.wiggle <- dgrp.wiggle.raw %>% wiggle.parse() %>% wiggle.rmr() %>% group_by(genotype,time,dose) %>% 
      summarize(rmr=mean(rmr),raw=mean(raw)) %>% arrange (desc(time),genotype)

initial.compare <- sum.wiggle %>% gather(mes,output,rmr,raw) %>% filter((time=="0" & mes=="raw") | (time=="60" & mes=="rmr")) 
initial.compare <- data.frame(genotype=initial.compare$genotype,initial.motility=filter(initial.compare,mes=="raw")$output,wiggle.response=filter(initial.compare,mes=="rmr")$output)

pdf("Figure S2: Correlation_of_Initial_Motility_and_RMR_Value.pdf")
gp <- ggplot(initial.compare,aes(x=initial.motility,y=wiggle.response))
gp <- gp+geom_point(size=2)
gp <- gp+geom_smooth(method="lm")
gp <- gp+ylab("Wiggle Response\n")
gp <- gp+xlab("\nInitial Motility")
gp <- gp+theme_bw(base_size=14,base_family="serif")
gp <- gp+corr_theme
print(gp) 
dev.off()

##### Calculate Heritability

corrected.wiggle.data <- dgrp.wiggle.raw %>% wiggle.parse() %>% wiggle.rmr()  %>% filter(time==60)
heritability.list=list()
for (a in doses){
      heritability.list[[a]]=list()
      sub.dose <- filter(corrected.wiggle.data,dose==a)
      sub.rmr <- sub.dose$rmr
      test <- lmer(sub.rmr~1+(1|sub.dose$genotype))
      heritability.list[[a]] <- c(genetic.variance=NA,environmental.variance=NA,heritability=NA)
      heritability.list[[a]]["genetic.variance"] <- VarCorr(test)$'sub.dose$genotype'[1]
      heritability.list[[a]]["environmental.variance"] <- heritability.list[[a]]["genetic.variance"]+attr(VarCorr(test),"sc")^2
      heritability.list[[a]]["heritability"] <- heritability.list[[a]]['environmental.variance']/(heritability.list[[a]]["environmental.variance"]+heritability.list[[a]]["genetic.variance"])
} 

write.csv(heritability.list,file="Heritability_List.csv")


## Clean phenotype data to submit
for(d in unique(sum.wiggle.60$dose)){
      single.dose <- sum.wiggle.60 %>% filter(dose==d) %>% select(genotype,rmr)
      single.dose[['genotype']] <-gsub("RAL_","",single.dose[['genotype']]) 
      clean.matrix <- matrix(as.numeric(as.matrix(single.dose)),ncol=2)
      write.table(clean.matrix,col.names=F,row.names=F,sep=",",file=paste("GWAS_Submit_Mackay_Pipeline",d,".csv",sep=""))                  
}             


setwd("..")



##########Section 3.2: A GWAS for Imidacloprid Response Yields Many Neurological Candidate Genes ##############################################
gwas.all.associations.25ppm <- fread("gwas.all_mean_25ppm.assoc")
gwas.all.associations.100ppm <- fread("gwas.all_mean_100ppm.assoc")

gwas.top.candidates.25ppm <- readLines("gwas.top_25ppm.annot")
gwas.top.candidates.100ppm <- readLines("gwas.top_100ppm.annot")

mod.encode.raw <- fread("modencode.tsv",sep="\t",header=T)

dir.create("./Section_3.2_A_GWAS_for_Imidacloprid_Response_Yields_Many_Neurological_Candidate_Genes")
setwd("./Section_3.2_A_GWAS_for_Imidacloprid_Response_Yields_Many_Neurological_Candidate_Genes")


##### Generate Manhattan and QQplots for Figures S3 and S4
assoc.list <- vector('list',2)
names(assoc.list) <- c("25ppm","100ppm")

for (d in doses){  
      all.associations <- get(paste("gwas.all.associations.",d,sep=""))
      all.associations <- select(all.associations,ID,SinglePval)
      variant <- unlist(strsplit(all.associations$ID,split="_"))
      chromosome <- variant[seq(1, length(variant), 3)]
      base <- as.numeric(variant[seq(2, length(variant), 3)])
      chromosome <- multi.gsub(chromosome,find=c("X",'2L','2R','3L','3R'),replace=c(1,2,3,4,5))
      gwas.plot.data <- all.associations %>% mutate(CHR=as.numeric(chromosome),BP=base) %>% select(SinglePval,ID,CHR,BP)
      colnames(gwas.plot.data)=c('P','SNP','CHR','BP')
      variant.sig <- gwas.plot.data$SNP[which(gwas.plot.data$P<=.05/length(gwas.plot.data$SNP))]
      gwas.plot.data <- filter(gwas.plot.data,P<.0001)
      assoc.list[[d]] <- gwas.plot.data
}


par(mar=c(5,6,2,1),mgp=c(3,1,0),mfrow=c(1,2))
manhattan(assoc.list[["25ppm"]],col=c("black","orange"),genomewideline=-log10(.05/length(assoc.list[["25ppm"]]$SNP)),chrlabs=c("X",'2L','2R','3L','3R'),
          highlight=variant.sig,main="25ppm",cex.main=2,ylim=c(0,-log10(1e-8)),family="serif",font=2,cex.lab=1.5,cex.axis=1.5,
          xlab="",ylab=list(expression(bold(paste("-log"[10]*"(p)"))),font=2))  

manhattan(assoc.list[["100ppm"]],col=c("black","orange"),genomewideline=-log10(.05/length(assoc.list[["100ppm"]]$SNP)),chrlabs=c("X",'2L','2R','3L','3R'),
          highlight=variant.sig,main="100ppm",cex.main=2,ylim=c(0,-log10(1e-8)),family="serif",font=2,cex.lab=1.5,cex.axis=1.5,
          xlab="",ylab=list(expression(bold(paste("-log"[10]*"(p)"))),font=2))

par(mar=c(5,6,2,1),mgp=c(3,1,0),mfrow=c(1,2),bty="l") 

qq(assoc.list[["25ppm"]]$P, main="25ppm", family="serif",cex.main=2,font.axis=2,cex.lab=1.5,cex.axis=1.5,
   ylab=list(expression(bold(paste("Observed -log"[10]*"(p)"))),font=2), 
   xlab=list(expression(bold(paste("Expected -log"[10]*"(p)"))),font=2))

qq(assoc.list[["100ppm"]]$P, main="100ppm", family="serif",cex.main=2,font.axis=2,cex.lab=1.5,cex.axis=1.5,
   ylab=list(expression(bold(paste("Observed -log"[10]*"(p)"))),font=2), 
   xlab=list(expression(bold(paste("Expected -log"[10]*"(p)"))),font=2)) 




##### Generate Table 1

low.dose.gwas <- dgrp.gwas.parse(gwas.top.candidates.25ppm) %>% mutate(dose="25ppm")
high.dose.gwas <- dgrp.gwas.parse(gwas.top.candidates.100ppm) %>% mutate(dose="100ppm")
write.csv(rbind(low.dose.gwas,high.dose.gwas),file="Table_1_GWAS_summary_output.csv")


##### Calculate enrichement in L3 CNS
unique.gwas.list <- unique(c(low.dose.gwas$Gene.Name,high.dose.gwas$Gene.Name))
unique.gwas.list <- unique.gwas.list[unique.gwas.list!="intergenic"]
l3.enrichment <- mod.encode.enrich(raw.data.file=mod.encode.raw,
                  gene.list = unique.gwas.list, 
                  parent.library="modENCODE_mRNA-Seq_tissues",
                  stage.search.string="L3",
                  tissue.search.string="CNS")
write.csv(l3.enrichment,"GWAS_L3_CNS_Enrichment.csv")


rm(list=c("mod.encode.raw","gwas.all.associations.25ppm","gwas.all.associations.100ppm"))

setwd("..")
##########Section 3.3: A Transcriptome Wide Association Study Suggests Roles for Cyp6g1 and Cyp6g2 ########################################################

## Read Relevant Data into R
g1.haplotypes <- read.csv("cyp6g1_haplotypes.csv",header=T)
dgrp.male.transcriptome <- fread("dgrp.array.exp.male.txt",sep=" ",header=T) ## Available from http://dgrp2.gnets.ncsu.edu/data.html
dgrp.female.transcriptome <- fread("dgrp.array.exp.female.txt",sep=" ",header=T) ## Available from http://dgrp2.gnets.ncsu.edu/data.html
fbgid.symbol.key <- fread("GeneID_Symbol_Key.csv" ,sep=",",header=T)
dir.create("Section_3.3_A_Transcriptome_Wide_Association_Study_Suggests_Roles_for_Cyp6g1_and_Cyp6g2")
setwd("Section_3.3_A_Transcriptome_Wide_Association_Study_Suggests_Roles_for_Cyp6g1_and_Cyp6g2")

##### Perform TWAS and generate Table 2
combined.transcriptomes <- transcript.combine(male=dgrp.male.transcriptome,female=dgrp.female.transcriptome)
combined.transcriptomes <- merge(combined.transcriptomes,fbgid.symbol.key,by="gene")#full_join(combined.transcriptomes,fbgid.symbol.key,by.x="gene")

clean.transcriptome.data <- combined.transcriptomes %>% 
      gather(key="Genotype",value="RPKM",RAL_21:RAL_913) %>%
      group_by(symbol,Genotype) %>% 
      summarize(avg=mean(RPKM,na.rm=T)) %>%
      spread(key=Genotype,value=avg)

# Subset the genotypes which are in both the phenotyic and genotypic datasets
transcriptome.genotypes <- colnames(clean.transcriptome.data)[-1]
wiggle.genotypes <- as.character(sum.wiggle.60$genotype)
common.genotypes <- intersect(transcriptome.genotypes,wiggle.genotypes)
phenotype.data.subset <- sum.wiggle.60[wiggle.genotypes %in% common.genotypes,]
transcriptome.subset.symbol <- data.table(data.frame(clean.transcriptome.data)[,c(TRUE,transcriptome.genotypes %in% common.genotypes)])

## Run a linear model predicting RMR value with each transcript in the transcriptome
symbol <- transcriptome.subset.symbol$symbol
transcriptome.subset <- select(transcriptome.subset.symbol,-symbol)
transcriptome.test <- sapply(select(transcriptome.subset,starts_with("RAL")),as.numeric)
tab.2 <- vector('list',2)
names(tab.2) <- doses
for (d in doses){
      sub.phen.dose <- as.numeric(subset(phenotype.data.subset,dose==d)$rmr)
      output=matrix(NA,nrow=length(transcriptome.test[,1]),ncol=3)
      colnames(output)=c("estimate","pval","adjrsq")
      for (i in 1:nrow(transcriptome.test)) {
            test <- summary(lm(sub.phen.dose~transcriptome.test[i,]))
            est <- test$coefficients[2,1]
            pval <- test$coefficients[2,4]
            rsq <- test$adj.r.squared
            output[i,] <- c(est,pval,rsq)
      }
      
      rownames(output)=symbol
      twas.summary <- data.frame(output) %>%
            mutate(gene=symbol,dose=d) %>%
            filter(pval<1e-3) %>% arrange(pval) %>%
            select(dose,gene,pval,estimate,adjrsq)
      colnames(twas.summary) <- c("dose","gene","p-value","effect size","adjusted R squared")
      table_2[[d]] <- twas.summary
}
rbindlist(table_2,use.names=T)
write.csv(rbindlist(table_2),file="Table_2_TWAS_Candidates.csv")

##### Create Figure 2: Cyp6g1 and Cyp6g2 Correlates with Imidacloprid Response
# Create data frame contaning RMR values and P450 expression
g1.haplotype.sub <- g1.haplotypes %>% filter(genotype%in%common.genotypes)
individual.transcripts <- transcriptome.subset.symbol %>% filter(symbol=="Cyp6g1"|symbol=="Cyp6g2") %>%
      gather(genotype, val, 2:(length(common.genotypes)+1)) %>%  
      spread(symbol, val) %>% arrange(genotype) %>% group_by(genotype) %>% summarise_each(funs(sum(., na.rm=TRUE))) %>%
      full_join(select(filter(phenotype.data.subset,dose=="25ppm"),genotype,rmr),by="genotype") %>% 
      full_join(g1.haplotype.sub,by="genotype") %>% gather(key="transcript",value="expression",Cyp6g1,Cyp6g2)

# Prepare data frame for plotting
individual.transcripts$transcript <- factor(individual.transcripts$transcript,levels=c("Cyp6g1","Cyp6g2"))
equation.list <- vector("character",2)
for (t in unique(individual.transcripts$transcript)) equation.list[which(t==unique(individual.transcripts$transcript))] <- equation.extract("rmr","expression",individual.transcripts,"transcript",t)
annot <- data.frame(x = rep(-Inf, length(unique(individual.transcripts$transcript))), 
                    y = rep(Inf, length(unique(individual.transcripts$transcript))), 
                    transcript=as.character(unique(individual.transcripts$transcript)),
                    equation=equation.list,
                    labs=LETTERS[1:length(unique(individual.transcripts$transcript))])
annot$transcript <- factor(annot$transcript,levels=annot$transcript)

#Plot Figure 2: Cyp6g1 and Cyp6g2 Correlates with Imidacloprid Response
pdf(file="Figure_2_Cyp6g1_and_Cyp6g2_Correlates_with_Imidacloprid_Response",width=16)
gp <- ggplot(individual.transcripts,aes(x=expression,y=rmr))
gp <- gp+facet_grid(facets=.~transcript,scales="free_x")
gp <- gp+geom_point(size=4,aes(shape=haplotype,colour=haplotype))
gp <- gp+scale_shape_identity()
gp <- gp+geom_smooth(method="lm")
gp <- gp+labs(y="RMR Value\n",x="\nTranscript Expression Expression (log2) (RPKM)")
gp <- gp+scale_colour_manual(values=c("red2","blue","green4","black"))
gp <- gp+ylim(c(0,1.2))
gp <- gp+geom_text(size=7,data=annot,aes(x, y, label=labs,fontface="bold",family="serif"), inherit.aes=FALSE,vjust=1.2,hjust=-.2)
gp <- gp+geom_text(size=4,data=annot,aes(Inf, Inf, label=equation,family="serif"), inherit.aes=FALSE,vjust=1.2,hjust=1)
gp <- gp+theme_bw(base_size=14,base_family="serif")
gp <- gp+cyp6g1.ind_theme
print(gp) 
dev.off()

rm(list=c("dgrp.male.transcriptome","dgrp.female.transcriptome"))

setwd("..")
########################Section 3.4 Imidacloprid and Metabolite Quantification in the DGRP ####

dgrp.hplc.data <- fread("DGRP_HPLC_data.csv",sep=",",header=T)

dir.create("3.4_DGRP_HPLC")
setwd("./3.4_DGRP_HPLC")


dgrp.hplc.data <- dgrp.hplc.data %>% gather("mes.type","mes.value",larvae,media) %>%
      filter(!is.na(mes.value)) %>%
      arrange(rank) 

dgrp.hplc.summary <- dgrp.hplc.data %>% group_by(genotype,metabolite,mes.type) %>%
      summarize(mes.value=mean(mes.value)) 
wiggle.hplc.subset <- sum.wiggle.60 %>% filter(dose=="25ppm") %>% filter(genotype %in% unique(dgrp.hplc.summary$genotype)) %>% select(genotype,rmr)

dgrp.hplc.summary <- full_join(wiggle.hplc.subset,dgrp.hplc.summary,by="genotype")

dgrp.hplc.plot.data <- dgrp.hplc.data[-intersect(which(dgrp.hplc.data$metabolite=="Imidacloprid"),which(dgrp.hplc.data$mes.type=="media")),]
dgrp.hplc.plot.data$genotype <- factor(dgrp.hplc.plot.data$genotype,levels=unique(dgrp.hplc.plot.data$genotype))
dgrp.hplc.plot.data$metabolite <- factor(dgrp.hplc.plot.data$metabolite,levels=c("Imidacloprid","5-OH","Olefin"))
dgrp.hplc.plot.data$mes.type <- factor(dgrp.hplc.plot.data$mes.type)

vars <- data.frame(expand.grid(levels(dgrp.hplc.plot.data$metabolite), levels(dgrp.hplc.plot.data$mes.type)))
colnames(vars) <- c("metabolite", "mes.type")

annot <- data.frame(x = rep(Inf, dim(vars)[1]), 
                    y = rep(Inf, dim(vars)[1]), 
                    vars,
                    labs=LETTERS[1:dim(vars)[1]])
equation.list <- vector("character",6)
for(m in unique(dgrp.hplc.plot.data$metabolite)){
      for(t in (unique(dgrp.hplc.plot.data$mes.type))){
            sub <- subset(dgrp.hplc.summary,metabolite==m & mes.type==t)
            equation.list[intersect(which(annot$metabolite==m),which(annot$mes.type==t))] <- equation.extract.hplc("mes.value","rmr",data=sub,sub=F)
      }}
equation.list[intersect(which(annot$metabolite=="Imidacloprid"),which(annot$mes.type=="media"))] <- ""
annot <- mutate(annot,equation.list=equation.list)

dot.col <- c(rep("steelblue2",180),rep("orangered2",175))
pdf(file=paste("Figure 3: DGRP_HPLC.pdf",sep="."))
gp <- ggplot(dgrp.hplc.plot.data,aes(x=genotype,y=mes.value)) 
gp <- gp+facet_grid(metabolite~mes.type,scales="free_y")
gp <- gp+geom_dotplot(binaxis='y', stackdir='center',colour="black",dotsize=.5,fill=dot.col) 
gp <- gp+labs(y="Quantity of Chemical",x="Genotype")
gp <- gp+annotate("rect",xmin=-Inf,xmax=9.5,ymin=-Inf,ymax=Inf,fill="steelblue2",alpha=.3)
gp <- gp+annotate("rect",xmin=9.5,xmax=Inf,ymin=-Inf,ymax=Inf,fill="orangered2",alpha=.3)
gp <- gp+geom_dotplot(binaxis='y', stackdir='center',colour="black",dotsize=.9,fill=dot.col) 
gp <- gp+stat_summary(geom="errorbar",colour='black',size=.75,width=.6,fun.data="mean_se")
gp <- gp+stat_summary(fun.y="mean",fun.ymax="mean",fun.ymin="mean", geom="crossbar", width=.6,size=.3)
gp <- gp+geom_text(size=7,data=annot,aes(x, y, label=labs,fontface="bold",family="serif"),vjust=1.2,hjust=1.2)
gp <- gp+geom_text(size=7,data=annot,aes(9.5, y, label=equation.list,fontface="bold",family="serif"),vjust=1.2)
gp <- gp+theme_bw(base_size=14,base_family="serif")
gp <- gp+dgrp.hplc_theme
print(gp)
dev.off()

setwd("..")





###############################Section 3.5: Cyp6g1KO ######################################

cyp6g1.ko.wiggle <- read.csv("total_cyp6g1_KO_wiggle.csv")
cyp6g1KO.hplc.data <- fread("RAL517KO_HPLC_data.csv",sep=",",header=T,data.table=F)

dir.create("3.5_Cyp6g1KO")
setwd("./3.5_Cyp6g1KO")

##### Generate Figure 4: Imidacloprid Resistance in various Cyp6g1KO backgrounds 

#Prepare data frame for plotting figure 4
cyp6g1.ko.wiggle.plot <- cyp6g1.ko.wiggle %>% filter(time==60) %>% select(genotype,dose,time,rmr,identifier)
cyp6g1.ko.wiggle.plot$genotype <- factor(cyp6g1.ko.wiggle.plot$genotype,levels=c("RAL_517","RAL_517-Cyp6g1KO","Canton-S","Canton-S-Cyp6g1KO","Wxac","Wxac-Cyp6g1KO"))
cyp6g1.ko.wiggle.plot$identifier <- factor(cyp6g1.ko.wiggle.plot$identifier,levels=c("RAL_517 (25ppm)","RAL_517 (100ppm)","Canton-S (5ppm)","Wxac (5ppm)"))
vars <- data.frame(expand.grid(levels(cyp6g1.ko.wiggle.plot$identifier)))
colnames(vars) <- c("identifier")
annot <- data.frame(x = rep(Inf, dim(vars)[1]), y = rep(Inf, dim(vars)[1]), vars,labs=LETTERS[1:dim(vars)[1]])

equation.list <- vector("character",4)
for(i in unique(levels(cyp6g1.ko.wiggle.plot$identifier))){
            sub <- subset(cyp6g1.ko.wiggle.plot,identifier==i)
            equation.list[which(annot$identifier==i)] <- equation.extract.hplc("rmr","genotype",data=sub,sub=F,bonf=4)
      }
annot$pval <- equation.list


#Plot Figure 4

pdf(file="Figure 4: Imidacloprid Resistance in various Cyp6g1KO backgrounds.pdf",width=10,height=5)
      gp <- ggplot(cyp6g1.ko.wiggle.plot,aes(x=genotype,y=rmr,fill=genotype))
      gp <- gp+facet_grid(.~identifier,scales="free_x")
      gp <- gp+geom_dotplot(binaxis='y',stackdir="center",dotsize=.75,binwidth=.05)
      gp <- gp+stat_summary(geom="errorbar",colour='black',size=1,width=.3,fun.data="mean_cl_normal")
      gp <- gp+stat_summary(fun.y="mean",fun.ymax="mean",fun.ymin="mean", geom="crossbar", width=.4,colour='black',size=.8)
      gp <- gp+scale_fill_manual(values=c("blue4","blue1","purple4","purple1","gold4","gold1")) 
      gp <- gp+labs(y="RMR Value\n",x="\nGenotype")
      gp <- gp+ylim(c(0,1.4))
      gp <- gp+geom_text(size=7,data=annot,aes(x, y, label=labs,fontface="bold",family="serif"),inherit.aes=F,vjust=1.2,hjust=1.2)
      gp <- gp+geom_text(size=7,data=annot,aes(1.5, y, label=pval,fontface="bold",family="serif"),inherit.aes=F,vjust=1.2)
      gp <- gp+theme_bw(base_size=14,base_family="serif")
      gp <- gp+g1KO.wiggle_theme
print(gp)
dev.off()




######Generate Figure 5: The HPLC profile of RAL517-Cyp6g1KO

#Prepare data frame
cyp6g1KO.hplc.data.plot <- cyp6g1KO.hplc.data[-intersect(which(cyp6g1KO.hplc.data$metabolite=="Imidacloprid"),which(cyp6g1KO.hplc.data$measurement_type=="Media")),]
cyp6g1KO.hplc.data.plot$genotype <- factor(cyp6g1KO.hplc.data.plot$genotype,levels=c("RAL_517","RAL_517-Cyp6g1KO"))
cyp6g1KO.hplc.data.plot$metabolite <- factor(cyp6g1KO.hplc.data.plot$metabolite,levels=c("Imidacloprid","5-OH","Olefin"))
vars <- data.frame(expand.grid(levels(cyp6g1KO.hplc.data.plot$metabolite), levels(as.factor(cyp6g1KO.hplc.data.plot$measurement_type))))
colnames(vars) <- c("metabolite", "measurement_type")
annot <- data.frame(x = rep(-Inf, dim(vars)[1]), y = rep(Inf, dim(vars)[1]), vars,labs=LETTERS[1:dim(vars)[1]])

equation.list <- vector("character",6)
for(m in unique(cyp6g1KO.hplc.data$metabolite)){
      for(t in (unique(cyp6g1KO.hplc.data$measurement_type))){
            sub <- subset(cyp6g1KO.hplc.data,metabolite==m & measurement_type==t)
            equation.list[intersect(which(annot$metabolite==m),which(annot$measurement_type==t))] <- equation.extract.hplc("measurement_value","genotype",data=sub,sub=F,bonf=6)
      }}
equation.list[intersect(which(annot$metabolite=="Imidacloprid"),which(annot$mes.type=="media"))] <- ""
annot <- mutate(annot,pval=equation.list)

#Plot Figure 5
pdf(file="Figure 5: HPLC Profile of RAL517-Cyp6g1KO.pdf",height=10,width=6)
      gp <- ggplot(cyp6g1KO.hplc.data.plot,aes(x=genotype,y=measurement_value,fill=genotype))
      gp <- gp+facet_grid(metabolite~measurement_type,scales="free_y")
      gp <- gp+geom_dotplot(binaxis='y', stackdir='center',dotsize=2) 
      gp <- gp+labs(y="Quantity of Chemical\n",x="\nGenotype")
      gp <- gp+stat_summary(geom="errorbar",colour='black',size=1,width=.4,fun.data="mean_se")
      gp <- gp+stat_summary(fun.y="mean",fun.ymax="mean",fun.ymin="mean", geom="crossbar", width=.4,size=.3)
      gp <- gp+scale_fill_manual(values=c("blue4","blue1"))
      gp <- gp+geom_text(size=7,data=annot,aes(x, y, label=labs,fontface="bold",family="serif"),inherit.aes=FALSE,vjust=1.2,hjust=-.2)
      gp <- gp+geom_text(size=7,data=annot,aes(1.5, y, label=pval,fontface="bold",family="serif"),inherit.aes=F,vjust=1.2)
      gp <- gp+theme_bw(base_size=14,base_family="serif")
      gp <- gp+g1KO.wiggle_theme
print(gp)
dev.off()


g1.ko.hplc.sig.list <- vector('list')
for (t in unique(cyp6g1KO.hplc.data$measurement_type)){ 
      for (m in unique(cyp6g1KO.hplc.data$metabolite)){
            sub <- subset(cyp6g1KO.hplc.data,measurement_type==t & metabolite==m) 
            #mod <- with(sub,aov(measurement_value~genotype))
            g1.ko.hplc.sig.list[[paste(t,m)]] <- with(sub,t.test(measurement_value~genotype)$p.val)
                  #TukeyHSD(x=mod, which='genotype', conf.level=0.95)$genotype[,"p adj"]
      }}
g1.ko.hplc.sig.list <- unlist(g1.ko.hplc.sig.list)*length(g1.ko.hplc.sig.list)
write.csv(g1.ko.hplc.sig.list,"T Test Bonferroni for Figure 5.csv")
setwd("..")


################Begin Section 6 Cyp6g2 is enriched in the Digestive Tissues of Resistant Lines and Metabolizes Imidacloprid#######################

uas.cyp6g2.hplc.data <- fread("UAS_cyp6g2_HPLC_data.csv",sep=",",header=T)
P450.real.time.data <- fread("P450_real_time_data.csv",sep=",",header=T)
uas.cyp6g12.wiggle <- read.csv("cyp6g12_wiggle.data.csv")
dir.create("Section_3.6_Cyp6g2")
setwd("./Section_3.6_Cyp6g2")

##### Generate Figure 6: Expression of Cyp6g1 and Cyp6g2 in the Digestive Tissue
parsed.rt.data <- rt.data.parse(data=P450.real.time.data)

#Remove biological samples 3 and 4 from RAL_360 due to abberations in the housekeeper expression levels
filtered.data <- parsed.rt.data[which(!((parsed.rt.data$biological.rep=="3" | parsed.rt.data$biological.rep=="4") & parsed.rt.data$Sample=="360")),]

sum.rt.data <- filtered.data %>% 
      group_by(Target,Sample,biological.rep) %>% 
      summarize(mean.ct=mean(Cq)) %>% 
      spread(key=Target,value=mean.ct) %>% 
      mutate(g1dct=-(RP49-CYP6G1),g2dct=-(RP49-CYP6G2),t3dct=-(RP49-CYP6T3)) %>% 
      mutate(genotype=as.character(paste("RAL_",Sample,sep=""))) %>%
      mutate(haplotype=allele.ad(genotype)) %>% data.table() %>%
      mutate(g1ddct=g1dct-mean(.[which(.[["genotype"]]=="RAL_776"),][["g1dct"]],na.rm=T)) %>% 
      mutate(g2ddct=g2dct-mean(.[which(.[["genotype"]]=="RAL_776"),][["g2dct"]],na.rm=T)) %>% 
      mutate(t3ddct=t3dct-mean(.[which(.[["genotype"]]=="RAL_776"),][["t3dct"]],na.rm=T)) %>% 
      mutate(Cyp6g1=2^-g1ddct,Cyp6g2=2^-g2ddct,Cyp6t3=2^-t3ddct) %>% 
      mutate(type="Digestive Tissues Larvae") %>%
      select(genotype,haplotype,type,Cyp6g1:Cyp6t3) %>%
      gather(key="gene",value="relative.expression",Cyp6g1:Cyp6t3)

dgrp.p450.trans <- combined.transcriptomes %>%
      select(symbol,RAL_138,RAL_360,RAL_776,RAL_843) %>%
      filter(symbol %in% c("Cyp6g1","Cyp6g2","Cyp6t3")) %>% 
      transpose.collapse(column="symbol") %>%
      mutate(haplotype=allele.ad(genotype)) %>%
      mutate(Cyp6g1=2^(Cyp6g1-mean(.[which(.[["genotype"]]=="RAL_776"),][["Cyp6g1"]],na.rm=T))) %>%
      mutate(Cyp6g2=2^(Cyp6g2-mean(.[which(.[["genotype"]]=="RAL_776"),][["Cyp6g2"]],na.rm=T))) %>%
      mutate(Cyp6t3=2^(Cyp6t3-mean(.[which(.[["genotype"]]=="RAL_776"),][["Cyp6t3"]],na.rm=T))) %>%
      mutate(type="Whole Body Adults") %>%
      select(genotype,haplotype,type,Cyp6g1:Cyp6t3) %>%
      gather(key="gene",value="relative.expression",Cyp6g1:Cyp6t3) 



combined.data <- rbindlist(list(sum.rt.data,dgrp.p450.trans),use.names=T)
combined.data <- data.frame(combined.data)[complete.cases(combined.data),]
combined.data$genotype <- factor(combined.data$genotype,levels=c("RAL_776","RAL_843","RAL_138","RAL_360"))
combined.data$haplotype <- factor(combined.data$haplotype,levels=c("M","AA"))

vars <- data.frame(expand.grid(unique(combined.data$gene), unique(combined.data$type)))
colnames(vars) <- c("gene", "type")

annot <- data.frame(x = rep(-Inf, dim(vars)[1]), y = rep(Inf, dim(vars)[1]), vars,labs=LETTERS[1:dim(vars)[1]])

pdf("Figure 6: Expression of Cyp6g1 and Cyp6g2 in the Digestive Tissue.pdf",width=12)
gp <- ggplot(data=combined.data,aes(x=genotype,y=relative.expression,fill=haplotype))
      gp <- gp+facet_grid(type~gene,scales="free_y")
      gp <- gp+geom_dotplot(binaxis='y',stackdir="center",dotsize=1)
      gp <- gp+annotate("rect",xmin=-Inf,xmax=2.5,ymin=-Inf,ymax=Inf,fill="steelblue2",alpha=.3)
      gp <- gp+annotate("rect",xmin=2.5,xmax=Inf,ymin=-Inf,ymax=Inf,fill="orangered2",alpha=.3)
      gp <- gp+geom_dotplot(binaxis='y',stackdir="center",dotsize=1.5)
      gp <- gp+stat_summary(geom="errorbar",colour='black',size=1,width=.4,fun.data="mean_se")
      gp <- gp+stat_summary(fun.y="mean",fun.ymax="mean",fun.ymin="mean", geom="crossbar", width=.4,size=.3)
      gp <- gp+scale_fill_manual(values=c("blue1","red1"))
      gp <- gp+labs(y="Expression Relative to RAL_776\n",x="\nGenotype")
      gp <- gp+geom_text(size=7,data=annot,aes(x, y, label=labs,fontface="bold",family="serif"),inherit.aes=F,vjust=1.2,hjust=-.2)
      gp <- gp+rt_theme
print(gp)
dev.off()

rt.sig.list <- list()
for(t in unique(combined.data$type)){
      for (g in unique(combined.data$gene)){
            sub <- subset(combined.data,type==t & gene==g)
            mod <- with(sub,aov(relative.expression~genotype))
            rt.sig.list[[paste(t,g)]] <- TukeyHSD(x=mod, which='genotype', conf.level=0.95)$genotype[,"p adj"]
      }}
write.csv(rt.sig.list,file="Table S3: Significance of qPCR Results.csv")
##### Generate Figure 7: Comparing Resistance Conferred by Cyp6g1 and Cyp6g2



uas.wiggle <- uas.cyp6g12.wiggle %>% wiggle.parse(write=F) %>% wiggle.rmr(write=F) %>% filter(time==60)

uas.wiggle$genotype <- multi.gsub(uas.wiggle$genotype,c("x"),c(" x "))
uas.wiggle$genotype <- factor(uas.wiggle$genotype,levels=c("HRG4 x UAS-NULL","HRG4 x UAS-Cyp6g1","HRG4 x UAS-Cyp6g2"))


pdf(file="Figure 7: Comparison of Cyp6g1 and Cyp6g2.pdf")
      gp <- ggplot(uas.wiggle,aes(x=genotype,y=rmr,fill=genotype))
      gp <- gp+facet_grid(.~dose,scales="fixed")
      gp <- gp+geom_dotplot(binaxis='y',stackdir="center",dotsize=.75,binwidth=.05)
      gp <- gp+stat_summary(geom="errorbar",colour='black',size=1,width=.3,fun.data="mean_cl_normal")
      gp <- gp+stat_summary(fun.y="mean",fun.ymax="mean",fun.ymin="mean", geom="crossbar", width=.4,colour='black',size=.8)
      gp <- gp+scale_fill_manual(values=c("black","red","blue")) 
      gp <- gp+labs(y="RMR Value\n",x="\nGenotype")
      gp <- gp+theme_bw(base_size=14,base_family="serif")
      gp <- gp+g1KO.wiggle_theme
      print(gp)
dev.off()
g12.sig.list <- list()
for (d in unique(uas.wiggle$dose)){
      sub <- subset(uas.wiggle,dose==d)
      mod <- with(sub,aov(rmr~genotype))
      g12.sig.list[[d]] <- TukeyHSD(x=mod, which='genotype', conf.level=0.95)$genotype[,"p adj"] 
      
}
write.csv(g12.sig,file="Table S4: Significance of Overexpression of Cyp6g1 and Cyp6g2.csv")
##### Generate Figure S5

uas.cyp6g2.hplc.data <- gather(uas.cyp6g2.hplc.data,key="measurement_type",value="measurement_value",larvae,media)
uas.cyp6g2.hplc.data <- uas.cyp6g2.hplc.data[complete.cases(uas.cyp6g2.hplc.data),]
plot.data <- uas.cyp6g2.hplc.data[-intersect(which(uas.cyp6g2.hplc.data$metabolite=="Imidacloprid"),which(uas.cyp6g2.hplc.data$measurement_type=="media")),]
plot.data$genotype <- factor(plot.data$genotype,levels=c("w1118","Cyp6g2-3a","Cyp6g2-3d"))
plot.data$metabolite <- factor(plot.data$metabolite,levels=c("Imidacloprid","5-OH","Olefin"))
plot.data$measurement_type <- factor(plot.data$measurement_type)
plot.data$measurement_value <- as.numeric(plot.data$measurement_value)

vars <- data.frame(expand.grid(levels(plot.data$metabolite), levels(plot.data$measurement_type)))
colnames(vars) <- c("metabolite", "measurement_type")
annot <- data.frame(x = rep(-Inf, dim(vars)[1]), y = rep(Inf, dim(vars)[1]), vars,labs=LETTERS[1:dim(vars)[1]])


equation.list <- vector("character",6)
for(m in unique(uas.cyp6g2.hplc.data$metabolite)){
      for(t in (unique(uas.cyp6g2.hplc.data$measurement_type))){
            sub <- subset(uas.cyp6g2.hplc.data,metabolite==m & measurement_type==t)
            equation.list[intersect(which(annot$metabolite==m),which(annot$measurement_type==t))] <- equation.extract.hplc("measurement_value","genotype",data=sub,sub=F,bonf=6)
      }}
equation.list[intersect(which(annot$metabolite=="Imidacloprid"),which(annot$mes.type=="media"))] <- ""
annot <- mutate(annot,pval=equation.list)

pdf(file="Figure S5: UAS-Cyp6g2 HPLC Output.pdf",height=10)
      gp <- ggplot(plot.data,aes(x=genotype,y=measurement_value,fill=genotype))
      gp <- gp+facet_grid(metabolite~measurement_type)
      gp <- gp+ylim(c(0,2200000))
      gp <- gp+geom_dotplot(binaxis='y', stackdir='center',dotsize=1.5) 
      gp <- gp+ylab("Quantity of Chemical")
      gp <- gp+stat_summary(geom="errorbar",colour='black',size=1,width=.4,fun.data="mean_se")
      gp <- gp+stat_summary(fun.y="mean",fun.ymax="mean",fun.ymin="mean", geom="crossbar", width=.4,size=.3)
      gp <- gp+geom_text(size=7,data=annot,aes(x, y, label=labs,fontface="bold",family="serif"),inherit.aes=FALSE,vjust=1.2,hjust=-.2)
      gp <- gp+geom_text(size=7,data=annot,aes(1.5, y, label=pval,fontface="bold",family="serif"),inherit.aes=F,vjust=1.2)
      gp <- gp+theme_bw(base_size=14,base_family="serif")
      gp <- gp+g1KO.wiggle_theme
print(gp)
dev.off()

setwd("..")
