wiggle.parse <- function(data,file.name="parsed.wiggle.data.csv",write=F){
      split.frame <- matrix(unlist(strsplit(as.character(data$Image.Name),"\\).")),ncol=8,byrow=T)[,c(1:6,8)]
      split.frame <- gsub("\\(","",split.frame)
      split.frame <- gsub(".tif","",split.frame)
      split.frame <- gsub("opRight","TopRight",split.frame)
      split.frame <- gsub("opLeft","TopLeft",split.frame)
      split.frame <- gsub("owRight","LowRight",split.frame)
      split.frame <- gsub("owLeft","LowLeft",split.frame)
      parsed.data <- cbind(split.frame,as.numeric(as.character(data$Wiggle.Index)))
      colnames(parsed.data) <- c("genotype","dose","time","plate","insecticide","date","position","Wiggle.Index")
      parsed.data <- data.frame(parsed.data)
      parsed.data$time <- as.numeric(gsub('min',"",parsed.data$time))
      parsed.data$Wiggle.Index <- as.numeric(as.character(parsed.data$Wiggle.Index))
      if(write==T){
            write.csv(x=parsed.data,file=file.name,row.names = F)
      }
      return(parsed.data) 
}





wiggle.rmr <- function(data,file.name="rmr.wiggle.data.csv",write=F){
      library(tidyr)
      library(dplyr)
      spread.data <- spread(data,time,Wiggle.Index)
      for(i in colnames(spread.data)[-1:-6]){
            rmr.values <- spread.data[[i]]/spread.data[["0"]]
            new.name <-paste(i,"rmr",sep="_")
            spread.data[[new.name]] <- rmr.values
      }
      raw.gather <- gather(spread.data,key="time",value="raw",
                           which(colnames(spread.data) %in% colnames(spread.data)[which(grepl("[0-9]",colnames(spread.data)) & !grepl("rmr",colnames(spread.data)))]))     
      rmr.gather <- gather(spread.data,key="time",value="rmr",
                           which(colnames(spread.data) %in% colnames(spread.data)[which(grepl("[0-9]",colnames(spread.data)) & grepl("rmr",colnames(spread.data)))]))     
      rmr.gather$time <- gsub("_rmr","",rmr.gather$time)
      rmr.data <- merge(raw.gather,rmr.gather) %>% select(genotype,dose,plate,insecticide,date,position,time,raw,rmr)
      if(write==T){
            write.csv(x=rmr.data,file=file.name,row.names = F)
      }
      return(rmr.data)
}



conf.diff <- function(data,level=.95){
      ci <- qnorm(1-((1-level)/2))*sd(data)/sqrt(length(data))
      return(ci)
}

se <- function(x){
      return(sd(x)/sqrt(length(x)))
}


dgrp.gwas.parse <- function(raw.gwas){
      split.gwas <- lapply(raw.gwas,strsplit," ")
      parsed.gwas <- vector("list",length(split.gwas))
      for(i in 1:length(split.gwas)){
            parsed.gwas[[i]] <- split.gwas[[i]][[1]][1:10] 
            parsed.gwas[[i]][11] <- unlist(strsplit(split.gwas[[i]][[1]][11],"\\|"))[2]
      }
      
      clean.gwas <- data.frame(matrix(ncol=11,nrow=length(parsed.gwas)-1))
      
      for(i in 1:nrow(clean.gwas)){
            clean.gwas[i,] <- parsed.gwas[[i]]
      }
      
      
      colnames(clean.gwas) <- c(clean.gwas[1,][1:10],"Gene.Name")
      clean.gwas <- clean.gwas[-1,]
      clean.gwas$SinglePval <- as.numeric(clean.gwas$SinglePval)
      parsed.gwas <- clean.gwas$SinglePval
      clean.gwas <- clean.gwas[which(parsed.gwas<1e-5),]
      for (i in 1:length(clean.gwas$Gene.Name)) if(clean.gwas$Gene.Name[i]=="") clean.gwas$Gene.Name[i] <- "intergenic"
      return(data.frame(clean.gwas))
}


mod.encode.enrich <- function(raw.data.file="/home/sdenecke/Documents/Large_Files/External_Datasets/Modencode.tsv",
                              gene.list=readLines("/home/sdenecke/Dropbox/PhD_Research/Research_Projects_Ongoing/DGRP_Paper/Results/GWAS/Unique_GWAS_Candidate_List.txt"),
                              parent.library="modENCODE_mRNA-Seq_tissues",
                              stage.search.string="L3",
                              tissue.search.string="CNS"){ 
      Upreg.subset <- raw.data.file %>% 
            filter(Parent_library_name==parent.library) %>% 
            filter(GeneSymbol %in% gene.list) %>%
            select(GeneSymbol,RNASource_name,RPKM_value) %>%
            filter(grepl(stage.search.string,RNASource_name)) %>%
            mutate(in.tissue.list=grepl(tissue.search.string,RNASource_name)) %>%
            group_by(in.tissue.list,GeneSymbol) %>% summarize(avg=mean(RPKM_value)) %>%
            spread(in.tissue.list,value=avg)
      prop.enriched.subset <- length(which(Upreg.subset$`TRUE`>2*Upreg.subset$`FALSE`))/length(Upreg.subset$`TRUE`)
      enriched.subset <- Upreg.subset %>% mutate(ratio=`TRUE`/`FALSE`) %>%
            filter(ratio>=2) %>%
            arrange(desc(ratio))
      colnames(enriched.subset) <- c("GeneSymbol","Other.Tissues","Your.Tissue","Ratio")
      
      Upreg.genome <- raw.data.file %>% 
            filter(Parent_library_name==parent.library) %>% 
            select(GeneSymbol,RNASource_name,RPKM_value) %>%
            filter(grepl(stage.search.string,RNASource_name)) %>%
            mutate(in.tissue.list=grepl(tissue.search.string,RNASource_name)) %>%
            group_by(in.tissue.list,GeneSymbol) %>% summarize(avg=mean(RPKM_value)) %>%
            spread(in.tissue.list,value=avg)  
      
      prop.enriched.genome <- length(which(Upreg.genome$`TRUE`>2*Upreg.genome$`FALSE`))/length(Upreg.genome$`TRUE`)
      enriched.genome <- Upreg.genome %>% mutate(ratio=`TRUE`/`FALSE`) %>% 
            filter(ratio>=2) %>%
            arrange(desc(ratio)) %>%
            mutate(in.gene.list=GeneSymbol %in% gene.list)
      
      
      write.csv(enriched.subset,file=paste(stage.search.string,tissue.search.string,"Gene_Enrichment.csv",sep="_"))
      comparison <- c(prop.enriched.subset,prop.enriched.genome) 
      names(comparison) <- c(paste(stage.search.string,tissue.search.string),"genome.wide") 
      
      
      enriched.genome2 <- data.frame(enriched.genome)
      for (i in enriched.genome2[,"GeneSymbol"]){
            if(enriched.genome2[which(enriched.genome2$GeneSymbol==i),"FALSE."]==0){ 
                  enriched.genome2 <- enriched.genome2[-which(enriched.genome2$GeneSymbol==i),]
            }
      }
      
      significant.enrichment <- with(enriched.genome2,
                                     summary(lm(ratio~in.gene.list)))$coefficients[2,4]
      names(significant.enrichment) <- "p.value"
      return(c(comparison,significant.enrichment))
} 


clean.all.associations <- function(all.assoc){
      ## Clean up your data frame containing all your associations 
      all.associations <- all.assoc
      all.associations <- select(all.associations,ID,SingleMixedPval)
      variant <- unlist(strsplit(all.associations$ID,split="_"))
      chromosome <- variant[seq(1, length(variant), 3)]
      base <- as.numeric(variant[seq(2, length(variant), 3)])
      chromosome<-gsub('X',1,chromosome)
      chromosome<-gsub('2L',2,chromosome)
      chromosome<-gsub('2R',3,chromosome)
      chromosome<-gsub('3L',4,chromosome)
      chromosome<-gsub('3R',5,chromosome)
      chromosome<-as.numeric(chromosome)
      gwas.plot.data <- all.associations %>% mutate(CHR=chromosome,BP=base) %>% select(SingleMixedPval,ID,CHR,BP)
      colnames(gwas.plot.data)=c('P','SNP','CHR','BP')
     
      return(gwas.plot.data)
}         




equation.extract <- function(var1,var2,data,subset.var="",subset.condition="",sub=T){
      if(sub==T){
            model <- lm(get(var1)~get(var2),data=data,subset=which(data[[subset.var]]==subset.condition))
            sum <- summary(model)$coefficients
            equation <- paste("y=",round(sum[2,1],2),"x + ",round(sum[1,1],2),"\n p value=",round(sum[2,4],6),sep="")
      }else{
            model <- lm(get(var1)~get(var2),data=data)
            sum <- summary(model)$coefficients
            equation <- paste("y=",round(sum[2,1],2),"x + ",round(sum[1,1],2),"\n p value=",round(sum[2,4],6),sep="")  
      
      }
      return(equation)
}


equation.extract.hplc <- function(var1,var2,data,subset.var="",subset.condition="",sub=T,bonf=1){
      if(sub==T){
            model <- lm(get(var1)~get(var2),data=data,subset=which(data[[subset.var]]==subset.condition))
            sum <- summary(model)$coefficients
            equation <- paste("y=",round(sum[2,1],2),"x + ",round(sum[1,1],2),"\n p value=",round(sum[2,4],6),sep="")
      }else{
            model <- lm(get(var1)~get(var2),data=data)
            sum <- summary(model)$coefficients
            p.val <- sum[2,4]*bonf
            if(p.val > .05){
                  equation <- ""#paste("p value = ",min(round(p.val,2),1),sep="")  
            }else if (p.val <.05 & p.val > .001){
                  equation <- "*"
            }else if (p.val <.001 & p.val > .0001){
                  equation <- "**"
            }else if (p.val<.0001){
                  equation <- "***"
            }
      }
      return(equation)
}




rt.data.parse <- function(data){
      colnames(data) <- gsub(" ",".",colnames(data))
      data.na.rm <- data[complete.cases(select(data,-Biological.Set.Name)),]
      clean.data <- data.na.rm %>%
            select(Target:Cq, -Biological.Set.Name)  %>%
            filter(Content!="NTC")
      l <- vector("list",length(unique(clean.data$Content)))
      names(l) <- unique(clean.data$Content)
      for (i in unique(clean.data$Content)){
            tech.reps <- clean.data %>% filter(Content==i)
            if(abs(range(tech.reps$Cq)[1]-range(tech.reps$Cq)[2])>1 | length(tech.reps$Cq)<3){
                  dodgy.reps <- tech.reps$Cq
                  names(dodgy.reps) <- as.character(1:length(dodgy.reps))
                  ord.dod.rep <- dodgy.reps[order(dodgy.reps)]
                  if(ord.dod.rep[2]-ord.dod.rep[1]>.5 & ord.dod.rep[3]-ord.dod.rep[2]>.5 | length(ord.dod.rep)<3){
                        l[[i]] <- NULL
                  }else if(ord.dod.rep[2]-ord.dod.rep[1] > ord.dod.rep[3]-ord.dod.rep[2]){
                        tech.reps <- tech.reps[-as.numeric(names(ord.dod.rep[1])),]
                        l[[i]] <- tech.reps
                  }else if(ord.dod.rep[2]-ord.dod.rep[1] < ord.dod.rep[3]-ord.dod.rep[2]){
                        tech.reps <- tech.reps[-as.numeric(names(ord.dod.rep[3])),]
                        l[[i]] <- tech.reps
                  }
            }else{
                  l[[i]] <- tech.reps
            }
      }
      filtered.data <- rbindlist(l,use.names=T)
      
      
      ## Vizualize the distribution of your data by sample gene and biolgoical replicate
      raw.number <- as.numeric(unlist(strsplit(filtered.data$Content,split="-"))[seq(from=2,to=nrow(filtered.data)*2,by=2)])
      
      biol.rep <- function(x){
            a <- x%/%5
            b <- x-(5*a)
            if(b==0){
                  b <- b+5
            }
            return(b)
      }
      filtered.data <- mutate(filtered.data,biological.rep=as.character(sapply(raw.number,biol.rep)))
      return(filtered.data)
      
}  


allele.ad <- function(x){
      y <- c()
      for (i in 1:length(x)){
            if(x[i]=="RAL_138" | x[i]=="RAL_360"){ 
                  y[i] <- "AA"
            }else{ 
                  y[i] <- "M"
            }
      }
      return(y)
}


multi.gsub <- function(vector,find,replace){
      for (i in 1:length(find)){
            vector <- gsub(find[i],replace[i],vector)
      }
      return(vector)
} 



transpose.collapse <- function(data,column){
      cols <- unique(data[[column]])
      multiplier <- length(data[[column]])/length(cols)
      rows <- colnames(data)[which(colnames(data)!=column)]
      emat <- matrix(ncol=length(cols),nrow=length(rows)*multiplier,dimnames=list(rep(rows,each=multiplier),cols))
      
      for(i in cols){
            for(j in rows){
                  df <- data[which(data[[column]]==i),j]
                  emat[which(rownames(emat)==j),i] <- df
            }
      }
      gen <- rownames(emat)
      rownames(emat) <- NULL
      final.data <- data.frame(emat,stringsAsFactors = F)
      final.data <- mutate(final.data,genotype=gen)
      return(final.data)
}




transcript.combine <- function(male,female){
      a <- data.frame(rbindlist(list(male,female),use.names=T,fill=T))
      names <- unlist(strsplit(colnames(a),split="\\."))[seq(2,(length(colnames(a))*2)-1,by=2)]
      v <- c()
      v2 <- c()
      for (i in unique(names)){
            v <- na.omit(c(v,grep(i,colnames(a))[1]))
            v2 <- na.omit(c(v2,grep(i,colnames(a))[2]))
      }
      
      a1 <- a[,c(1,v)]
      a2 <- a[,c(1,v2)]
      
      colnames(a1) <- c("gene",unique(names))
      a2.names <- unlist(strsplit(colnames(a2),split="\\."))[seq(2,(length(colnames(a2))*2)-1,by=2)]
      colnames(a2) <- c("gene",unique(a2.names))
      
      final <- data.frame(rbindlist(list(a1,a2),use.names=T,fill=T))
      colnames(final) <- multi.gsub(colnames(final),c("line"),c("RAL"))
      return(data.table(final))
}
