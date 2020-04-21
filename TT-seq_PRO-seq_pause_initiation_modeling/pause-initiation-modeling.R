

### libraries ###


library(Rsamtools)
library(foreach)
library(doParallel)
library(GenomicRanges)
library(GenomicAlignments)
library(rtracklayer)
library(DESeq2)


### number of cores ###


mc.cores = 12 #detectCores()


### basic objects ###


load(file.path("TTseqObjects","JaegerWinter2019TTseqKBM7","bam.files.RData"))
load(file.path("TTseqObjects","JaegerWinter2019TTseqKBM7","bam.files.mat.RData"))

basic.objects = c(ls(),"basic.objects","bam.file","bam.files","bam.files.mat")


### initiation rate estimation ###


if (TRUE){
  load(file.path("AnnotationObjects","human.refseq.anno.transcript.merge.RData"))
  load(file.path("AnnotationObjects","human.refseq.constitutive.exons.RData"))
  load(file.path("TTseqObjects","JaegerWinter2019TTseqKBM7","size.factors.RData"))
  load(file.path("TTseqObjects","JaegerWinter2019TTseqKBM7","conversion.factor.to.amount.per.cell.RData"))
  load(file.path("TTseqObjects","JaegerWinter2019TTseqKBM7","expressed.TR.2.20.RData"))
  
  dir.create(file.path("..","Output","Visualization","JaegerWinter2019TTseqKBM7"))
  dir.create(file.path("..","Output","Visualization","JaegerWinter2019TTseqKBM7","InitiationElongationPlots"))
  
  # subselection
  
  TR.subsel = intersect(expressed.TR.2.20,unique(human.refseq.constitutive.exons[which(human.refseq.constitutive.exons[,"first"] == TRUE),"TR_id"]))
  length(TR.subsel)
  
  dir.create(file.path("TTseqObjects","JaegerWinter2019TTseqKBM7"))
  
  save(TR.subsel,file=file.path("TTseqObjects","JaegerWinter2019TTseqKBM7","TR.subsel.RData"))
    
  ### up to 10 kb ###
  
  load(file.path("AnnotationObjects","human.refseq.anno.transcript.500.to.10000.merge.RData"))
  load(file.path("TTseqObjects","JaegerWinter2019TTseqKBM7","human.refseq.anno.transcript.500.to.10000.coverage.antisense.corrected.RData"))
  
  # transcript coverage and length
  
  human.refseq.anno.transcript.500.to.10000.coverage.antisense.corrected.TTseq = apply(t(t(human.refseq.anno.transcript.500.to.10000.coverage.antisense.corrected[TR.subsel,c("L_with_1_AACCAG_1_KBM7_DMSO1h.bam","L_with_2_ACCTCA_5_KBM7_DMSO1h.bam"),drop = FALSE])/size.factors[c("L_with_1_AACCAG_1_KBM7_DMSO1h.tabular","L_with_2_ACCTCA_5_KBM7_DMSO1h.tabular")]),1,sum)
  names(human.refseq.anno.transcript.500.to.10000.coverage.antisense.corrected.TTseq) = TR.subsel
  
  human.refseq.anno.transcript.500.to.10000.lengths = human.refseq.anno.transcript.500.to.10000.merge[TR.subsel,"length"]
  
  # initiation rate calculation
  
  initiation.rate.transcript = human.refseq.anno.transcript.500.to.10000.coverage.antisense.corrected.TTseq/(5*human.refseq.anno.transcript.500.to.10000.lengths)
  initiation.rate.transcript = initiation.rate.transcript[TR.subsel]
  save(initiation.rate.transcript,file=file.path("TTseqObjects","JaegerWinter2019TTseqKBM7","initiation.rate.transcript.RData"))
    
  initiation.rate.per.cell = initiation.rate.transcript/conversion.factor.to.amount.per.cell # conversion.factor.to.amount.per.cell is used to calibrate the initiation frequency to a literature estimate (size factor 𝜅)
  
  save(initiation.rate.per.cell,file=file.path("TTseqObjects","JaegerWinter2019TTseqKBM7","initiation.rate.per.cell.RData"))
  
  rm(list = setdiff(ls(),basic.objects))
  gc()
}


### pause site estimation ###


if (TRUE){
  load(file.path("AnnotationObjects","human.refseq.anno.transcript.merge.RData"))
  load(file.path("AnnotationObjects","human.refseq.constitutive.exons.RData"))
  load(file.path("TTseqObjects","JaegerWinter2019TTseqKBM7","expressed.TR.2.20.RData"))
  load(file.path("TTseqObjects","JaegerWinter2019TTseqKBM7","TR.subsel.RData"))
  
  dir.create(file.path("..","Output","Visualization","JaegerWinter2019TTseqKBM7"))
  dir.create(file.path("..","Output","Visualization","JaegerWinter2019TTseqKBM7","PausingPlots"))
  
  human.refseq.constitutive.exon.lengths.subsel = human.refseq.constitutive.exons[which(human.refseq.constitutive.exons[,"first"] == TRUE),]
  human.refseq.constitutive.exon.lengths = human.refseq.constitutive.exon.lengths.subsel[,"length"]
  names(human.refseq.constitutive.exon.lengths) = human.refseq.constitutive.exon.lengths.subsel[,"TR_id"]
  
  gene.anno = human.refseq.anno.transcript.merge[TR.subsel,]
  dim(gene.anno)
  
  human.refseq.anno.transcript.pausing.list = list()
  
  for (bam.file in <PRO-seq track>){
    index.subsets = split(1:nrow(gene.anno),paste(gene.anno[,"strand"],"_",gene.anno[,"chr"],sep = ""))
    
    pausing.list = list()
    
    build.pausing.list = function(j){
      from.transcript = strand.chr.gene.anno[j,"start"]
      to.transcript = strand.chr.gene.anno[j,"end"]
      strand.transcript = as.character(strand.chr.gene.anno[j,"strand"])
      
      if(strand.transcript == "+"){
        strand.length = length(strand.chr.pausing.from.bam)
        positions = from.transcript:to.transcript
        valid.positions = positions[positions > 0 & positions <= strand.length]				
        vec.template = Rle(rep(0,to.transcript - from.transcript + 1))
        if (length(valid.positions) > 0){
          cov.vec = strand.chr.pausing.from.bam[valid.positions]
          vec.template[positions > 0 & positions <= strand.length] = cov.vec
        }
      } else {
        strand.length = length(strand.chr.pausing.from.bam)
        positions = to.transcript:from.transcript
        valid.positions = positions[positions > 0 & positions <= strand.length]
        vec.template = Rle(rep(0,to.transcript - from.transcript + 1))
        if (length(valid.positions) > 0){
          cov.vec = strand.chr.pausing.from.bam[valid.positions]
          vec.template[positions > 0 & positions <= strand.length] = cov.vec
        }
        return(vec.template)
      }
    }
    
    for (index.subset in names(index.subsets)){
      chr.pausing.list = get(load(file.path("PROseqRleTracks",paste0(bam.file,".track.list.",unlist(strsplit(index.subset,split = "_"))[2],".RData"))))
      
      strand.chr.pausing.from.bam = chr.pausing.list[[unlist(strsplit(index.subset,split = "_"))[1]]]
      strand.chr.gene.anno = gene.anno[index.subsets[[index.subset]],c("start","end","strand","length")]
      
      registerDoParallel(cores = mc.cores)
      pausing.list[[index.subset]] = foreach(n = 1:nrow(strand.chr.gene.anno),.noexport = setdiff(ls(),c("strand.chr.pausing.from.bam","strand.chr.gene.anno","build.pausing.list"))) %dopar% build.pausing.list(n)
    }
    
    pausing.list = unlist(pausing.list,recursive = FALSE,use.names = FALSE)
    names(pausing.list) = as.character(gene.anno[as.numeric(unlist(index.subsets)),"id"])
    human.refseq.anno.transcript.pausing.list[[bam.file]] = pausing.list
  }	
  
  human.refseq.anno.transcript.pausing.list = Reduce('+',lapply(human.refseq.anno.transcript.pausing.list,RleList))
  
  find.polymerase.pausing.sites = function(j){
    rle.vec = human.refseq.anno.transcript.pausing.list[[j]][1:(human.refseq.constitutive.exon.lengths[j] - 5)]
    try({rle.vec[rle.vec <= quantile(rle.vec[rle.vec != 0],0.5,na.rm = TRUE)*5] = 0},silent = TRUE)
    return(c(which.max(rle.vec),runValue(rle.vec[which.max(rle.vec)])))
  }
  
  polymerase.pausing.sites = foreach(n = names(human.refseq.anno.transcript.pausing.list),.combine = "rbind",.noexport = setdiff(ls(),c("human.refseq.anno.transcript.pausing.list","find.polymerase.pausing.sites"))) %dopar% find.polymerase.pausing.sites(n)
  rownames(polymerase.pausing.sites) = names(human.refseq.anno.transcript.pausing.list)
  polymerase.pausing.sites = polymerase.pausing.sites[polymerase.pausing.sites[,1] > 1,1]
  summary(polymerase.pausing.sites)
  length(polymerase.pausing.sites)
  
  ### human.refseq.anno.transcript.pausing.annotation ###
  
  human.refseq.anno.transcript.pausing.annotation = human.refseq.anno.transcript.merge[names(polymerase.pausing.sites),c("chr","strand","start","end","type","id")]
  human.refseq.anno.transcript.pausing.annotation = cbind(human.refseq.anno.transcript.pausing.annotation,"pausing.position" = polymerase.pausing.sites[rownames(human.refseq.anno.transcript.pausing.annotation)],"pausing.distance" = polymerase.pausing.sites[rownames(human.refseq.anno.transcript.pausing.annotation)])
  human.refseq.anno.transcript.pausing.annotation[which(human.refseq.anno.transcript.pausing.annotation[,"strand"] == "+"),"pausing.position"] = human.refseq.anno.transcript.pausing.annotation[which(human.refseq.anno.transcript.pausing.annotation[,"strand"] == "+"),"start"] + human.refseq.anno.transcript.pausing.annotation[which(human.refseq.anno.transcript.pausing.annotation[,"strand"] == "+"),"pausing.position"]
  human.refseq.anno.transcript.pausing.annotation[which(human.refseq.anno.transcript.pausing.annotation[,"strand"] == "-"),"pausing.position"] = human.refseq.anno.transcript.pausing.annotation[which(human.refseq.anno.transcript.pausing.annotation[,"strand"] == "-"),"end"] - human.refseq.anno.transcript.pausing.annotation[which(human.refseq.anno.transcript.pausing.annotation[,"strand"] == "-"),"pausing.position"]
  human.refseq.anno.transcript.pausing.annotation = human.refseq.anno.transcript.pausing.annotation[,c("chr","strand","start","pausing.position","pausing.distance","end","type","id")]
  human.refseq.anno.transcript.pausing.annotation[,c("type")] = "pausing_position"
  human.refseq.anno.transcript.pausing.annotation = cbind(human.refseq.anno.transcript.pausing.annotation,"length" = abs(human.refseq.anno.transcript.pausing.annotation[,"start"] - human.refseq.anno.transcript.pausing.annotation[,"end"]) + 1)  
  human.refseq.anno.transcript.pausing.annotation = cbind(human.refseq.anno.transcript.pausing.annotation,"pausing.region.length" = NA)
  human.refseq.anno.transcript.pausing.annotation[which(human.refseq.anno.transcript.pausing.annotation[,"strand"] == "+"),"pausing.region.length"] = abs(human.refseq.anno.transcript.pausing.annotation[which(human.refseq.anno.transcript.pausing.annotation[,"strand"] == "+"),"start"] - human.refseq.anno.transcript.pausing.annotation[which(human.refseq.anno.transcript.pausing.annotation[,"strand"] == "+"),"pausing.position"]) + 1  
  human.refseq.anno.transcript.pausing.annotation[which(human.refseq.anno.transcript.pausing.annotation[,"strand"] == "-"),"pausing.region.length"] = abs(human.refseq.anno.transcript.pausing.annotation[which(human.refseq.anno.transcript.pausing.annotation[,"strand"] == "-"),"pausing.position"] - human.refseq.anno.transcript.pausing.annotation[which(human.refseq.anno.transcript.pausing.annotation[,"strand"] == "-"),"end"]) + 1
  head(human.refseq.anno.transcript.pausing.annotation)
  
  save(human.refseq.anno.transcript.pausing.annotation,file=file.path("TTseqObjects","JaegerWinter2019TTseqKBM7","human.refseq.anno.transcript.pausing.annotation.RData"))
  
  human.refseq.anno.transcript.pausing.annotation.ranges = as(human.refseq.anno.transcript.pausing.annotation,"GRanges")
  save(human.refseq.anno.transcript.pausing.annotation.ranges,file=file.path("TTseqObjects","JaegerWinter2019TTseqKBM7","human.refseq.anno.transcript.pausing.annotation.ranges.RData"))
  
  rm(list = setdiff(ls(),basic.objects))
  gc()
}


### pause duration estimation ###


if (TRUE){
  load(file.path("TTseqObjects","JaegerWinter2019TTseqKBM7","human.refseq.anno.transcript.pausing.annotation.RData"))
  
  human.refseq.anno.transcript.pausing.annotation[which(human.refseq.anno.transcript.pausing.annotation[,"strand"] == "+"),c("start","end")] = cbind(human.refseq.anno.transcript.pausing.annotation[which(human.refseq.anno.transcript.pausing.annotation[,"strand"] == "+"),c("pausing.position")] - 100 + abs(pmin(human.refseq.anno.transcript.pausing.annotation[which(human.refseq.anno.transcript.pausing.annotation[,"strand"] == "+"),c("pausing.distance")] - 100,0)),human.refseq.anno.transcript.pausing.annotation[which(human.refseq.anno.transcript.pausing.annotation[,"strand"] == "+"),c("pausing.position")] + 100 + abs(pmin(human.refseq.anno.transcript.pausing.annotation[which(human.refseq.anno.transcript.pausing.annotation[,"strand"] == "+"),c("pausing.distance")] - 100,0)))
  human.refseq.anno.transcript.pausing.annotation[which(human.refseq.anno.transcript.pausing.annotation[,"strand"] == "-"),c("start","end")] = cbind(human.refseq.anno.transcript.pausing.annotation[which(human.refseq.anno.transcript.pausing.annotation[,"strand"] == "-"),c("pausing.position")] - 100 - abs(pmin(human.refseq.anno.transcript.pausing.annotation[which(human.refseq.anno.transcript.pausing.annotation[,"strand"] == "-"),c("pausing.distance")] - 100,0)),human.refseq.anno.transcript.pausing.annotation[which(human.refseq.anno.transcript.pausing.annotation[,"strand"] == "-"),c("pausing.position")] + 100 - abs(pmin(human.refseq.anno.transcript.pausing.annotation[which(human.refseq.anno.transcript.pausing.annotation[,"strand"] == "-"),c("pausing.distance")] - 100,0)))
  
  ### human.refseq.anno.transcript.mNETseq.pausing.coverage ###
  
  human.refseq.anno.transcript.mNETseq.pausing.coverage = list()
  for (bam.file in <PRO-seq track>){
    index.subsets = split(1:nrow(human.refseq.anno.transcript.pausing.annotation),paste(human.refseq.anno.transcript.pausing.annotation[,"strand"],"_",human.refseq.anno.transcript.pausing.annotation[,"chr"],sep = ""))
    
    coverage.list = list()
    
    build.coverage.list = function(j){
      from.transcript = strand.chr.human.refseq.anno.transcript.pausing.annotation[j,"start"]
      to.transcript = strand.chr.human.refseq.anno.transcript.pausing.annotation[j,"end"]
      return(sum(as.vector(strand.chr.coverage.from.bam[from.transcript:to.transcript])))
    }
    
    for (index.subset in names(index.subsets)){
      chr.coverage.list = get(load(file.path("PROseqRleTracks",paste0(bam.file,".track.list.",unlist(strsplit(index.subset,split = "_"))[2],".RData"))))
      
      strand.chr.coverage.from.bam = chr.coverage.list[[unlist(strsplit(index.subset,split = "_"))[1]]]
      strand.chr.human.refseq.anno.transcript.pausing.annotation = human.refseq.anno.transcript.pausing.annotation[index.subsets[[index.subset]],c("start","end")]
      
      registerDoParallel(cores = mc.cores)
      coverage.list[[index.subset]] = foreach(n = 1:nrow(strand.chr.human.refseq.anno.transcript.pausing.annotation),.noexport = setdiff(ls(),c("strand.chr.coverage.from.bam","strand.chr.human.refseq.anno.transcript.pausing.annotation","build.coverage.list"))) %dopar% build.coverage.list(n)
    }
    human.refseq.anno.transcript.mNETseq.pausing.coverage[[bam.file]] = unlist(coverage.list,recursive = TRUE,use.names = FALSE)
  }	
  human.refseq.anno.transcript.mNETseq.pausing.coverage = sapply(human.refseq.anno.transcript.mNETseq.pausing.coverage,c)
  rownames(human.refseq.anno.transcript.mNETseq.pausing.coverage) = as.character(human.refseq.anno.transcript.pausing.annotation[as.numeric(unlist(index.subsets,recursive = TRUE,use.names = FALSE)),"id"])
  colnames(human.refseq.anno.transcript.mNETseq.pausing.coverage) = c("jaeger.winter.2019.PROseq.KBM7.DMSO1h")
  
  ### human.refseq.anno.transcript.mNETseq.pausing.antisense.coverage ###
  
  human.refseq.anno.transcript.mNETseq.pausing.antisense.coverage = list()
  for (bam.file in <PRO-seq track>){
    index.subsets = split(1:nrow(human.refseq.anno.transcript.pausing.annotation),paste(Vectorize(strand.switch)(as.character(human.refseq.anno.transcript.pausing.annotation[,"strand"])),"_",human.refseq.anno.transcript.pausing.annotation[,"chr"],sep = ""))
    
    coverage.list = list()
    
    build.coverage.list = function(j){
      from.transcript = strand.chr.human.refseq.anno.transcript.pausing.annotation[j,"start"]
      to.transcript = strand.chr.human.refseq.anno.transcript.pausing.annotation[j,"end"]
      return(sum(as.vector(strand.chr.coverage.from.bam[from.transcript:to.transcript])))
    }
    
    for (index.subset in names(index.subsets)){
      chr.coverage.list = get(load(file.path("PROseqRleTracks",paste0(bam.file,".track.list.",unlist(strsplit(index.subset,split = "_"))[2],".RData"))))
      
      strand.chr.coverage.from.bam = chr.coverage.list[[unlist(strsplit(index.subset,split = "_"))[1]]]
      strand.chr.human.refseq.anno.transcript.pausing.annotation = human.refseq.anno.transcript.pausing.annotation[index.subsets[[index.subset]],c("start","end")]
      
      registerDoParallel(cores = mc.cores)
      coverage.list[[index.subset]] = foreach(n = 1:nrow(strand.chr.human.refseq.anno.transcript.pausing.annotation),.noexport = setdiff(ls(),c("strand.chr.coverage.from.bam","strand.chr.human.refseq.anno.transcript.pausing.annotation","build.coverage.list"))) %dopar% build.coverage.list(n)
    }
    human.refseq.anno.transcript.mNETseq.pausing.antisense.coverage[[bam.file]] = unlist(coverage.list,recursive = TRUE,use.names = FALSE)
  }	
  human.refseq.anno.transcript.mNETseq.pausing.antisense.coverage = sapply(human.refseq.anno.transcript.mNETseq.pausing.antisense.coverage,c)
  rownames(human.refseq.anno.transcript.mNETseq.pausing.antisense.coverage) = as.character(human.refseq.anno.transcript.pausing.annotation[as.numeric(unlist(index.subsets,recursive = TRUE,use.names = FALSE)),"id"])
  colnames(human.refseq.anno.transcript.mNETseq.pausing.antisense.coverage) = <PRO-seq track>
  
  ### correct human.refseq.anno.transcript.mNETseq.pausing.coverage for antisense bias ###
  
  load(file.path("PROseqObjects","JaegerWinter2019TTseqKBM7","antisense.bias.ratio.RData"))
  
  human.refseq.anno.transcript.mNETseq.pausing.coverage.antisense.corrected = t(t(human.refseq.anno.transcript.mNETseq.pausing.coverage - t(t(human.refseq.anno.transcript.mNETseq.pausing.antisense.coverage)*antisense.bias.ratio[colnames(human.refseq.anno.transcript.mNETseq.pausing.coverage)]))/(1 - antisense.bias.ratio[colnames(human.refseq.anno.transcript.mNETseq.pausing.coverage)]^2))
  human.refseq.anno.transcript.mNETseq.pausing.coverage.antisense.corrected[human.refseq.anno.transcript.mNETseq.pausing.coverage.antisense.corrected < 0] = 0
  
  load(file.path("TTseqObjects","JaegerWinter2019TTseqKBM7","initiation.rate.transcript.RData"))
  mNETseq.size.factor = ... # mNETseq.size.factor is used to calibrate the pause duration to a literature estimate (size factor σ)
  save(mNETseq.size.factor,file=file.path("TTseqObjects","JaegerWinter2019TTseqKBM7","mNETseq.size.factor.RData"))
  load(file.path("TTseqObjects","JaegerWinter2019TTseqKBM7","conversion.factor.to.amount.per.cell.RData"))
  
  common = rownames(human.refseq.anno.transcript.pausing.annotation)
  
  pause.duration = apply(human.refseq.anno.transcript.mNETseq.pausing.coverage.antisense.corrected[common,,drop = FALSE],1,sum)/(initiation.rate.transcript[common]*mNETseq.size.factor/conversion.factor.to.amount.per.cell)
  summary(pause.duration)
  
  save(pause.duration,file=file.path("TTseqObjects","JaegerWinter2019TTseqKBM7","pause.duration.RData"))
  
  # calibrate magnitude to literature estimate (with mNETseq.size.factor)
  
  if (FALSE){
    pause.duration.raji = get(load(file.path("TTseqObjects","GresselCramer2016TTseqRaji","1stModel","pause.duration.RData")))
    load(file.path("TTseqObjects","JaegerWinter2019TTseqKBM7","pause.duration.RData"))
    common = intersect(names(pause.duration.raji),names(pause.duration))
    median(pause.duration[common]/pause.duration.raji[common],na.rm = TRUE)
  }
  
  checkvector(pause.duration)
  checkvector(pause.duration[pause.duration > 0 & pause.duration < 20])
  
  plotsfkt = function(){
    hist(pause.duration[pause.duration > 0 & pause.duration < 20],main = "",xlab = expression(paste("Pause duration [min]")),col = "lightgrey",breaks = seq(0,20,length.out = 13))
    mtext(paste("mean:",round(mean(pause.duration[pause.duration > 0 & pause.duration < 20]),2)," median:",round(median(pause.duration[pause.duration > 0 & pause.duration < 20]),2)),3,1,col="darkgrey")
    mtext(expression(paste("Pause duration [min]")),3,2,cex=1.25)
  }
  plotit(filename = file.path("..","Output","Visualization","JaegerWinter2019TTseqKBM7","InitiationElongationPlots","pause.duration.distribution"),sw = 2,sh = 2,sres = 2,plotsfkt = plotsfkt,ww = 7,wh = 7,saveit = TRUE,notinR = TRUE,addformat = "pdf")
  
  load(file.path("TTseqObjects","JaegerWinter2019TTseqKBM7","initiation.rate.per.cell.RData"))
  
  common = intersect(names(initiation.rate.per.cell),names(pause.duration))
  
  cor(initiation.rate.per.cell[common],pause.duration[common],method = "spearman",use = "na.or.complete")
  
  checkvector(initiation.rate.per.cell[common])
  checkvector(pause.duration[common])
  
  plotsfkt = function(){
    heatscatter(initiation.rate.per.cell[common],pause.duration[common],main = "",xlab = expression(paste("Initiation rates [",cell^-1*min^-1,"]")),ylab = expression(paste("Pause duration [min]")),xlim = c(0,6),ylim = c(0,6))
    mtext(paste("Ehrensberger inequality: v/I >= 50 [bp]"),3,1,col="darkgrey")
    polygon(c(seq(200/(10*50),10,length.out = 100),10),c((200/50)/seq(200/(10*50),10,length.out = 100),10),col = convertcolor("lightgrey",50))
    mtext(paste("Ehrensberger zone"),3,2,cex=1.25)
  }
  plotit(filename = file.path("..","Output","Visualization","JaegerWinter2019TTseqKBM7","InitiationElongationPlots","ehrensberger.zone.pause.duration"),sw = 2,sh = 2,sres = 2,plotsfkt = plotsfkt,ww = 7,wh = 7,saveit = TRUE,notinR = TRUE,addformat = "pdf")
    
  rm(list = setdiff(ls(),basic.objects))
  gc()
}



