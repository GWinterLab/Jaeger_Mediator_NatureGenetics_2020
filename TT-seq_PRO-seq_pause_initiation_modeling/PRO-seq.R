

### create PRO-seq coverage track ###


if (TRUE){
  PROseqBED_hg38_plus = import(con=file.path(prewd,"<plus strand>.bw"))
  strand(PROseqBED_hg38_plus) = "+"
  names(PROseqBED_hg38_plus) = NULL
  
  PROseqBED_hg38_minus = import(con=file.path(prewd,"<minus strand>.bw"))
  strand(PROseqBED_hg38_minus) = "-"
  names(PROseqBED_hg38_minus) = NULL
  
  PROseq.anno = rbind(as.data.frame(PROseqBED_hg38_plus,stringsAsFactors = FALSE),as.data.frame(PROseqBED_hg38_minus,stringsAsFactors = FALSE))
  head(PROseq.anno)
  dim(PROseq.anno)
  
  PROseq.anno = cbind(PROseq.anno,"length" = abs(PROseq.anno[,"end"] - PROseq.anno[,"start"] + 1))
  unique(PROseq.anno[,"length"])
  
  colnames(PROseq.anno) = c("chr",setdiff(colnames(PROseq.anno),"seqnames"))
  
  PROseq.anno[,"score"] = abs(PROseq.anno[,"score"])
  
  dm.PRO.seq.DMSO1h = sapply(c("dm6_chr2L","dm6_chr2R","dm6_chr3L","dm6_chr3R","dm6_chr4"),function(x){sum(PROseq.anno[which(PROseq.anno[,"chr"] %in% x),"score"])})
  
  save(dm.PRO.seq.DMSO1h,file=file.path("TTseqObjects","JaegerWinter2019TTseqKBM7","dm.PRO.seq.DMSO1h.RData"))
  
  jaeger.winter.2019.PROseq.KBM7.DMSO1h.track.list = list()
  jaeger.winter.2019.PROseq.KBM7.DMSO1h.track.list[["+"]] = list()
  jaeger.winter.2019.PROseq.KBM7.DMSO1h.track.list[["-"]] = list()
  
  registerDoParallel(cores = mc.cores)
  build.list = function(chr,which.strand){
    pro.seq.coverage.vec = Rle(0,human.chrs.lengths[chr])
    PROseq.anno.part = PROseq.anno[which(PROseq.anno[,"strand"] == which.strand & PROseq.anno[,"chr"] == chr),]
    if (dim(PROseq.anno.part)[1] > 0){
      for (j in 1:max(unique(PROseq.anno[,"length"]))){
        PROseq.anno.part.etj = PROseq.anno.part[which(PROseq.anno.part[,"length"] >= j),]
        if (dim(PROseq.anno.part.etj)[1] > 0){pro.seq.coverage.vec[PROseq.anno.part.etj[,"start"]+(j-1)] = PROseq.anno.part.etj[,"score"]}
      }
    }		
    return(pro.seq.coverage.vec)
  }
  
  jaeger.winter.2019.PROseq.KBM7.DMSO1h.track.list[["+"]] = foreach(n = human.chrs) %dopar% build.list(n,"+")
  names(jaeger.winter.2019.PROseq.KBM7.DMSO1h.track.list[["+"]]) = human.chrs
  
  jaeger.winter.2019.PROseq.KBM7.DMSO1h.track.list[["-"]] = foreach(n = human.chrs) %dopar% build.list(n,"-")
  names(jaeger.winter.2019.PROseq.KBM7.DMSO1h.track.list[["-"]]) = human.chrs
  
  for (which.chr in human.chrs){
    jaeger.winter.2019.PROseq.KBM7.DMSO1h.track.list.chr = list()
    
    jaeger.winter.2019.PROseq.KBM7.DMSO1h.track.list.chr[["+"]] = jaeger.winter.2019.PROseq.KBM7.DMSO1h.track.list[["+"]][[which.chr]]
    jaeger.winter.2019.PROseq.KBM7.DMSO1h.track.list.chr[["-"]] = jaeger.winter.2019.PROseq.KBM7.DMSO1h.track.list[["-"]][[which.chr]]
    
    save(jaeger.winter.2019.PROseq.KBM7.DMSO1h.track.list.chr,file = file.path("PROseqRleTracks",paste0("jaeger.winter.2019.PROseq.KBM7.DMSO1h.track.list.",which.chr,".RData")))
  }
  
  rm(list = setdiff(ls(),basic.objects))
  gc()
}



