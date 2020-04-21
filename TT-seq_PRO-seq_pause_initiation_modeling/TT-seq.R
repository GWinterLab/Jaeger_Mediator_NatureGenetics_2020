

### libraries ###


library(Rsamtools)
library(foreach)
library(doParallel)
library(GenomicRanges)
library(GenomicAlignments)
library(rtracklayer)
library(DESeq2)


### strand.switch ###


strand.switch = function(which.strand){switch(which.strand,"+" = "-","-" = "+")}


### create bam.files and bam.files.mat ###


if (TRUE){
	bam.files = list.files(file.path(prewd,"RawData","JaegerWinter2019TTseqKBM7","BamFiles"))
	bam.files = bam.files[which(substr(bam.files,nchar(bam.files)-3,nchar(bam.files)) == ".bam")]
	
	dir.create(file.path("TTseqObjects"))
	dir.create(file.path("TTseqObjects","JaegerWinter2019TTseqKBM7"))
	
	save(bam.files,file=file.path("TTseqObjects","JaegerWinter2019TTseqKBM7","bam.files.RData"))
	
	bam.files.mat = matrix(NA,nrow = length(bam.files),ncol = 9)
	rownames(bam.files.mat) = substr(bam.files,1,nchar(bam.files)-4)
	colnames(bam.files.mat) = c("name","fraction","fragmentation","replicate","barcode","sample","cellline","condition")
	
	split.mat = t(sapply(strsplit(rownames(bam.files.mat),"_"),c))
	
	bam.files.mat[,"name"] = bam.files
	bam.files.mat[,c("fraction","fragmentation","replicate","barcode","sample","cellline","condition")] = split.mat
	head(bam.files.mat)

	save(bam.files.mat,file=file.path("TTseqObjects","JaegerWinter2019TTseqKBM7","bam.files.mat.RData"))
}


### basic objects ###


load(file.path("TTseqObjects","JaegerWinter2019TTseqKBM7","bam.files.RData"))
load(file.path("TTseqObjects","JaegerWinter2019TTseqKBM7","bam.files.mat.RData"))

basic.objects = c(ls(),"basic.objects","bam.file","bam.files","bam.files.mat")


### create TT-seq coverage track ###


if (TRUE){
	dir.create(file.path("TTseqRleTracks"))
	dir.create(file.path("TTseqRleTracks","JaegerWinter2019TTseqKBM7"))
	dir.create(file.path("TTseqRleTracks","JaegerWinter2019TTseqKBM7","UniquePairedNonSplicedFragmentCoverageRleTracks"))

	load(file.path("AnnotationObjects","human.refseq.extended.RData"))
	human.refseq.intron.anno.ranges = as(human.refseq.extended[which(human.refseq.extended[,"type"] == "intron"),],"GRanges")
	
	for (bam.file in bam.files){
		dir.create(file.path("TTseqRleTracks","JaegerWinter2019TTseqKBM7","UniquePairedNonSplicedFragmentCoverageRleTracks",rownames(bam.files.mat[which(bam.files.mat[,"name"] == bam.file),,drop = FALSE])))
		
		registerDoParallel(cores = mc.cores)
		build.jaeger.winter.2019.TTseq.KBM7.unique.paired.non.spliced.fragment.coverage.track.list = function(which.chr){
			jaeger.winter.2019.TTseq.KBM7.unique.paired.non.spliced.fragment.coverage.track.list.chr = list()
			param = ScanBamParam(which=GRanges(seqnames = which.chr,ranges = IRanges(0,human.chrs.lengths[which.chr])))
			bam = readGAlignmentPairsFromBam(file = file.path(prewd,"RawData","JaegerWinter2019TTseqKBM7","BamFiles",bam.file),param = param)
			bam = bam[start(left(bam)) <= end(right(bam))]
			bam.inner.mate = bam[(start(right(bam)) - end(left(bam))) >= 2]

			inner.mate.granges = GRanges(seqnames = which.chr,strand = strand(bam.inner.mate),ranges = IRanges(start = end(left(bam.inner.mate)) + 1,end = start(right(bam.inner.mate)) - 1))
			refseq.intron.anno.ranges = human.refseq.intron.anno.ranges[seqnames(human.refseq.intron.anno.ranges) == which.chr]
			inner.mate.intron.overlaps = findOverlaps(refseq.intron.anno.ranges,inner.mate.granges,maxgap = 0L,minoverlap = 1L,type = "within",select = "all",ignore.strand = FALSE)
			bam.inner.mate = bam.inner.mate[setdiff(1:length(bam.inner.mate),subjectHits(inner.mate.intron.overlaps))]
			
			rle.vec = Rle(0,human.chrs.lengths[which.chr])
			coverage.vec = coverage(bam[strand(bam) == "+"])[[which.chr]]
			rle.vec[1:length(coverage.vec)] = coverage.vec
			coverage.vec = coverage(GRanges(seqnames = which.chr,ranges = IRanges(start = end(left(bam.inner.mate[strand(bam.inner.mate) == "+"])) + 1,end = start(right(bam.inner.mate[strand(bam.inner.mate) == "+"])) - 1)))[[which.chr]]
			rle.vec[1:length(coverage.vec)] = rle.vec[1:length(coverage.vec)] + coverage.vec
			jaeger.winter.2019.TTseq.KBM7.unique.paired.non.spliced.fragment.coverage.track.list.chr[["+"]] = rle.vec
			
			rle.vec = Rle(0,human.chrs.lengths[which.chr])
			coverage.vec = coverage(bam[strand(bam) == "-"])[[which.chr]]
			rle.vec[1:length(coverage.vec)] = coverage.vec
			coverage.vec = coverage(GRanges(seqnames = which.chr,ranges = IRanges(start = end(left(bam.inner.mate[strand(bam.inner.mate) == "-"])) + 1,end = start(right(bam.inner.mate[strand(bam.inner.mate) == "-"])) - 1)))[[which.chr]]
			rle.vec[1:length(coverage.vec)] = rle.vec[1:length(coverage.vec)] + coverage.vec
			jaeger.winter.2019.TTseq.KBM7.unique.paired.non.spliced.fragment.coverage.track.list.chr[["-"]] = rle.vec
			
			save(jaeger.winter.2019.TTseq.KBM7.unique.paired.non.spliced.fragment.coverage.track.list.chr,file = file.path("TTseqRleTracks","JaegerWinter2019TTseqKBM7","UniquePairedNonSplicedFragmentCoverageRleTracks",rownames(bam.files.mat[which(bam.files.mat[,"name"] == bam.file),,drop = FALSE]),paste0("jaeger.winter.2019.TTseq.KBM7.unique.paired.non.spliced.fragment.coverage.track.list.",which.chr,".RData")))
			return()
		}
		
		jaeger.winter.2019.TTseq.KBM7.unique.paired.non.spliced.fragment.coverage.track.list = foreach(n = human.chrs,.noexport = setdiff(ls(),c("human.chrs.lengths"))) %dopar% build.jaeger.winter.2019.TTseq.KBM7.unique.paired.non.spliced.fragment.coverage.track.list(n)
	}
	
	rm(list = setdiff(ls(),basic.objects))
	gc()
}



