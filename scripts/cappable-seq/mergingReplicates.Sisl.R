require(rtracklayer)

samples <- c("E234","cbp1")

for(i in 1:2){
  r1.plus <- import.bw(paste0("data/Cappable-seq_S.islandicus/", samples[i], "_r1.5primeCov.plus.bw"), as="NumericList")
  r1.min <- import.bw(paste0("data/Cappable-seq_S.islandicus/", samples[i], "_r1.5primeCov.minus.bw"), as="NumericList")
  r2.plus <- import.bw(paste0("data/Cappable-seq_S.islandicus/", samples[i], "_r2.5primeCov.plus.bw"), as="NumericList")
  r2.min <- import.bw(paste0("data/Cappable-seq_S.islandicus/", samples[i], "_r2.5primeCov.minus.bw"), as="NumericList")

  #calculate geometric mean
  plus.gm <- unlist(sqrt(r1.plus * r2.plus))
  minus.gm <- unlist(sqrt(r1.min * r2.min))

  #scale to tpm values
  sc.fac <- 10^6/(sum(plus.gm) + sum(minus.gm))

  plus.gr <- makeGRangesFromDataFrame(data.frame(score=plus.gm*sc.fac, 
                                                start=1:length(plus.gm), 
                                                end=1:length(plus.gm),
                                                seqname = rep(seqinfo(r1.plus)@seqnames, times=length(plus.gm))),
                                      ignore.strand = T, 
                                      start.field="start", 
                                      end.field="end",
                                      seqnames.field = "seqname",
                                      keep.extra.columns = T
                                      )
  minus.gr <- makeGRangesFromDataFrame(data.frame(score=minus.gm*sc.fac, 
                                                 start=1:length(minus.gm), 
                                                 end=1:length(minus.gm),
                                                 seqname = rep(seqinfo(r1.min)@seqnames, times=length(minus.gm))),
                                      ignore.strand = T, 
                                      start.field="start", 
                                      end.field="end",
                                      seqnames.field = "seqname",
                                      keep.extra.columns = T
                                      )
  seqlengths(plus.gr) <- length(plus.gm)
  seqlengths(minus.gr) <- length(minus.gm)
  export.bw(plus.gr, paste0("data/Cappable-seq_S.islandicus/", samples[i], "_gm.5primeCov.plus.bw"))
  export.bw(minus.gr, paste0("data/Cappable-seq_S.islandicus/", samples[i], "_gm.5primeCov.minus.bw"))
}
