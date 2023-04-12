require(rtracklayer)

#for loop to process Cbp1 [1] and Cren7 [2] data
for(i in 1:2){
  if(i==1){
    #Cbp1
    r1.plus <- import.bw("data/ChIP-exo_S.solfataricus/GSM7061727_X1_L76_Cbp1_exo.plus.bw", as="NumericList")
    r1.min <- import.bw("data/ChIP-exo_S.solfataricus/GSM7061727_X1_L76_Cbp1_exo.minus.bw", as="NumericList")
    r2.plus <- import.bw("data/ChIP-exo_S.solfataricus/GSM7061728_X2_L77_Cbp1_exo.plus.bw", as="NumericList")
    r2.min <- import.bw("data/ChIP-exo_S.solfataricus/GSM7061728_X2_L77_Cbp1_exo.minus.bw", as="NumericList")
  }
  if(i==2){
  #CreN7
    r1.plus <- import.bw("data/ChIP-exo_S.solfataricus/GSM7061729_X3_L78_CreN7_exo.plus.bw", as="NumericList")
    r1.min <- import.bw("data/ChIP-exo_S.solfataricus/GSM7061729_X3_L78_CreN7_exo.minus.bw", as="NumericList")
    r2.plus <- import.bw("data/ChIP-exo_S.solfataricus/GSM7061730_X4_L79_CreN7_exo.plus.bw", as="NumericList")
    r2.min <- import.bw("data/ChIP-exo_S.solfataricus/GSM7061730_X4_L79_CreN7_exo.minus.bw", as="NumericList")
  }
  
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
  
if(i==1){
  #Cbp1  
  export.bw(plus.gr, "data/ChIP-exo_S.solfataricus/Cbp1_sso_gm.plus.bw")
  export.bw(minus.gr, "data/ChIP-exo_S.solfataricus/Cbp1_sso_gm.minus.bw")
  }
if(i==2){
  #CreN7
  export.bw(plus.gr, "data/ChIP-exo_S.solfataricus/CreN7_sso_gm.plus.bw")
  export.bw(minus.gr, "data/ChIP-exo_S.solfataricus/CreN7_sso_gm.minus.bw")
  }
}
