require(rtracklayer)

#for loop to process Cbp1 [1] and Cren7 [2] data filtered for 2nd or 1st read 5' occupancy and minimum 100 bp insert size
for(i in 1:3){
  if(i==1){
    #Cbp1
    r1.plus <- import.bw("data/ChIP-exo_S.solfataricus/X1_L76_Cbp1_exo.min100.2ndRead.cov.plus.bw", as="NumericList")
    r1.min <- import.bw("data/ChIP-exo_S.solfataricus/X1_L76_Cbp1_exo.min100.2ndRead.cov.minus.bw", as="NumericList")
    r2.plus <- import.bw("data/ChIP-exo_S.solfataricus/X2_L77_Cbp1_exo.min100.2ndRead.cov.plus.bw", as="NumericList")
    r2.min <- import.bw("data/ChIP-exo_S.solfataricus/X2_L77_Cbp1_exo.min100.2ndRead.cov.minus.bw", as="NumericList")
  }
  if(i==2){
  #CreN7
    r1.plus <- import.bw("data/ChIP-exo_S.solfataricus/X3_L78_Cren7_exo.min100.2ndRead.cov.plus.bw", as="NumericList")
    r1.min <- import.bw("data/ChIP-exo_S.solfataricus/X3_L78_Cren7_exo.min100.2ndRead.cov.minus.bw", as="NumericList")
    r2.plus <- import.bw("data/ChIP-exo_S.solfataricus/X4_L79_Cren7_exo.min100.2ndRead.cov.plus.bw", as="NumericList")
    r2.min <- import.bw("data/ChIP-exo_S.solfataricus/X4_L79_Cren7_exo.min100.2ndRead.cov.minus.bw", as="NumericList")
  }
  if(i==3){
    #Cbp1
    r1.plus <- import.bw("data/ChIP-exo_S.solfataricus/X1_L76_Cbp1_exo.min100.1stRead.cov.plus.bw", as="NumericList")
    r1.min <- import.bw("data/ChIP-exo_S.solfataricus/X1_L76_Cbp1_exo.min100.1stRead.cov.minus.bw", as="NumericList")
    r2.plus <- import.bw("data/ChIP-exo_S.solfataricus/X2_L77_Cbp1_exo.min100.1stRead.cov.plus.bw", as="NumericList")
    r2.min <- import.bw("data/ChIP-exo_S.solfataricus/X2_L77_Cbp1_exo.min100.1stRead.cov.minus.bw", as="NumericList")
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
  #Cbp1 2nd read 
  export.bw(plus.gr, "data/ChIP-exo_S.solfataricus/Cbp1_sso_gm.2ndRead.plus.bw")
  export.bw(minus.gr, "data/ChIP-exo_S.solfataricus/Cbp1_sso_gm.2ndRead.minus.bw")
  }
if(i==2){
  #CreN7 2nd read
  export.bw(plus.gr, "data/ChIP-exo_S.solfataricus/CreN7_sso_gm.2ndRead.plus.bw")
  export.bw(minus.gr, "data/ChIP-exo_S.solfataricus/CreN7_sso_gm.2ndRead.minus.bw")
}
if(i==3){
  #Cbp1 1st read  
  export.bw(plus.gr, "data/ChIP-exo_S.solfataricus/Cbp1_sso_gm.1stRead.plus.bw")
  export.bw(minus.gr, "data/ChIP-exo_S.solfataricus/Cbp1_sso_gm.1stRead.minus.bw")
  }
}
