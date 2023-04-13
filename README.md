# Cbp1_public

This repository contains code and processed data related out study of the archaeal chromatin protein Cbp1 entitled:
Cbp1-Cren7 chromatinization of CRISPR arrays favours transcription from leader- over cryptic promoters.

The preprint of the manuscript is available here:<br />
[@ResearchSquare](https://www.researchsquare.com/article/rs-2781205/v1)<br />
[@Biorxiv](https://www.biorxiv.org/content/10.1101/2023.03.24.534125v1)<br />
<br />
While [processed data for individual replicates](https://github.com/fblombach/Cbp1_public/tree/main/data) are included in this repository, some of the code
relies on alignment files (bam and bed format) that are too large to be included. The raw sequencing data from which these bam files were generated are available at
 NCBI GEO under superseries [GSE226026](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE226026)<br />
<br />
Code is mostly in R markdown format including R and bash code chunks. Where relevant, the html output for the R markdown is included. 
The html output cannot be visualised directly in github.
<br />

# Replicate handling
Basic code for the combining of replicates for Cappable-seq and ChIP-exo data and normalisation is available here:<br />
[Cappable-seq data for *S. islandicus* REY15A](https://github.com/fblombach/Cbp1_public/blob/main/scripts/cappable-seq/mergingReplicates.Sisl.R)<br />
[Cappable-seq data for *S. solfataricus* P2](https://github.com/fblombach/Cbp1_public/blob/main/scripts/cappable-seq/mergingReplicates.Sso.R)<br />
[Cbp1 and Cren7 ChIP-exo data *S. solfataricus* P2](https://github.com/fblombach/Cbp1_public/blob/main/scripts/mergingReplicates.ChIP-exo.R)<br />
<br />

# Peak calling
MACS2 peak calling for Cbp1 ChIP-seq data is documented here:<br />
[*S. solfataricus* P2](https://github.com/fblombach/Cbp1_public/blob/main/scripts/PeakCalling_SsoCbp1.Rmd)<br />
[*S. islandicus* REY15A](https://github.com/fblombach/Cbp1_public/blob/main/scripts/PeakCalling_SislCbp1.Rmd)<br />
[*S. islandicus* LAL14/1](https://github.com/fblombach/Cbp1_public/blob/main/scripts/PeakCalling_SisLAL141_Cbp1.Rmd)<br />
<br />

# Genome-wide comparison of Cbp1 and Cren7 ChIP-seq occupancy in *S. solfataricus* P2
[R markdown](https://github.com/fblombach/Cbp1_public/blob/main/scripts/Cbp1_Cren7_occupancy.Rmd)

# Identification of a CRISPR repeat-like binding motif in non-canonical Cbp1 binding sites
[*S. solfataricus* P2](https://github.com/fblombach/Cbp1_public/blob/main/scripts/nonCRISPR_peaks.Sso.Rmd)<br />
[*S. islandicus* REY15A](https://github.com/fblombach/Cbp1_public/blob/main/scripts/nonCRISPR_peaks.SislREY15A.Rmd)<br />

# Differential binding analysis for Cren7 in *S.islandicus* REY15A *cbp1* deletion strain
[R markdown](https://github.com/fblombach/Cbp1_public/blob/main/scripts/DiffBind_CreN7.Rmd) and [html output](https://github.com/fblombach/Cbp1_public/blob/main/scripts/DiffBind_CreN7.html)

# Aggregate plots of Cbp1 and Cren7 ChIP-exo profiles on CRISRP A/B and non-canonical binding sites in *S. solfataricus* P2
[CRISPR A/B repeats](https://github.com/fblombach/Cbp1_public/blob/main/scripts/metaplots.chip-exo.Rmd) 
and [html output](https://github.com/fblombach/Cbp1_public/blob/main/scripts/metaplots.chip-exo.html)<br />
[non-canonical binding sites](https://github.com/fblombach/Cbp1_public/blob/main/scripts/metaplots.chip-exo.nonCRISPR.Rmd)
and [html output](https://github.com/fblombach/Cbp1_public/blob/main/scripts/metaplots.chip-exo.nonCRISPR.html)<br />

# Differential expression analysis for Cappable-seq data from *S.islandicus* REY15A *cbp1* deletion strain vs WT
[R markdown](https://github.com/fblombach/Cbp1_public/blob/main/scripts/TSScallingFromCappableSeqData.v2.Rmd) and
[html output](https://github.com/fblombach/Cbp1_public/blob/main/scripts/TSScallingFromCappableSeqData.v2.html)

# The alignment of *IS*110 family transposons in *S. solfataricus* P2
[R markdown](https://github.com/fblombach/Cbp1_public/blob/main/scripts/ISC1229bindingByCbp1.v2.Rmd) and
[html output](https://github.com/fblombach/Cbp1_public/blob/main/scripts/ISC1229bindingByCbp1.v2.html)
