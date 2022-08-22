## ----include=FALSE------------------------------------------------------------
NOT_CRAN <- identical(tolower(Sys.getenv("NOT_CRAN")), "true")
knitr::opts_chunk$set(purl = NOT_CRAN)

## ----eval=FALSE---------------------------------------------------------------
#  install.packages("geneHapR")

## ----startup, eval=NOT_CRAN, echo=TRUE, message=FALSE, warning=FALSE----------
library(geneHapR)

## ----setups, eval=NOT_CRAN, include=FALSE-------------------------------------
data("geneHapR_test")

## ----import vcf, eval=FALSE, include=TRUE-------------------------------------
#  # import vcf file
#  vcf <- import_vcf("your_vcf_file_path.vcf")
#  
#  # import gziped vcf file
#  vcf <- import_vcf("your_vcf_file_path.vcf.gz")

## ----import gff, eval=FALSE, include=TRUE-------------------------------------
#  # import GFFs
#  gff <- import_gff("your_gff_file_path.gff", format = "GFF")

## ----import DNA seqs, eval=FALSE, include=TRUE--------------------------------
#  # import DNA sequences in fasta format
#  seqs <- import_seqs("your_DNA_seq_file_path.fa", format = "fasta")

## ----import phenotype and accession groups, eval=FALSE, include=TRUE----------
#  # import phynotype data
#  pheno <- import_AccINFO("your_pheno_file_path.txt")
#  pheno

## ----eval=NOT_CRAN, echo=FALSE------------------------------------------------
head(pheno)

## ----import pheno, eval=FALSE-------------------------------------------------
#  
#  # import accession group/location information
#  AccINFO <- import_AccINFO("accession_group_file_path.txt")

## ----eval=NOT_CRAN, echo=FALSE------------------------------------------------
head(AccINFO)

## ----import pheno and accessions_groups, eval=FALSE, include=TRUE-------------
#  # import pheno from space ' ' delimed table
#  pheno <- read.table("your_pheno_file_path.csv", header = TRUE, row.names = 1, comment.char = "#")
#  
#  # import pheno from ',' delimed table
#  pheno <- read.csv("your_pheno_file_path.csv", header = TRUE, comment.char = "#")

## ----VCF filtration, eval=FALSE, include=TRUE---------------------------------
#  # filter VCF by position
#  vcf_f1 <- filter_vcf(vcf, mode = "POS",
#                       Chr = "scaffold_1",
#                       start = 4300, end = 5890)
#  
#  # filter VCF by annotation
#  vcf_f2 <- filter_vcf(vcf, mode = "type",
#                       gff = gff,
#                       type = "CDS")
#  
#  # filter VCF by position and annotation
#  vcf_f3 <- filter_vcf(vcf, mode = "both",
#                       Chr = "scaffold_1",
#                       start = 4300, end = 5890,
#                       gff = gff,
#                       type = "CDS")

## ----preprocess large VCF, eval=FALSE, include=TRUE---------------------------
#  # new VCF file will be saved to disk
#  # extract a single gene/range from a large vcf
#  
#  filterLargeVCF(VCFin = "Ori.vcf.gz",
#                 VCFout = "filtered.vcf.gz",
#                 Chr = "scaffold_8",
#                 POS = c(19802,24501),
#                 override = TRUE)
#  
#  # extract multi genes/ranges from large vcf
#  filterLargeVCF(VCFin = "Ori.vcf.gz",          # surfix should be .vcf.gz or .vcf
#                 VCFout = c("filtered1.vcf.gz", # surfix should be .vcf.gz or .vcf
#                            "filtered2.vcf.gz",
#                            "filtered3.vcf.gz"),
#                 Chr = c("scaffold_8",
#                         "scaffold_8",
#                         "scaffold_7"),
#                 POS = list(c(19802,24501),
#                            c(27341,28949),
#                            c(38469,40344)),
#                 override = TRUE)               # if TRUE, existed file will be override without warning

## ----DNA alignment and trim, eval=FALSE, message=FALSE, warning=FALSE, include=TRUE----
#  # sequences alignment
#  seqs <- allignSeqs(seqs, quiet = TRUE)
#  
#  # sequences trim
#  seqs <- trimSeqs(seqs, minFlankFraction = 0.1)
#  seqs

## ----haplotype calculation from vcf, eval=NOT_CRAN, include=TRUE, message=FALSE, warning=FALSE----
hapResult <- vcf2hap(vcf,
                     hapPrefix = "H",
                     hyb_remove = TRUE,
                     na.drop = TRUE)
hapResult

## ----haplotype calculation from seqs, eval=FALSE------------------------------
#  hapResult <- seqs2hap(seqs,
#                        Ref = names(seqs)[1],
#                        hapPrefix = "H",
#                        hyb_remove = TRUE,
#                        na.drop = TRUE,
#                        maxGapsPerSeq = 0.25)

## ----check site numbers, eval=NOT_CRAN----------------------------------------
# Chech number of sites conclude in hapResult
sites(hapResult)

## ----replace INFO, eval=NOT_CRAN, include=TRUE, message=FALSE, warning=FALSE----
# add annotations to INFO field
hapResult <- addINFO(hapResult,
                     tag = "PrChange",
                     values = rep(c("C->D", "V->R", "G->N"),3),
                     replace = TRUE)

## ----add INFO, eval=NOT_CRAN, include=TRUE, message=FALSE, warning=FALSE------
# To replace the origin INFO by set 'replace' as TRUE
hapResult <- addINFO(hapResult,
                     tag = "CDSChange",
                     values = rep(c("C->A", "T->C", "G->T"),3),
                     replace = FALSE)

## ----ATG position, eval=FALSE-------------------------------------------------
#  # set ATG position as zero in gff
#  newgff <- gffSetATGas0(gff = gff, hap = hapResult,
#                         geneID = "test1G0387",
#                         Chr = "scaffold_1",
#                         POS = c(4300, 7910))
#  
#  # set position of ATG as zero in hapResult/hapSummary
#  newhap <- hapSetATGas0(gff = gff, hap = hapResult,
#                         geneID = "test1G0387",
#                         Chr = "scaffold_1",
#                         POS = c(4300, 7910))

## ----hapSummary, eval=NOT_CRAN, fig.height=4, fig.width=6, message=FALSE, warning=FALSE, paged.print=FALSE----
hapSummary <- hap_summary(hapResult)
hapSummary

## ---- eval=NOT_CRAN, fig.height=4, fig.width=6, message=FALSE, warning=FALSE, paged.print=FALSE----
plotHapTable(hapSummary)

## ---- eval=NOT_CRAN, fig.height=4, fig.width=6, message=FALSE, warning=FALSE, paged.print=FALSE----
# add one annotation
plotHapTable(hapSummary,
             hapPrefix = "H",
             INFO_tag = "CDSChange", 
             tag_name = "CDS",
             displayIndelSize = 1, 
             angle = 45,
             replaceMultiAllele = TRUE,
             ALLELE.color = "grey90")

## ---- eval=NOT_CRAN, fig.height=4, fig.width=6, message=FALSE, warning=FALSE, paged.print=FALSE----
# add multi annotation
plotHapTable(hapSummary,
             hapPrefix = "H",
             INFO_tag = c("CDSChange", "PrChange"),
             displayIndelSize = 1, 
             angle = 45,
             replaceMultiAllele = TRUE,
             ALLELE.color = "grey90")

## ---- eval=NOT_CRAN, fig.height=4, fig.width=6, message=FALSE, warning=FALSE, paged.print=FALSE----
# add multi annotation
plotHapTable(hapSummary,
             hapPrefix = "H",
             INFO_tag = c("CDSChange", "PrChange"),
             tag_name = c("CDS", "Pr"),
             displayIndelSize = 1, 
             angle = 45,
             replaceMultiAllele = TRUE,
             ALLELE.color = "grey90")

## ----display variations on gene model, eval=NOT_CRAN, fig.height=4, fig.width=6, message=FALSE, warning=FALSE, paged.print=FALSE----
displayVarOnGeneModel(gff,
                      hapSummary,
                      Chr = "scaffold_1",
                      startPOS = 4300, endPOS = 7910,
                      type = "pin", cex = 0.7,
                      CDS_h = 0.05, fiveUTR_h = 0.02, threeUTR_h = 0.01)

## ----hapNet, eval=NOT_CRAN----------------------------------------------------
hapNet <- get_hapNet(hapSummary,
                     AccINFO = AccINFO,
                     groupName = "Type")

## ---- eval=NOT_CRAN, fig.height=4, fig.width=6, message=FALSE, warning=FALSE, paged.print=FALSE----
# plot haploNet
plotHapNet(hapNet,
           size = "freq",                   # circle size
           scale = "log2",                 # scale circle with 'log10(size + 1)'
           cex = 0.8,                       # size of hap symbol
           col.link = 2,                    # link colors
           link.width = 2,                  # link widths
           show.mutation = 2,               # mutation types one of c(0,1,2,3)
           legend = c(-13,-2))        # legend position

## ----hapVsPheno merged, eval=NOT_CRAN, fig.height=4, fig.width=10, message=FALSE, warning=FALSE, paged.print=FALSE----
results <-hapVsPheno(hapResult,
                     hapPrefix = "H",
                     title = "This is title",
                     mergeFigs = TRUE,
                     pheno = pheno,
                     phenoName = "GrainWeight.2021",
                     minAcc = 3)
plot(results$figs)

## ----hapVsPheno separated, eval=NOT_CRAN, fig.height=4, fig.width=6, message=FALSE, warning=FALSE, paged.print=FALSE----
results <- hapVsPheno(hap = hapResult,
                      hapPrefix = "H",
                      title = "This is title",
                      pheno = pheno,
                      phenoName = "GrainWeight.2021",
                      minAcc = 3,
                      mergeFigs = FALSE)
plot(results$fig_pvalue)
plot(results$fig_Violin)

## ----eval=FALSE---------------------------------------------------------------
#  hapVsPhenos(hapResult,
#              pheno,
#              outPutSingleFile = TRUE,
#              hapPrefix = "H",
#              title = "Seita.0G000000",
#              file = "mypheno.tiff",
#              width = 12,
#              height = 8,
#              res = 300)

## ----geo distribution, eval=NOT_CRAN, fig.height=4, fig.width=6, message=FALSE, warning=FALSE, paged.print=FALSE----
# library(mapdata)
# library(maptools)
hapDistribution(hapResult,
                AccINFO = AccINFO,
                LON.col = "longitude",
                LAT.col = "latitude", 
                hapNames = c("H001", "H002", "H003"), 
                legend = TRUE)
