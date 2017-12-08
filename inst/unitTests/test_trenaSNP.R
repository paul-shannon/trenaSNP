library(trenaSNP)
library(RUnit)
#------------------------------------------------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
#------------------------------------------------------------------------------------------------------------------------
if(!exists("mtx")){
   load("~/s/work/priceLab/cory/brainExpressionData/mtx.rosmap.normalized.RData")
   mtx <- mtx.rosmap.normalized
   }
if(!exists("tbl.variants"))
   tbl.variants <- read.table(system.file(package="trenaSNP", "extdata", "mef2c.2snps.tsv"),
                              sep="\t", header=TRUE, as.is=TRUE)
#------------------------------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_constructor()
   test_createModel.allDNA()

} # runTests
#------------------------------------------------------------------------------------------------------------------------
test_constructor <- function()
{
   printf("--- test_constructor")
   targetGene <- "MEF2C"
   colnames(tbl.variants)[1] <- "broken"
   checkException(
      tsnp <- trenaSNP("hg38", tbl.variants),
      "caught!", silent=TRUE)

   colnames(tbl.variants)[1] <- "id"
   regionsOfInterestSpec <- "allDNA"
   tsnp <- trenaSNP("hg38", tbl.variants)
   checkTrue(is(tsnp, "trenaSNP"))

   regionsOfInterestSpec <- c("allDNA", "encodeHumanDHS",
                              "postgres://bddsrds.globusgenomics.org/skin_wellington_20")
   tsnp <- trenaSNP("hg38", tbl.variants)
   checkTrue(is(tsnp, "trenaSNP"))
   checkEquals(sort(getVariantTable(tsnp)$id), c("rs13158247", "rs244761"))

} # test_constructor
#------------------------------------------------------------------------------------------------------------------------
test_createModel.allDNA <- function()
{
   printf("--- test_createModel.allDNA")

   targetGene <- "MEF2C"
   targetGene.tss <- 88904257

   motifSource <- "allDNA"

   tsnp <- trenaSNP("hg38", tbl.variants)
   showGenomicRegion(tsnp, targetGene)
   displayVariants(tsnp)
   pfms <- as.list(query(query(MotifDb, "jaspar2016"), "hsapiens"))

   modelName <- "MEF2C-allDNA-90"
   showGenomicRegion(tsnp, "chr5:88,734,414-88,735,443")

   roi <- parseChromLocString(getGenomicRegion(tsnp))
   x <- createModel(tsnp, modelName, targetGene, targetGene.tss,
                    roi, motifSource,
                    pfms, pfmMatchThreshold=90, mtx,
                    motifToTFmapper="MotifDb", display=TRUE)
   checkEquals(names(x), c("model", "regions"))
   tbl.model <- x$model
   tbl.regions <- x$regions

      # this cherry-picked region has cherry-picked snp rs13158247, at chr5:88734853, IUPAC_CODE_MAP[["Y"]] ->  "CT"
   tbl.snps <- tbl.variants[, c("chrom", "pos", "pos")]
   colnames(tbl.snps) <- c("chrom", "start", "end")
   tbl.overlaps <- as.data.frame(findOverlaps(GRanges(tbl.snps[, 1:3]), GRanges(tbl.regions[, 2:4])))
   colnames(tbl.overlaps) <- c("snp", "fp")
   tbl.snpFps <- cbind(tbl.variants[tbl.overlaps$snp,], tbl.regions[tbl.overlaps$fp,])
   checkEquals(tbl.snpFps$geneSymbol, c("FOXC1", "FOXC1"))
   checkEquals(tbl.snpFps$shortMotif, c("MA0032.1", "MA0032.1"))
   checkEquals(tbl.snpFps$id,         c("rs13158247", "rs13158247"))
   pfm.ma0032.1 <- as.list(query(query(MotifDb, "MA0032.1"), "jaspar2016"))
   tbl.damage <- assessSnp(tsnp@trena, pfm.ma0032.1, "rs13158247", 7, pwmMatchMinimumAsPercentage=85)

     #  MA0032.1   chr5 88734851 88734858      -     4.3750          0.9090909 TACTGAAA
     # rs13158247  chr5 88734853

} # test_createModel.allDNA
#------------------------------------------------------------------------------------------------------------------------
test_createModel.brain_hint_20 <- function()
{
   printf("--- test_createModel.brain_hint_20")

   targetGene <- "MEF2C"
   targetGene.tss <- 88904257

   motifSouce <-  "postgres://bddsrds.globusgenomics.org/brain_hint_20"

   tsnp <- trenaSNP("hg38", tbl.variants)
   showGenomicRegion(tsnp, "MEF2C")
   displayVariants(tsnp)
   pfms <- as.list(query(query(MotifDb, "jaspar2016"), "hsapiens"))

   modelName <- "MEF2C-allDNA-90"
   showGenomicRegion(tsnp, "chr5:88,734,414-88,735,443")

   roi <- parseChromLocString(getGenomicRegion(tsnp))
   x <- createModel(tsnp, modelName, targetGene, roi, motifSource,
                            pfms, pfmMatchThreshold=90,
                            mtx,
                            motifToTFmapper="MotifDb",
                            display=TRUE)

   #x <- createModel(tsnp, targetGene

} # test_createModel.brain_hint_20
#------------------------------------------------------------------------------------------------------------------------
