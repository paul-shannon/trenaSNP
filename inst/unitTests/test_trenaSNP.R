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
test_createModel <- function()
{
   printf("--- test_with.snp.table")

   targetGene <- "MEF2C"
   regionsOfInterestSpec <- "allDNA"
   tsnp <- trenaSNP("hg38", tbl.variants)
   checkTrue(is(tsnp, "trenaSNP"))
   checkEquals(sort(getVariantTable(tsnp)$id), c("rs13158247", "rs244761"))
   #x <- createModel(tsnp, targetGene,

} # test_with.snp.table
#------------------------------------------------------------------------------------------------------------------------
