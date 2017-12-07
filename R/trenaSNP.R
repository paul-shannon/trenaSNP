.trenaSNP <- setClass ("trenaSNP",
                    representation = representation(
                       genomeName="character",
                       tbl.variants="data.frame",
                       state="environment",                # for dynamic slot data
                       trena="Trena",
                       tv="trenaViz",
                       quiet="logical")
                    )

#------------------------------------------------------------------------------------------------------------------------
setValidity('trenaSNP', function(object)
{
   msg = NULL

      # any number of columns are supported, but these three are required
      # the id should be and rsid or, eg, "chrN:loc:ref:alt"

   missing.columns.in.variants.table <- setdiff(c("id", "chrom", "pos"), colnames(object@tbl.variants))
   if(length(missing.columns.in.variants.table) > 0)
      msg <- c(msg, sprintf("missing columns in variant table: %s",
                            paste(missing.columns.in.variants.table, collapse=",")))
   if (is.null(msg))
      return(TRUE)
   else
      return(msg)

}) # setValditity
#------------------------------------------------------------------------------------------------------------------------
setGeneric('getVariantTable',  signature='obj', function(obj) standardGeneric("getVariantTable"))
setGeneric('createModel',      signature='obj', function(obj, modelName, targetGene, pfms, motifMatchThreshold,
                                                         expressionMatrix, chrom, start, end)
                                       standardGeneric("createModel"))
#------------------------------------------------------------------------------------------------------------------------
trenaSNP = function(genomeName, tbl.variants, trenaVizPortRange=12000:12030, quiet=TRUE)
{
    stopifnot(genomeName %in% c("hg38", "mm10"))

    state <- new.env(parent=.GlobalEnv)
    trena <- Trena(genomeName)
    tv <- trenaViz(trenaVizPortRange)

    obj <- .trenaSNP(genomeName=genomeName,
                     tbl.variants=tbl.variants,
                     trena=trena,
                     tv=tv,
                     state=state,
                     quiet=quiet)

    obj

} # constructor
#------------------------------------------------------------------------------------------------------------------------
setMethod('getVariantTable', 'trenaSNP',

  function (obj) {
     obj@tbl.variants
     })

#------------------------------------------------------------------------------------------------------------------------
setMethod('createModel', 'trenaSNP',

  function (obj, modelName, targetGene, pfms, motifMatchThreshold, expressionMatrix, chrom, start, end) {

     })

#------------------------------------------------------------------------------------------------------------------------
#.modelWithParameters <- function(pfms, tbl.regulatoryRegions, motifMatchThreshold, mtx)
#{
#   tbl.motifsInRegulatoryRegions  <- findMotifs(pfms, tbl.regulatoryRegions, motifMatchThreshold)
#   tbl.motifs <- associateTranscriptionFactors(MotifDb, tbl.motifsInRegulatoryRegions, source="MotifDb", expand.rows=TRUE)
#   tbl.motifs$geneSymbol <- toupper(tbl.motifs$geneSymbol)
#
#   #solver.names <- c("lasso", "pearson", "randomForest", "ridge", "spearman", "sqrtlasso", "lassopv")
#   solver.names <- c("lasso", "pearson", "randomForest", "ridge", "spearman")
#   candidates <- unique(intersect(tbl.motifs$geneSymbol, rownames(mtx)))
#   #save(targetGene, mtx, candidates, solver.names, file=sprintf("solver.bug.%s.RData", gsub(" ", ".", Sys.time(), fixed=TRUE)))
#   #suppressWarnings(
#      tbl.geneModel <- createGeneModel(trena, targetGene, solver.names, tbl.motifs, mtx)
#   #   )
#
#   if(nrow(tbl.geneModel) > 0)
#      tbl.geneModel <- tbl.geneModel[order(tbl.geneModel$rfScore, decreasing=TRUE),]
#
#   targetGene.tss <- with(tbl.promoters, ifelse(strand[1] == "+",
#                                                max(transcription_start_site),
#                                                min(transcription_start_site)))
#
#   tbl.motifs.strong <- subset(tbl.motifs, geneSymbol %in% tbl.geneModel$gene)
#   distance <- tbl.motifs.strong$motifStart - targetGene.tss
#   direction <- rep("upstream", length(distance))
#   direction[which(distance < 0)] <- "downstream"
#   tbl.motifs.strong$distance.from.tss <- distance
#   tbl.motifs.strong$id <- sprintf("%s.fp.%s.%06d.%s", targetGene, direction, abs(distance), tbl.motifs.strong$motifName)
#
#   list(model=tbl.geneModel, regions=tbl.motifs.strong)
#
#} # modelWithParameters
#------------------------------------------------------------------------------------------------------------------------
