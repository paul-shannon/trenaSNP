.trenaSNP <- setClass ("trenaSNP",
                    representation = representation(
                       genomeName="character",
                       tbl.variants="data.frame",
                       state="environment",                # for dynamic slot data
                       trena="Trena",
                       trenaViz="trenaViz",
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
setGeneric('displayVariants',  signature='obj', function(obj, title="variants", color="red") standardGeneric("displayVariants"))
setGeneric('createModel',      signature='obj', function(obj, modelName, targetGene, targetGene.tss,
                                                         roi, motifSources, pfms, pfmMatchThreshold=90,
                                                         mtx, motifToTFmapper="MotifDb", display=TRUE)
                                   standardGeneric("createModel"))
setGeneric('findMotifs',       signature='obj', function(obj, pfms, tbl.regions, pwmMatchMinimumAsPercentage,
                                                         source, display=FALSE, trackName=NA_character_)
                                   standardGeneric('findMotifs'))
#------------------------------------------------------------------------------------------------------------------------
trenaSNP = function(genomeName, tbl.variants, trenaVizPortRange=12000:12030, quiet=TRUE)
{
    stopifnot(genomeName %in% c("hg38", "mm10"))

    state <- new.env(parent=.GlobalEnv)
    trena <- Trena(genomeName)
    trenaViz <- trenaViz(trenaVizPortRange)
    setGenome(trenaViz, genomeName)

    obj <- .trenaSNP(genomeName=genomeName,
                     tbl.variants=tbl.variants,
                     trena=trena,
                     trenaViz=trenaViz,
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
setMethod('displayVariants', 'trenaSNP',

  function (obj, title="variants", color="red") {
     tbl <- obj@tbl.variants
     tbl.bed <- tbl[, c("chrom", "pos", "pos", "id")]
     colnames(tbl.bed) <- NULL
     removeTracksByName(obj@trenaViz, title)
     addBedTrackFromDataFrame(obj@trenaViz, title, tbl.bed, displayMode="squished", color=color, trackHeight=50)
     })

#------------------------------------------------------------------------------------------------------------------------
setMethod('showGenomicRegion', 'trenaSNP',

  function (obj, regionString) {
     showGenomicRegion(obj@trenaViz, regionString)
     })

#------------------------------------------------------------------------------------------------------------------------
setMethod('getGenomicRegion', 'trenaSNP',

  function (obj) {
     getGenomicRegion(obj@trenaViz)
     })

#------------------------------------------------------------------------------------------------------------------------
setMethod('createModel', 'trenaSNP',

  function(obj, modelName, targetGene, targetGene.tss, roi, motifSources, pfms, pfmMatchThreshold=90,
           mtx, motifToTFmapper="MotifDb", display=TRUE){

     browser()
     stopifnot(.recognizedMotifSource(motifSources))
     tbl.regions <- as.data.frame(roi, stringsAsFactors=FALSE)

     motif.list <- list()

     for(motifSource in motifSources){

        if(motifSource == "allDNA"){
          tbl.newMotifs <- findMotifs(obj, pfms, tbl.regions, pfmMatchThreshold, motifSource, display, modelName)
          if(nrow(tbl.newMotifs) == 0){
             printf(" trenaSNP::createModel: no motifs found in region");
             } # no motifs
          else{
             motif.list[[motifSource]] <- tbl.newMotifs
             }
          } # allDNA

        if(grepl("postgres://", motifSource)){
           browser()
           tbl.newMotifs <- getRegulatoryChromosomalRegions(obj@trena, roi$chrom, roi$start, roi$end, motifSource, targetGene, tss)
           motif.list[[motifSource]] <- tbl.newMotifs
           }
        } # for motifSource

     browser()
     tbl.motifs <- do.call(rbind, motif.list)
     colnames(tbl.motifs)[grep("motifStart", colnames(tbl.motifs))] <- "start"
     colnames(tbl.motifs)[grep("motifEnd", colnames(tbl.motifs))] <- "end"
     solver.names <- c("lasso", "pearson", "randomForest", "ridge", "spearman")
     candidates <- unique(intersect(tbl.motifs$geneSymbol, rownames(mtx)))
     tbl.geneModel <- createGeneModel(obj@trena, targetGene, solver.names, tbl.motifs, mtx)
     tbl.geneModel.strong <- subset(tbl.geneModel, rfScore > 1 | pcaMax > 1)
     tbl.motifs.strong <- subset(tbl.motifs, geneSymbol %in% tbl.geneModel.strong$gene)
     browser()
     if(display){
        addBedTrackFromDataFrame(tsnp@trenaViz, "motifsInModel",
                                 tbl.motifs.strong[, c("chrom", "start", "end", "motifName", "motifRelativeScore")],
                                 displayMode="squished", color="lightGreen", trackHeight=30);
        }
     list(model=tbl.geneModel.strong, regions=tbl.motifs.strong)
     }) # createModel

#------------------------------------------------------------------------------------------------------------------------
setMethod('findMotifs', 'trenaSNP',

   function(obj, pfms, tbl.regions, pwmMatchMinimumAsPercentage, source, display=FALSE, trackName=NA_character_){
      mm <- MotifMatcher(genomeName="hg38", pfms)
      tbl.regions.uniq <- unique(tbl.regions[, 1:3])
      tbl.motifs <- findMatchesByChromosomalRegion(mm, tbl.regions.uniq, pwmMatchMinimumAsPercentage=pwmMatchMinimumAsPercentage)
      if(nrow(tbl.motifs) == 0){
         printf("--- no match of pfms in supplied regions at %d%%", pwmMatchMinimumAsPercentage)
         return(data.frame())
         }

      shortMotifs <- unlist(lapply(strsplit(tbl.motifs$motifName, "-"), function(tokens) tokens[length(tokens)]))
      tbl.motifs$shortMotif <- shortMotifs
      tbl.motifs <- associateTranscriptionFactors(MotifDb, tbl.motifs, source="MotifDb", expand.rows=TRUE)

      motifs.without.tfs <- which(is.na(tbl.motifs$geneSymbol))

      if(length(motifs.without.tfs) > 0){
         printf("%d/%d motifs had no TF/geneSymbol, removing", length(motifs.without.tfs), nrow(tbl.motifs))
         tbl.motifs <- tbl.motifs[-motifs.without.tfs,]
         }


      motif.tfs <- sort(unique(tbl.motifs$geneSymbol))
      if(display){
         addBedTrackFromDataFrame(obj@trenaViz, trackName,
                                  tbl.motifs[, c("chrom", "motifStart", "motifEnd", "motifName", "motifRelativeScore")],
                                  color="darkGreen")
         }

      invisible(tbl.motifs)

      }) # findMotifs
#----------------------------------------------------------------------------------------------------
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
.recognizedMotifSource <- function(motifSources)
{
   all.recognized <- TRUE;  # be optimisitc

   for(motifSource in motifSources){
      if(motifSource == "allDNA")
         next;
      if(grepl("postgres://", motifSource))
         next;
      all.recognized <- FALSE;
      printf("  unrecognized motifSource: %s", motifSource)
      }

   return(all.recognized)

} # .recognizedMotifSource
#------------------------------------------------------------------------------------------------------------------------


