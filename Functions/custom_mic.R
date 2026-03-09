

## Microbial ------------------------------------------------------------------------------------------
## Functions to index in phloseq objects

otu_mat <- function (ps)
  as(otu_table(ps), "matrix")

tax_mat <- function (ps)
  as(tax_table(ps), "matrix")
sample_df <- function (ps)
  as(sample_data(ps), "data.frame")

## Function to import physeq into metagenomeSeq, modified
## (original gave errors, matrix conversion issues)
physeq_to_metagenomeSeq_mod <- function (physeq, ...) 
{
  if (!taxa_are_rows(physeq)) {
    physeq <- t(physeq)
  }
  countData = round(as(otu_table(physeq), "matrix"), digits = 0)
  if (!is.null(sample_data(physeq, FALSE))) {
    ADF = AnnotatedDataFrame(data.frame(sample_data(physeq)))
  }
  else {
    ADF = NULL
  }
  if (!is.null(tax_table(physeq, FALSE))) {
    TDF = AnnotatedDataFrame(data.frame(OTUname = taxa_names(physeq), 
                                        data.frame(as(tax_table(physeq), "matrix")), 
                                        row.names = taxa_names(physeq)))
  }
  else {
    TDF = AnnotatedDataFrame(data.frame(OTUname = taxa_names(physeq), 
                                        row.names = taxa_names(physeq)))
  }
  if (requireNamespace("metagenomeSeq")) {
    mrobj = metagenomeSeq::newMRexperiment(counts = countData, 
                                           phenoData = ADF, featureData = TDF, ...)
    if (sum(colSums(countData > 0) > 1) < ncol(countData)) {
      p = suppressMessages(metagenomeSeq::cumNormStat(mrobj))
    }
    else {
      p = suppressMessages(metagenomeSeq::cumNormStatFast(mrobj))
    }
    mrobj = metagenomeSeq::cumNorm(mrobj, p = p)
    return(mrobj)
  }
}

dada2.rrndb <- function(x, rrndb) {
  if (all(names(x) %in% levels(factor(rrndb$rank)) ==  F)) {
    stop("Rank names do not match.")
  }
  
  
  sub <- x %>% select_if(~ !any(is.na(.)))
  if(ncol(sub) == 0L){
    out <- data.frame(
      OTU = row.names(x),
      x,
      class.rank = "no.match",
      copy.number = NA,
      stringsAsFactors = F
    )
  } else {
    #names(sub) <- names(x)[1:length(sub)]
    class.rank <- colnames(sub)[ncol(sub)]
    match.rrndb <-
      rrndb[rrndb$rank == class.rank &
              rrndb$name == sub[,paste(class.rank)],]

    if (nrow(match.rrndb) == 0L) {
      while (nrow(match.rrndb) == 0L) {
        class.rank <- colnames(sub)[which(colnames(sub) == class.rank) - 1]
        match.rrndb <-
          rrndb[rrndb$rank == class.rank &
                  rrndb$name == sub[,paste(class.rank)],]
        
        if(class.rank == "domain" & nrow(match.rrndb) == 0L){
          out <- data.frame(
            OTU = row.names(x),
            x,
            class.rank = "no.match",
            copy.number = NA,
            stringsAsFactors = F
          )
          break
        }
      }
    }
  }
  
  if(nrow(match.rrndb) > 0L){
    out <- data.frame(
      OTU = row.names(x),
      x,
      class.rank = match.rrndb$rank,
      copy.number = match.rrndb$mean,
      stringsAsFactors = F
    )
  }
  
  return(out)
}

dada2_match_rrndb <- function(tax, rrndb, .parallel = F) {
  if (.parallel == T) {
    if (exists("cl") == F) {
      stop(
        "Define cluster for parallel computing with `detectCores()` and `makeCluster()`."
      )
    }
    
    do.call("rbind",
            pbapply::pbapply(
              cl = cl,
              X = tax,
              MARGIN = 1,
              FUN = dada2.rrndb,
              rrndb = rrndb
            ))
  } else {
    # normal version
    do.call("rbind",
            pbapply::pbapply(
              X = tax,
              MARGIN = 1,
              FUN = dada2.rrndb,
              rrndb = rrndb
            ))
  }
}


## Big data functions -------------------------------------------------------------------------------------------

minhead <- function(x){x[1:5,1:5]}


