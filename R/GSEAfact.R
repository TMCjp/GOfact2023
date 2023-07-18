##' @importMethodsFrom fgsea fgsea
##' @importMethodsFrom qvalue qvalue

### declarative class
setClass("gseaResult",
         representation   = representation(
           result          = "data.frame",
           organism        = "character",
           setType         = "character",
           geneSets        = "list",
           geneList        = "numeric",
           keytype         = "character",
           permScores      = "matrix",
           params          = "list",
           gene2Symbol     = "character",
           readable        = "logical"
         )
)


calculate_qvalue <- function(pvals){
  qobj <- tryCatch(qvalue(pvals, lambda = 0.05, pi0.method = "bootstrap"),
                   error = function(e) NULL)
  if (class(qobj) == "qvalue") {
    qvalues <- qobj$qvalues
  }
  else {
    qvalues <- NA
  }
  return(qvalues)
}


GSEAfact <- function(geneList,
                     org = "human",
                     ont = "KEGG",
                     ver = "latest",
                     dir = system.file(package = "GOfact"),
                     exponent = 1,   #weight of each step
                     nPerm = 1000,   #permutation numbers
                     minGSSize = 1,
                     maxGSSize = 20000,
                     pvalueCutoff = 0.05,
                     pAdjustMethod = "BH",
                     verbose = TRUE,
                     seed=FALSE) {

  #Check parameter
  if (ver == "latest") {
    dirs <- dir(paste0(dir, "/extdata/version"))
    #      dirs<-append(dirs,"20190422")
    dirs <- sort(dirs, decreasing = T)
    ver <- dirs[1]
  }
  else if (!all(grepl("^\\d{8}$", ver, perl = "T"))) {
    stop(
      "\nError! GOfact version error for ver:",
      ver,
      ". It should be with the similar format as: 20190421!\n"
    )
  }
  if (!dir.exists(paste0(dir, "/extdata/version/", ver))) {
    stop(
      "\nError! GOfact version set error. No dir for:",
      paste0(dir, "/extdata/", ver),
      ". Dir not exists!\n"
    )
  }



  if(verbose)
    message("preparing geneSet collections...")

  if(ont == "MSigDB"){
    ### import msigdb
    if(org == "human"){specie = "Homo sapiens"}
    if(org == "mouse"){specie = "Mus musculus"}
    if(org == "rat"){specie = "Rattus norvegicus"}
    m = msigdbr(species = specie)

    df.gene2go <- unique(data.frame(gene = m$entrez_gene, go = m$gs_id))
    term2gene <- with(df.gene2go, split(as.character(gene), as.character(go)))

    m <- m[, c('gs_id', 'gs_cat', 'gs_subcat')]
    m$hie <- paste(m$gs_cat, m$gs_subcat, sep = ':')
    m <- m[, c('gs_id', 'gs_subcat', 'hie')]
    ontology_hie <- m
    colnames(ontology_hie) <- c('id', 'gs_subcat', 'hie')
    rm(list = 'm')

  } else if(ont == "KEGG"){
    ontology <- "ko"
    oa.Rdata.file <- paste0(dir, "/extdata/version/", ver, "/", org, ".", ontology, "a.Rdata")
    load(oa.Rdata.file)

    o.Rdata.file <- paste0(dir, "/extdata/version/", ver, "/", ontology, ".Rdata")
    load(o.Rdata.file)
    ontology_hie <- ko
    rm(list = 'ko')

  } else {
    #Three cases to be added to GO
    ontology <- "go"
    oa.Rdata.file <- paste0(dir, "/extdata/version/", ver, "/", org, ".", ontology, "a.Rdata")
    load(oa.Rdata.file)

    o.Rdata.file <- paste0(dir, "/extdata/version/", ver, "/", ontology, ".Rdata")
    load(o.Rdata.file)
    ontology_hie <- go
    colnames(ontology_hie)[2] <- 'id'
    rm(list = 'go')

  }

  if(verbose)
    message("GSEA analysis...")

  tmp_res <- fgsea(pathways=term2gene,
                   stats=geneList,
                   nperm=nPerm,
                   minSize=minGSSize,
                   maxSize=maxGSSize,
                   gseaParam=exponent,
                   nproc = 0)

  p.adj <- p.adjust(tmp_res$pval, method=pAdjustMethod)
  qvalues <- calculate_qvalue(tmp_res$pval)

  # Description <- TERM2NAME(tmp_res$pathway, USER_DATA)

  params <- list(pvalueCutoff = pvalueCutoff,
                 nPerm = nPerm,
                 pAdjustMethod = pAdjustMethod,
                 exponent = exponent,
                 minGSSize = minGSSize,
                 maxGSSize = maxGSSize
  )

  res <- data.frame(
    ID = as.character(tmp_res$pathway),
    # Description = Description,
    setSize = tmp_res$size,
    enrichmentScore = tmp_res$ES,
    NES = tmp_res$NES,
    pvalue = tmp_res$pval,
    p.adjust = p.adj,
    qvalues = qvalues,
    stringsAsFactors = FALSE
  )

  res <- res[!is.na(res$pvalue),]
  res <- res[ res$pvalue <= pvalueCutoff, ]
  res <- res[ res$p.adjust <= pvalueCutoff, ]
  idx <- order(res$pvalue, decreasing = FALSE)
  res <- res[idx, ]

  ### add hie view
  idx <- na.omit(match(res$ID, ontology_hie$id))
  if('gs_subcat' %in% colnames(ontology_hie)){
    n <- which(ontology_hie$gs_subcat[idx] == '')
    ontology_hie$hie[idx][n] <- gsub(':', '', ontology_hie$hie[idx][n], fixed = T)
  }
  res$hie <- ontology_hie$hie[idx]


  if (nrow(res) == 0) {
    message("no term enriched under specific pvalueCutoff...")
    # return(
    #   new("gseaResult",
    #       result     = res,
    #       geneSets   = term2gene,
    #       geneList   = geneList,
    #       params     = params,
    #       readable   = FALSE
    #   )
    # )
  }

  # row.names(res) <- res$ID
  # observed_info <- lapply(term2gene[res$ID], function(gs)
  #   gseaScores(geneSet=gs,
  #              geneList=geneList,
  #              exponent=exponent)
  # )
  #
  # if (verbose)
  #   message("leading edge analysis...")
  #
  # ledge <- leading_edge(observed_info)
  #
  # res$rank <- ledge$rank
  # res$leading_edge <- ledge$leading_edge
  # res$core_enrichment <- sapply(ledge$core_enrichment, paste0, collapse='/')

  if (verbose)
    message("done...")

  new("gseaResult",
      result     = res,
      geneSets   = term2gene,
      geneList   = geneList,
      params     = params,
      readable   = FALSE
  )
}
