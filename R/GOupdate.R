##' @importMethodsFrom DBI dbGetQuery
##' @importMethodsFrom GO.db GO_dbconn


#update GO
updateGO <- function(ver = format(Sys.time(), "%Y%m%d"), dir = system.file(package = "GOfact")) {

  ptc1 <- proc.time()
  # installPackages(c('DBI', 'GO.db'))
  # library(DBI)
  library(GO.db)
  GOterm <- dbGetQuery(GO_dbconn(), "SELECT * FROM go_term")
  rownames(GOterm) <- GOterm$go_id

  go <- GOterm[-which(GOterm$go_id == 'all'), ]
  go$hie <- NA
  go$level <- NA
  go$relation <- NA
  GOroots <- c("GO:0008150", "GO:0005575", "GO:0003674")
  go[GOroots, 'hie'] <- c("BP", "CC", "MF")
  go[GOroots, 'relation'] <- c("root", "root", "root")
  go[GOroots, 'level'] <- 1

  #  GOroots<-GOroots[2]
  Nontology <- length(GOroots)

  go2child <-
    c(as.list(GOCCCHILDREN),
      as.list(GOBPCHILDREN),
      as.list(GOMFCHILDREN))

  next.goids <- GOroots
  Level <- 2
  while (length(next.goids) >= 1) {
    cat("Now analyze Level:",
        Level,
        " ( ",
        Nontology,
        " ontologies )\n")
    new.next.goids <- character(0)
    for (ng in next.goids) {
      if (is.na(ng)) {
        stop("one of the next GO term is NA!\n", ng)
      }

      goids <- go2child[[ng]]

      goids <-
        goids[names(goids) == 'is_a' | names(goids) == 'part_of']

      goids <- goids[is.na(go[goids, 'hie'])]

      if (length(goids) > 0) {

        new.next.goids <- c(new.next.goids, goids)

        goterms <- go[goids, 'term']
        goids <- goids[order(goterms)]

        if (length(goids) == 0) {
          stop("Error in children of goid, no children!")
        } else if (length(goids) < 10) {
          code <- paste0("0", 1:length(goids))
        } else{
          code <- c(paste0("0", 1:9), 10:length(goids))
        }

        hies <-
          paste0(go[ng, 'hie'], ".", code)
        go[goids, 'relation'] <- names(goids)
        go[goids, 'hie'] <- hies
        go[goids, 'level'] <- Level
      }
      else{
        next
      }
    }
    cat("Now level",
        Level,
        " for #GOterm",
        length(next.goids),
        "finished!\n")

    next.goids <- new.next.goids
    Level = Level + 1
  }

  if(!dir.exists(paste0(dir, "/extdata/version/", ver))){
    dir.create(paste0(dir, "/extdata/version/", ver), recursive = TRUE)
  }
  save(go, file = paste0(dir, "/extdata/version/", ver, "/go.Rdata"))
  write.csv(go, file = paste0(dir, "/extdata/version/", ver, "/go.csv"))
  cat("\n\nGO file updated!\n\n version: ", ver, "\n\n")
  cat("Time passed", (proc.time() - ptc1)[1], " seconds\n")
  cat("\n\nProgram run finished!\n")
}



#updateGOfact####
updateGOA <- function (org = "human", ver = format(Sys.time(), "%Y%m%d"),
                       dir = system.file(package = "GOfact"))
{
  org.lst <- alist()
  org.lst[[1]] <- 'Hs'
  org.lst[[2]] <- 'Mm'
  org.lst[[3]] <- 'Rn'
  names(org.lst) <- c("human", "mouse", "rat")


  if (is.null(org.lst[[org]])) {
    stop("Update GOfact set error for org:", org, " not allowed!\n")
  }
  if (!all(grepl("^\\d{8}$", ver, perl = "T"))) {
    stop("Update GOfact set error for ver:", ver, ". It should be 20190421!\n")
  }
  if (!dir.exists(dir)) {
    stop("Update GOfact set error for dir:", dir, ". Not exists!\n")
  }
  version.dir <- paste0(dir, "/extdata/version/", ver)
  if (!dir.exists(version.dir)) {
    if (!dir.create(version.dir, recursive = TRUE)) {
      stop(paste0("dir.create error for ", version.dir))
    }
    else {
      cat(paste0("dir.created for ", version.dir))
    }
  }


  if (org %in% c("human", "mouse", "rat")) {
    cat("\nNow begin to process the GOA data (internet needed!). \n           Need 110 or 5 seconds(130 or 240 seconds in total)......\n")
    ptc0 <- proc.time()
    goa.tmp.Rdata.file <- paste0(dir, "/extdata/version/",
                                 ver, "/tmp.goa.Rdata")
    need.download <- T
    if (need.download) {
      OrgDb.name = paste0("org.", org.lst[[org]], ".eg.db")
      library(OrgDb.name, character.only = T)
      OrgDb <- eval(parse(text = OrgDb.name))
      keytype <- "ENTREZID"
      kk <- keys(OrgDb, keytype = keytype)
      goAnno <- select(OrgDb, keys = kk, keytype = keytype,
                       columns = c("GOALL", "ONTOLOGYALL", "SYMBOL"))
      goAnno <- unique(goAnno[!is.na(goAnno$GOALL), ])
      save(goAnno, file = goa.tmp.Rdata.file)
      cat("\nGOA download and parse success:", (proc.time() -
                                                  ptc0)[1], " seconds\n")
      cat("\nTime version is", toString(Sys.time()), "\n")
    }
    else {
      if (!exists(goa.tmp.Rdata.file)) {
        stop("No goa tmp Rdata file! change the set!")
      }
      load(goa.tmp.Rdata.file)
      cat("\ngoAnno Rdata file read finished:", (proc.time() -
                                                   ptc0)[1], " seconds\n")
    }
    cat("\nNow begin to ananlyze the GOA datasets! About 130 seconds.....\n")
    ptc1 <- proc.time()
    df.gene2go <- unique(data.frame(gene = goAnno$ENTREZID,
                                    go = goAnno$GOALL))
    df.geneid2genesym <- unique(data.frame(id = goAnno$ENTREZID,
                                           sym = goAnno$SYMBOL))
    rownames(df.geneid2genesym) <- df.geneid2genesym$id
    colnames(df.gene2go)
    term2gene <- split(as.character(df.gene2go$gene), df.gene2go$go)
    u.t2Ng <- unlist(lapply(term2gene, length))
    u.t2g <- lapply(term2gene, function(x) paste(x, collapse = ";"))
    cat("\nProcess goa for:", (proc.time() - ptc1)[1], " seconds\n")
  }
  else {
    stop("Error for org: ", org, "! only human/mouse/rat allowed\n")
  }

  org.goa.file <- paste0(dir, "/extdata/version/", ver, "/",
                         org, ".goa.Rdata")
  save(df.geneid2genesym, u.t2Ng, u.t2g, term2gene, list = c("df.geneid2genesym",
                                                             "u.t2Ng", "u.t2g", "term2gene"), file = org.goa.file)
}
