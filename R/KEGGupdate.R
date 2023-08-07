##' @importMethodsFrom rjson fromJSON
##' @importMethodsFrom magrittr %<>%


#update KO####
updateKO <- function(ver = format(Sys.time(), "%Y%m%d"), dir = system.file(package = "GOfact")) {
  # library(rjson)
  ko.in <- fromJSON(file = "https://www.kegg.jp/kegg-bin/download_htext?htext=br08901.keg&format=json&filedir=")

  ko <- data.frame()
  ko["K", "term"] <- "all"
  ko["K", "id"] <- "K"
  ko["K", "hie"] <- "K"
  ko["K", "level"] <- "0"

  #level 1
  for (i in 1:length(ko.in[["children"]])) {
    name <- ko.in[["children"]][[i]][["name"]]
    hie.1 <- ifelse(i<10,paste0("K", ".0", as.character(i)),paste0("K", ".", as.character(i)))
    ko[hie.1, "term"] <- name
    ko[hie.1, "hie"] <-hie.1
    ko[hie.1, "level"] <- 1
    ko[hie.1, "id"] <- hie.1
    for (j in 1:length(ko.in[["children"]][[i]][["children"]])) {
      name<-
        ko.in[["children"]][[i]][["children"]][[j]][["name"]]
      hie.2<-ifelse(j<10,paste0(hie.1,".0",j),paste0(hie.1,".",j))
      ko[hie.2, "hie"] <- hie.2
      ko[hie.2, "term"] <- name
      ko[hie.2, "level"] <- 2
      ko[hie.2, "id"] <- ko[hie.2, "hie"]
      for (k in 1:length(ko.in[["children"]][[i]][["children"]][[j]][["children"]])) {
        name <-
          ko.in[["children"]][[i]][["children"]][[j]][["children"]][[k]][["name"]]
        hie.3<-ifelse(k<10,paste0(hie.2,".0",k),paste0(hie.2,".",k))
        id.tmp <- unlist(strsplit(name, " "))
        ko[hie.3, "hie"] <- hie.3
        ko[hie.3, "id"] <- id.tmp[1]
        ko[hie.3, "term"] <-
          paste(id.tmp[-c(1, 2)], collapse = " ")
        ko[hie.3, "level"] <- 3
      }
    }
  }

  ko$relation <- "is_a"
  ko$ontology <- "KEGG"
  ko<-ko[,c("id","term","ontology","hie","level","relation")]

  if(!dir.exists(paste0(dir, "/extdata/version/", ver))){
    dir.create(paste0(dir, "/extdata/version/", ver), recursive = TRUE)
  }

  save(ko, file = paste0(dir, "/extdata/version/", ver, "/ko.Rdata"))
  write.csv(ko, file = paste0(dir, "/extdata/version/", ver, "/ko.csv"))


  cat("\n\nKO file updated!\n\n version: ", ver, "\n\n")
  cat("\n\nProgram run finished!\n")
}


#update KOA####
updateKOA <- function(org = "human",
                      ver = format(Sys.time(), "%Y%m%d"),
                      dir = system.file(package = "GOfact")) {
  ### judging species
  if(!(org %in% c("human", "mouse", "rat"))){
    stop("Error for org: ", org, "! only human/mouse/rat allowed\n")
  }

  ### org information conversion
  if (org == "anopheles") {
    species <- "aga"
  } else if (org == "arabidopsis") {
    species <- "ath"
  } else if (org == "bovine") {
    species <- "bta"
  } else if (org == "canine") {
    species <- "cfa"
  } else if (org == "chicken") {
    species <- "gga"
  } else if (org == "chipm") {
    species <- "ptr"
  } else if (org == "ecolik12") {
    species <- "eco"
  } else if (org == "ecsakai") {
    species <- "ecs"
  } else if (org == "fly") {
    species <- "dme"
  } else if (org == "human") {
    species <- "hsa"
  } else if (org == "malaria") {
    species <- "pfa"
  } else if (org == "mouse") {
    species <- "mmu"
  } else if (org == "pig") {
    species <- "ssc"
  } else if (org == "rat") {
    species <- "rno"
  } else if (org == "rhesus") {
    species <- "mcc"
  } else if (org == "worm" || org == "celegans") {
    species <- "cel"
  } else if (org == "xenopus") {
    species <- "xla"
  } else if (org == "yeast") {
    species <- "sce"
  } else if (org == "zebrafish") {
    species <- "dre"
  } else {
    species <- org
  }

  org.lst <- alist()
  org.lst[[1]] <- 'Hs'
  org.lst[[2]] <- 'Mm'
  org.lst[[3]] <- 'Rn'
  names(org.lst) <- c("human", "mouse", "rat")


  ### KEGG api parse
  kegg_api <- paste0("http://rest.kegg.jp/link/", species, "/pathway")
  content <- tryCatch(suppressWarnings(readLines(kegg_api)),
                      error = function(e) NULL)

  # library(magrittr)
  content %<>% strsplit(., "\t") %>% do.call("rbind", .)
  content[,1] %<>% gsub("[^:]+:", "", .)
  content[,1] %<>% gsub(species, "", .)
  content[,2] %<>% gsub("[^:]+:", "", .)
  res <- data.frame(from = content[, 1], to = content[, 2])


  ### term annotation (No parent-child backtracking)
  simple_term2gene <- with(res,
                           split(as.character(to), as.character(from)))
  ### kegg term list & hierachy
  library(rjson)
  ko.in <- fromJSON(file = "https://www.kegg.jp/kegg-bin/download_htext?htext=br08901.keg&format=json&filedir=")
  ko <- data.frame()
  ko["K", "term"] <- "all"
  ko["K", "id"] <- "K"
  ko["K", "hie"] <- "K"
  ko["K", "level"] <- "0"

  for (i in 1:length(ko.in[["children"]])) {
    name <- ko.in[["children"]][[i]][["name"]]
    hie.1 <- ifelse(i<10,paste0("K", ".0", as.character(i)),paste0("K", ".", as.character(i)))
    ko[hie.1, "term"] <- name
    ko[hie.1, "hie"] <-hie.1
    ko[hie.1, "level"] <- 1
    ko[hie.1, "id"] <- hie.1
    for (j in 1:length(ko.in[["children"]][[i]][["children"]])) {
      name<-
        ko.in[["children"]][[i]][["children"]][[j]][["name"]]
      hie.2<-ifelse(j<10,paste0(hie.1,".0",j),paste0(hie.1,".",j))
      ko[hie.2, "hie"] <- hie.2
      ko[hie.2, "term"] <- name
      ko[hie.2, "level"] <- 2
      ko[hie.2, "id"] <- ko[hie.2, "hie"]
      for (k in 1:length(ko.in[["children"]][[i]][["children"]][[j]][["children"]])) {
        name <-
          ko.in[["children"]][[i]][["children"]][[j]][["children"]][[k]][["name"]]
        hie.3<-ifelse(k<10,paste0(hie.2,".0",k),paste0(hie.2,".",k))
        id.tmp <- unlist(strsplit(name, " "))
        ko[hie.3, "hie"] <- hie.3
        ko[hie.3, "id"] <- id.tmp[1]
        ko[hie.3, "term"] <-
          paste(id.tmp[-c(1, 2)], collapse = " ")
        ko[hie.3, "level"] <- 3
      }
    }
  }
  ko<-ko[,c("id","hie")]

  ### get term2gene information
  term2gene <- alist()
  for(i in 1:nrow(ko)){
    n <- grep(ko$hie[i], ko$hie)

    p2c_term2gene <- simple_term2gene[ko$id[n]]
    names(p2c_term2gene) <- NULL
    p2c_term2gene <- unique(unlist(p2c_term2gene))

    term2gene[[i]] <- p2c_term2gene
    if(!is.null(p2c_term2gene)){names(term2gene)[i] <- ko$id[i]}
  }

  term2gene = term2gene[-which(sapply(term2gene, is.null))]
  u.t2Ng <- unlist(lapply(term2gene, length))
  u.t2g <- lapply(term2gene, function(x) paste(x, collapse = ";"))


  ### entrezID to symbols
  OrgDb.name = paste0("org.", org.lst[[org]], ".eg.db")
  library(OrgDb.name, character.only = T)
  OrgDb <- eval(parse(text = OrgDb.name))

  keytype <- "ENTREZID"
  kk <- keys(OrgDb, keytype = keytype)

  goAnno <-
    AnnotationDbi::select(
      OrgDb,
      keys = kk,
      keytype = keytype,
      columns = "SYMBOL"
    )

  df.geneid2genesym <- unique(data.frame(id = goAnno$ENTREZID, sym = goAnno$SYMBOL))
  rownames(df.geneid2genesym) <- df.geneid2genesym$id


  ### save Rdata
  if(!dir.exists(paste0(dir, "/extdata/version/", ver))){
    dir.create(paste0(dir, "/extdata/version/", ver), recursive = TRUE)
  }

  org.koa.file <- paste0(dir, "/extdata/version/", ver, "/",
                         org, ".koa.Rdata")
  save(df.geneid2genesym, u.t2Ng, u.t2g, term2gene, list = c("df.geneid2genesym",
                                                             "u.t2Ng", "u.t2g", "term2gene"), file = org.koa.file)
}
