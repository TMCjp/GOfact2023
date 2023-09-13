#gofact#####
gofact <- function(gene,
           slim = "generic",
           org = "human",
           ver = "latest",
           dir = system.file(package = "GOfact")){
    options(digits = 4)

    # #org
    # if (is.null(org.lst[[org]])) {
    #   stop("\nError! GOfact run error for org:", org, " not allowed!\n")
    # }
    #ver
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
    #dir
    if (!dir.exists(dir)) {
      stop("\nError! GOfact error for dir:", dir, ". Not exists!\n")
    }

    TimeName <- gsub(":", "-", sub("\\s", "\\_", Sys.time()))
    cat(
      paste0(
        "\n\n--Welcome to use the GOfact2019. The time stamp is :",
        TimeName,
        ".\n     Program is running......."
      )
    )
    ptm <- proc.time()


    ### if length(gene) > 1, gene must is character; but if length(gene) == 1, gene maybe character/file.
    ### so you must verification. date:20190701
    if (is.data.frame(gene)) {
      stop("\nGene file is data.frame\n")
    } else if (length(gene) > 1) {
      cat("\nGene input is character array!\n")
      genes = gene
      gene = "Vector"
      cat("\n", length(genes), " genes read sucess!\n")
    } else {
      if (file.exists(gene)) {
        genefile = gene
        genes = readLines(genefile)
        cat("\nGene file: ",
            genefile,
            " read success (with ",
            length(genes),
            " genes)!\n")
      } else if (length(grep(pattern = '^[0-9]*$', gene, perl = T)) > 0){   #Regular matching number, entrezID
        cat("\nGene input is character array!\n")
        genes = gene
        gene = "Vector"
        cat("\n", length(genes), " genes read sucess!\n")
      } else {
        stop("\nFile: ", toString(gene), "open error!\n")
      }
    }


    # ontology <- "go"
    if (slim == "KEGG" | slim == 'kegg') {
      cat("\n\n----Will call KEGG analyses")
      ontology <- "ko"
    } else {
      ontology <- "go"
    }


    cat("\n\n--Now load the ",ontology," ontology Rdata......\n")
    ptm1 <- proc.time()
    o.Rdata.file <-
      paste0(dir, "/extdata/version/", ver, "/", ontology, ".Rdata")
    oa.Rdata.file <-
      paste0(dir, "/extdata/version/", ver, "/", org, ".", ontology, "a.Rdata")

    if (file.exists(o.Rdata.file)) {
      load(o.Rdata.file)
    } else{
      stop("\n\nError!!!\n There is still no ",
           ontology,
           " data for version:",
           ver,
           "\n\n")
    }
    if (file.exists(oa.Rdata.file)) {
      load(oa.Rdata.file)
    } else{
      stop("\n\nError!!!\n There is still no ", ontology, "a data for version:", ver, "\n\n")
    }

    pt <- proc.time() - ptm1
    cat("\nLoad",ontology," ontology Rdata success. Time: ", pt[1], " seconds!\n")


    cat("\nNow load the gene datasets\n")
    genes <- unique(genes)



    cat("\n--Now read ", ontology, " slim files......")
    if (slim == "KEGG" | slim == "kegg") {
      slim.all<-ko
    } else{

      if (all(grepl("^GO(\\d+)$", slim, ignore.case = T) |
              all(grepl("^all$", slim, ignore.case = T)))) {
        cat("\nGOfact will use the GO hierarchy for slim:",
            toupper(slim))
      } else{
        #slim files(user files!)
        slim.file <- paste0(dir, "/extdata/slim/", slim, ".slim")
        if (!file.exists(slim.file)) {
          stop("\n\n\nError!!!! Slim file read error: =>",
               slim.file,
               "<=!")
        } else{
          slim.r <- readLines(slim.file)
          #      slim.r<-slim.r[grep("^#.*",slim.r)]
          slims <- grep("^GO:\\d{7}$",
                        slim.r,
                        perl = T,
                        value = T)

          if (length(slim.r) != length(slims)) {
            stop(
              "\n\n\n\nDuring slim making,\nerror GO term format in slim files.\n\t=>",
              slim.file,
              ".\n\n\n\n"
            )
          }
          else{
            cat(
              "\n\nSlim file read success read success (with ",
              length(slims),
              " GO terms):",
              slim.file,
              "!\n"
            )
          }
        }
      }


      slim.Rdata <-
        paste0(dir, "/extdata/version/", ver, "/", slim, "-slim.all.Rdata")
      #go.slim <- paste0("GO:000", 1000:1200)
      if (!file.exists(slim.Rdata)) {
        cat("\nSlim: ", slim.Rdata, " not exists!  Now create......\n")
        ptm1 <- proc.time()
        if (all(grepl("^GO\\d+$", slim, ignore.case = T)) |
            (all(grepl("^all$", slim, ignore.case = T)))) {
          if ((all(grepl("^GO\\d+$", slim, ignore.case = T)))) {
            N.cutoff <- sub("GO", "", slim, ignore.case = T)
            N.cutoff <- as.integer(N.cutoff)
          } else{
            N.cutoff <- 1000
          }
          slims <- go[go$level <= N.cutoff, "go_id"]
        }

        slim.leave <-
          go[which(go$go_id %in% slims), c("hie", "go_id", "term", "ontology")]

        go$hie.c <- paste0(gsub("\\.", "#", go$hie), "#")
        slim.leave$hie.c <-
          paste0(gsub("\\.", "#", slim.leave$hie), "#")
        slim.all <-
          go[sapply(go$hie.c, function(x)
            length(grep(x, slim.leave$hie.c)) > 0), c("hie", "go_id", "term", "ontology")]
        slim.all$id<-slim.all$go_id
        slim.all$go_id<-NULL
        go$hie.c <- NULL
        save(slim.all, file = slim.Rdata)
        pt <- proc.time() - ptm1
        cat("\nFinished!  Time:", pt[1], " seconds\n")
      } else{
        load(slim.Rdata)
        cat("\nDirectly use the slim.Rdata: ", slim.Rdata, "\n")
      }
    }

    cat("\n\nNow begin the data analysis! About 4 seconds.......")

    result <- slim.all[order(slim.all$hie), ]
    rownames(result) <- result$id


    result$all_genes <- u.t2g[result$id]
    result$nall <- u.t2Ng[result$id]



    result <- result[which(result$nall > 0),]


    #ontology root
    if(slim=="KEGG"|slim=="kegg"){
      o2g = list(KEGG = "K", kegg = "K")
      u.o2g <- unlist(o2g)

    }else{
      o2g = list(BP = "GO:0008150", CC = "GO:0005575", MF = "GO:0003674")
      u.o2g <- unlist(o2g)
    }


    ptm <- proc.time()

    Qterm2gene <-
      lapply(term2gene, function(x)
        x[x %in% genes])

    u.N.Qterm2gene <-
      unlist(lapply(Qterm2gene, length))

    Qterm2gene = Qterm2gene[which(u.N.Qterm2gene != 0)]

    u.N.Qterm2gene <- u.N.Qterm2gene[u.N.Qterm2gene != 0]

    Qterm2gene.id <-
      lapply(Qterm2gene, function(x)
        paste(x, collapse = ";"))

    u.Qterm2gene.id <- unlist(Qterm2gene.id)


    Qterm2gene.sym <-
      lapply(Qterm2gene, function(x)
        paste(df.geneid2genesym[x, 'sym'], collapse = ";"))

    u.Qterm2gene.sym <- unlist(Qterm2gene.sym)
    u.N.Qterm2gene <- u.N.Qterm2gene[u.N.Qterm2gene > 0]

    cat(
      "\nClustering finished. Time passed: ",
      as.character(proc.time() - ptm)[1],
      "seconds!\n"
    )


    result <-
      result[result$id %in% names(u.N.Qterm2gene), ]

    result$npro <- u.N.Qterm2gene[result$id]
    result$genes.id <- u.Qterm2gene.id[result$id]
    result$genes.sym <- u.Qterm2gene.sym[result$id]

    result$Npro <-  result[u.o2g[result$ontology], 'npro']
    result$Nall <-  result[u.o2g[result$ontology], 'nall']



    result$ratio <- round(result$npro / result$Npro, digits = 4)
    result$bg.ratio <-
      round(result$nall / result$Nall, digits = 4)
    # result$ER <- round(result$ratio / result$bg.ratio, digits = 4)
    result$ER <- round((result$npro / result$Npro) / (result$nall / result$Nall), digits = 4)
    result$exp <-
      round(result$nall / result$Nall * result$Npro, digits = 4)



    for (goid in result$id) {
      npro <- result[goid, 'npro']
      nall <- result[goid, 'nall']
      Nall <- result[goid, 'Nall']
      Npro <- result[goid, 'Npro']
      Exp <- result[goid, 'exp']
      ER <- result[goid, 'ER']
      p.h <-
        phyper(npro - 1, nall, Nall - nall, Npro, lower.tail = F)
      p.h <- round(p.h, digits = 4)
      if (ER >= 1) {
        result[goid, 'p.value'] <- p.h
      } else{
        result[goid, 'p.value'] <- 1 - p.h
      }
      if ((ER > 1) &
          (result[goid, 'p.value'] < 0.05)) {
        result[goid, 'P.direction'] = "++"
      } else if ((ER > 1) &
                 (result[goid, 'p.value'] >= 0.05)) {
        result[goid, 'P.direction'] = "+"
      } else if ((ER < 1) &
                 (result[goid, 'p.value'] < 0.05)) {
        result[goid, 'P.direction'] = "--"
      } else if ((ER < 1) &
                 (result[goid, 'p.value'] >= 0.05)) {
        result[goid, 'P.direction'] = "-"
      } else{
        result[goid, 'P.direction'] = "0"
      }
    }

    p.adj <-
      p.adjust(result$p.value, method = 'BH')  #BH
    p.adj <- round(p.adj, digits = 4)
    result$p.value.adj <- p.adj


    for (goid in result$id) {
      ER <- result[goid, 'ER']
      if ((ER > 1) & (result[goid, 'p.value.adj'] < 0.05)) {
        result[goid, 'P.adj.direction'] = "++"
      } else if ((ER > 1) &
                 (result[goid, 'p.value.adj'] >= 0.05)) {
        result[goid, 'P.adj.direction'] = "+"
      } else if ((ER < 1) &
                 (result[goid, 'p.value.adj'] < 0.05)) {
        result[goid, 'P.adj.direction'] = "--"
      } else if ((ER < 1) &
                 (result[goid, 'p.value.adj'] >= 0.05)) {
        result[goid, 'P.adj.direction'] = "-"
      } else{
        result[goid, 'P.adj.direction'] = "0"
      }
    }

    ### add direction, if enrichment result is null! date:20190701
    #if(!("P.direction" %in% colnames(result))){
    if(nrow(result) < 1){
      result[1, ] <- NA
      result$p.value <- NA
      result$P.direction <- NA
      result$P.adj.direction <- NA
      result$genes.id <- NA
      result$genes.sym <- NA
      result <- na.omit(result)
    }


    result <-
      result[, c(
        "ontology",
        "hie",
        "id",
        "term",
        "npro",
        "ratio",
        "ER",
        "exp",
        "p.value",
        "P.direction",
        "p.value.adj",
        "P.adj.direction",
        "nall",
        "bg.ratio",
        "Npro",
        "Nall",
        "genes.id",
        "genes.sym"
      )]

    cat(
      "\nProgram finished.  Time passed: ",
      as.character(proc.time() - ptm)[1],
      "seconds!\n\n\n"
    )

    return(result)
  }




#multicol gofact#####
gofact_MultiCol <- function(csvfile,
                            slim = "generic",
                            org = "human",
                            ver = "latest",
                            dir = system.file(package = "GOfact")){

  ###csv header names
  genematrix <- read.table(csvfile, sep = ',', header = T, na.strings = NA)

  ngenelist <- ncol(genematrix)
  result <- alist()
  for(i in 1:ngenelist){
    cat('run genelist ', i)
    genelist <- as.character(na.omit(genematrix[, i]))
    result[[i]] <- gofact(genelist, slim = slim, org = org, ver = ver, dir = dir)
  }
  names(result) <- colnames(genematrix)

  return(result)
}


