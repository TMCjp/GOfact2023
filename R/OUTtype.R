##' @importMethodsFrom DT datatable
##' @importMethodsFrom networkD3 saveNetwork
##' @importMethodsFrom xlsx createWorkbook,createSheet,createRow,createCell,setColumnWidth,setCellValue,setCellStyle,addHyperlink,saveWorkbook,CellStyle


#html####
exportHtml <- function(result, htmlfile) {

  Venrich_style <- "<span class='Venrich_type'><b>"
  Nenrich_style <- "<span class='Nenrich_type'><b>"
  Vmiss_style <- "<span class='Vmiss_type'><b>"
  Nmiss_style <- "<span class='Nmiss_type'><b>"
  back_style <- "</b></span>"


  enrichResult <- result
  # enrichResult$pvalue <- signif(enrichResult$pvalue, digits = 3)
  # enrichResult$p.adjust <- signif(enrichResult$p.adjust, digits = 3)
  # enrichResult$qvalue <- signif(enrichResult$qvalue, digits = 3)


  for (i in 1:nrow(enrichResult)) {
    if (enrichResult$P.direction[i] == "++") {
      enrichResult$P.direction[i] = paste0(Venrich_style, enrichResult$P.direction[i], back_style)
    }
    if (enrichResult$P.direction[i] == "+") {
      enrichResult$P.direction[i] = paste0(Nenrich_style, enrichResult$P.direction[i], back_style)
    }

    if (enrichResult$P.direction[i] == "--") {
      enrichResult$P.direction[i] = paste0(Vmiss_style, enrichResult$P.direction[i], back_style)
    }
    if (enrichResult$P.direction[i] == "-") {
      enrichResult$P.direction[i] = paste0(Nmiss_style, enrichResult$P.direction[i], back_style)
    }
    #      enrichResult$Count[i] <- paste0('<a href="#" onclick="alert(\'', enrichResult$geneID[i], '\');">', enrichResult$Count[i], '</a>')
  }

  # enrichResult <-
  #   enrichResult[,-which(colnames(enrichResult) == 'geneID')]

  enrichResult$P.direction <-
    unlist(lapply(enrichResult$P.direction, function(x)
      paste(
        '<style type="text/css"> .Vmiss_type { background: darkgreen;padding: 2px; } .Nmiss_type { background: lightgreen;padding: 2px; } .Venrich_type { background: darkred;padding: 2px; } .Nenrich_type { background: lightcoral;padding: 2px; } </style>',
        x
      )))


  network_text <-
    datatable(
      enrichResult,
      escape = FALSE,
      rownames = FALSE,
      fillContainer = FALSE,
      options = list(
        pageLength = 100,
        autoWidth = TRUE,
        dom = '<"top"lf>rt<"bottom"ip><"clear">'
      )
    )

  saveNetwork(network_text, file = htmlfile, selfcontained = FALSE)
  cat("Success: write HTML file....\nHTML:", htmlfile, "\n")

  # return(network_text)
}



#Excel####
exportExcel <- function(result, xls.file, sheetname, multisheet = FALSE) {

  wb <- createWorkbook(type = "xlsx")


  if(multisheet){

    sumResult <- alist()
    for(i in 1:length(result)){
      sumResult[[i]] <- result[[i]][, c("hie","id","term","p.value","P.direction","p.value.adj","P.adj.direction")]
      colnames(sumResult[[i]])[4:7] <- paste(colnames(sumResult[[i]])[4:7], names(result)[i], sep = "_")
    }
    for(i in 1:(length(sumResult)-1)){
      y <- i + 1
      sumResult[[y]] <- merge(sumResult[[i]], sumResult[[y]], by = c("hie", "id","term"), all = TRUE)
    }
    sumResult <- sumResult[[length(sumResult)]]
    sumResult <- sumResult[, c("hie","id","term", paste("p.value", names(result), sep = "_"), paste("P.direction", names(result), sep = "_"), paste("p.value.adj", names(result), sep = "_"), paste("P.adj.direction", names(result), sep = "_"))]


    ###summary
    sheet <- createSheet(wb, "summary")

    Nrow <- nrow(sumResult)
    Ncol <- ncol(sumResult)
    rows  <- createRow(sheet, rowIndex = 1:(Nrow + 2))
    cells <- createCell(rows, colIndex = 1:Ncol)

    setColumnWidth(sheet, 3, 24)

    col.name <- colnames(sumResult)
    for (cid in 1:Ncol) {
      setCellValue(cells[1, cid][[1]], col.name[cid])
      setCellStyle(cells[1, cid][[1]], style(wb, sig = "0", is.bold = T))
    }

    addMergedRegion(sheet, 1, 2, 1, 1)
    addMergedRegion(sheet, 1, 2, 2, 2)
    addMergedRegion(sheet, 1, 2, 3, 3)
    addMergedRegion(sheet, 1, 1, 4, 3+length(result))
    addMergedRegion(sheet, 1, 1, 4+length(result), 3+2*length(result))
    addMergedRegion(sheet, 1, 1, 4+2*length(result), 3+3*length(result))
    addMergedRegion(sheet, 1, 1, 4+3*length(result), 3+4*length(result))

    setCellValue(cells[1, 4][[1]], "p.value")
    setCellValue(cells[1, 4+length(result)][[1]], "P.direction")
    setCellValue(cells[1, 4+2*length(result)][[1]], "p.value.adj")
    setCellValue(cells[1, 4+3*length(result)][[1]], "P.adj.direction")

    for(i in 4:ncol(sumResult)){
      namevalue <- rep(names(result), 4)
      n <- match(i, 4:ncol(sumResult))
      setCellValue(cells[[2, i]], namevalue[n])
    }



    for (rid in 1:Nrow) {
      for (cid in 1:Ncol) {
        setCellValue(cells[rid + 2, cid][[1]], sumResult[rid, cid])
        if (strsplit(col.name[cid], split = '_')[[1]][1] == "P.direction" |
            strsplit(col.name[cid], split = '_')[[1]][1] == "P.adj.direction") {
          if(!is.na(sumResult[rid, cid])){
            setCellStyle(cells[rid + 2, cid][[1]], style(wb, sig = sumResult[rid, cid]))
          }
        }
      }

      ### KEGG mapping pathway API ###
      if(substr(as.character(sumResult[rid, 2]), start = 1, stop = 2) != "GO"){
        if(!("K" %in% strsplit(as.character(sumResult[rid, 2]), split = '')[[1]])){
          math_1 <- match(sumResult[rid, 2], result[[1]][, 3])
          math_2 <- match(sumResult[rid, 2], result[[2]][, 3])
          math_3 <- match(sumResult[rid, 2], result[[3]][, 3])

          maplist <- paste(result[[1]][math_1, 17], result[[2]][math_2, 17], result[[3]][math_3, 17], sep = "@")
          maplist <- gsub('NA', '', maplist)
          address <- paste0("http://111.198.139.89/pathview?pathway_id=hsa", sumResult[rid, 2], "&symbols=", maplist, "&org=hsa")
          addHyperlink(cells[[rid + 2, 3]], address)
        }
      }
    }


    wb <- gofact_sheet(Workbook = wb, content = result)

  } else {

    wb <- gofact_sheet(wb, sheetname, result)
  }


  saveWorkbook(wb, xls.file)
  cat("Success: write Excel file....\nExcel:", xls.file, "\n")
  # return(xls.file)
}





###set excel style###
style <- function(Workbook, sig = '0', is.bold = F) {
  #set the color of cells
  CellColor <-
    list("#ED2024", "#FDCC99", "#D1E8C5", "#69BD45", "white")
  enr.cla <- c("++", "+", "-", "--", "0")
  if (!sig %in% enr.cla) {
    cat(sig, " is an illeage sig code!")
    stop()
  }
  names(CellColor) <- enr.cla
  cell.style <- CellStyle(Workbook) +
    Font(
      Workbook,
      heightInPoints = 10,
      isBold = is.bold,
      isItalic = F,
      #text font
      name = "Arial",
      color = "#050505"
    ) +
    Fill(
      backgroundColor = "#050505",
      foregroundColor = CellColor[[sig]],
      #cell, last color is important
      pattern = "SOLID_FOREGROUND"
    ) +
    Alignment(h = "ALIGN_CENTER")
  return(cell.style)
}





###creat sheet###
gofact_sheet <- function(Workbook, sheetname, content){

  if(!is.data.frame(content)){

    for(i in 1:length(content)){
      sheet <- createSheet(Workbook, names(content)[i])

      #how to understand the rows and cells!!!!
      Nrow <- nrow(content[[i]])  #with the header for xlsx
      Ncol <- ncol(content[[i]])

      rows  <- createRow(sheet, rowIndex = 1:(Nrow + 1))
      cells <- createCell(rows, colIndex = 1:Ncol)

      #To set the width for working sheets
      # setColumnWidth(sheet, 3, 14)
      setColumnWidth(sheet, 4, 24)


      col.name <- colnames(content[[i]])
      for (cid in 1:Ncol) {
        setCellValue(cells[1, cid][[1]], col.name[cid])
        setCellStyle(cells[1, cid][[1]], style(Workbook, sig = "0", is.bold = T))
      }


      for (rid in 1:Nrow) {
        #    cat(as.character(gofact.result[rid,]))
        for (cid in 1:Ncol) {
          setCellValue(cells[rid + 1, cid][[1]], content[[i]][rid, cid])
          if (col.name[cid] == "P.direction" |
              col.name[cid] == "P.adj.direction") {
            setCellStyle(cells[rid + 1, cid][[1]], style(Workbook, sig = content[[i]][rid, cid]))
          }
        }
        ### KEGG mapping pathway API ###
        if(substr(as.character(content[[i]][rid, 3]), start = 1, stop = 2) != "GO"){
          if(!("K" %in% strsplit(as.character(content[[i]][rid, 3]), split = '')[[1]])){
            maplist <- content[[i]][rid, 17]
            address <- paste0("http://111.198.139.89/pathview?pathway_id=hsa", content[[i]][rid, 3], "&symbols=", maplist, "&org=hsa")
            addHyperlink(cells[[rid + 1, 4]], address)
          }
        }

      }
    }

  } else {

    sheet <- createSheet(Workbook, sheetname)

    Nrow <- nrow(content)
    Ncol <- ncol(content)
    rows  <- createRow(sheet, rowIndex = 1:(Nrow + 1))
    cells <- createCell(rows, colIndex = 1:Ncol)

    setColumnWidth(sheet, 4, 24)

    col.name <- colnames(content)
    for (cid in 1:Ncol) {
      setCellValue(cells[1, cid][[1]], col.name[cid])
      setCellStyle(cells[1, cid][[1]], style(Workbook, sig = "0", is.bold = T))
    }


    for (rid in 1:Nrow) {
      for (cid in 1:Ncol) {
        setCellValue(cells[rid + 1, cid][[1]], content[rid, cid])
        if (col.name[cid] == "P.direction" | col.name[cid] == "P.adj.direction") {
          setCellStyle(cells[rid + 1, cid][[1]], style(Workbook, sig = content[rid, cid]))
        }
      }

      ### KEGG mapping pathway API ###
      if(substr(as.character(content[rid, 3]), start = 1, stop = 2) != "GO"){
        if(!("K" %in% strsplit(as.character(content[rid, 3]), split = '')[[1]])){
          maplist <- content[rid, 17]
          address <- paste0("http://111.198.139.89/pathview?pathway_id=hsa", content[rid, 3], "&symbols=", maplist, "&org=hsa")
          addHyperlink(cells[[rid + 1, 4]], address)
        }
      }
    }



  }

  return(Workbook)
}

