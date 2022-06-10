library(shiny)
library(ggplot2)
library(ggrepel)
library(nashvillePlot)

options(shiny.maxRequestSize = 1024^3)

## Read a specific metaxcan file from a zip file
processMetaxcanFile <- function(file, zipFile) {
    dat <- read.csv(unz(zipFile, file))
    dat$dtype <- file
    return(dat)
}

## Create the nessary dataset for the nashville plot
processMetaxcanFiles <- function(zipFile, mapVersion){
    ## Get list of metaxcan files contained in zip file and only select those that are csv files
    files <- unzip(zipFile, list=TRUE)$Name
    files <- labels <- files[grepl("\\.csv$", files)]

    map <- switch(mapVersion,
                  `37` = gene.build.37,
                  `38` = gene.build.38,
                  stop("error"))

    ## Create a dataset with all the sub-metaxcan dataset from the zip file combined
    metaxcan <- do.call(rbind, lapply(files, processMetaxcanFile, zipFile=zipFile))

    names(metaxcan)[names(metaxcan) == "gene"] <- "ENSG"

    return(merge(metaxcan, map, by="ENSG"))
}

function(input, output) {
    labelsProcess <- reactive({
        labelFile <- input$labels

        cat("Trying Labels Processing ...\n")
        if(is.null(labelFile))
            return(NULL)
        cat("Processing Labels File ...\n")

        return(read.table(labelFile$datapath, header=FALSE, stringsAsFactors=FALSE))
    }) |>
    bindCache(input$labels)
    
    gwasProcess <- reactive({
        gwasFile <- input$gwas

        cat("Trying GWAS Processing ...\n")
        req(gwasFile)
        cat("Processing GWAS File ...\n")

        gwas <- read.table(unz(gwasFile$datapath, unzip(gwasFile$datapath, list=TRUE)$Name[1]), header= TRUE)

        return(make.valid.object(CHR = gwas$CHR,
                                 P = gwas$P,
                                 BP = gwas$BP))
    }) |>
    bindCache(input$gwas)

    metaProcess <- reactive({
        metaxcanFile <- input$metaxcan
        mapVersion <- input$version

        cat("Trying metaXcan processing ...\n")
        req(metaxcanFile, mapVersion)
        cat("Processing metaXcan File ...\n")

        meta <- processMetaxcanFiles(metaxcanFile$datapath, mapVersion)
        
        return(make.valid.object(CHR = as.numeric(meta$CHR),
                                 P = meta$pvalue,
                                 BP = meta$MID_POS,
                                 gene.start = meta$START_POS,
                                 gene.end = meta$END_POS,
                                 group = meta$dtype,
                                 gene.name = meta$Gene))
    }) |>
    bindCache(input$metaxcan, input$version)

    plotGenomeWide <- reactive({
        gwas <- gwasProcess()
        metaxcan <- metaProcess()
        mapVersion <- input$version
        labels <- labelsProcess()
        
        return(nashville.plot(data1=gwas,
                              data2=metaxcan,
                              map_df=mapVersion,
                              config=labels))
    }) |>
    bindCache(gwasProcess, metaProcess, input$version, labelsProcess)

    plotChromosome <- reactive({
        gwas <- gwasProcess()
        metaxcan <- metaProcess()
        mapVersion <- input$version
        labels <- labelsProcess()
        chr <- input$chr
        
        return(nashville.plot(data1=gwas,
                              data2=metaxcan,
                              map_df=mapVersion,
                              config=labels,
                              chr=chr))
    }) |>
    bindCache(gwasProcess, metaProcess, input$version, labelsProcess, input$chr)
    
    plotGene <- reactive({
        gwas <- gwasProcess()
        metaxcan <- metaProcess()
        mapVersion <- input$version
        labels <- labelsProcess()
        gene <- input$gene
        geneEnsg <- input$gene_ensg
        
        return(switch(geneEnsg,
                      Gene=nashville.plot(data1=gwas,
                              data2=metaxcan,
                              map_df=mapVersion,
                              config=labels,
                              zoom_gene=gene),
                      Ensg=nashville.plot(data1=gwas,
                                          data2=metaxcan,
                                          map_df=mapVersion,
                                          config=labels,
                                          zoom_ensg=gene),
                      stop("error")))
    }) |>
    bindCache(gwasProcess, metaProcess, input$version, labelsProcess, input$gene, input$gene_ensg)

    output$distPlot <- renderPlot({
        plotType <- input$plottype
        cat("Running ...\n")

        ptc <- proc.time()
        
        

        p <- switch(plotType,
                    `Genome Wide`=plotGenomeWide(),
                    `Zoom to Chromosome`=plotChromosome(),
                    `Zoom to Gene`=plotGene(),
                    stop("error"))
        
        print(proc.time() - ptc)
        return(p)
    }) |>
    bindCache(plotGenomeWide, input$plottype) |>
    bindEvent(input$go)
}
