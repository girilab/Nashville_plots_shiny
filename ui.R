library(shiny)

fluidPage(

    ## Application title
    titlePanel("Nashville Plot"),

    ## Sidebar with a slider input for number of bins
    verticalLayout(
        inputPanel(
            verticalLayout(
                fileInput("gwas",
                          "GWAS data file",
                          accept = "application/zip",
                          buttonLabel="Browse ..."),
                fileInput("metaxcan",
                          "Metaxcan Data Set",
                          accept="application/zip",
                          buttonLabel="Browse ..."),
                fileInput("labels",
                          "Metaxcan file -> tissue name mapping file (Optional)",
                          accept="text/plain",
                          buttonLabel="Browse ..."),
                radioButtons("version",
                             "Select Ensemble Gene Map Version",
                             list(37, 38),
                             character(0)),
                actionButton("go", "Go")
            ),
            verticalLayout(
                radioButtons("plottype",
                             "Type of plot to be made",
                             list("Genome Wide", "Zoom to Chromosome", "Zoom to Gene"),
                             selected="Genome Wide"),
                numericInput("chr",
                             "Zoom to Chromosome (1 - 25)",
                             1,
                             min=1,
                             max=25,
                             step=1),
                textInput("gene",
                          "Zoom to Gene"),
                radioButtons("gene_ensg",
                             "Type of Gene name given",
                             list("Gene Name", "Ensemble Gene Name"),
                             character(0))
            )
        ),

        ## Show a plot of the generated distribution
        plotOutput("distPlot")
    )
)
