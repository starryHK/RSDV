# This is a Shiny web application. 
# http://shiny.rstudio.com/

# packages used in this app
library(shiny)
library(ggplot2)
library(gplots)
library(DESeq2)
library(RColorBrewer)
library(shinythemes)
library(DT)
#source("DT")
source("mydds.R")
source("cmcdistance.R")


# Define UI for application 
ui <- fluidPage(
  
  theme = shinytheme("cerulean"),
  # Application title
  titlePanel("Shiny-RDAV"),
  helpText("A Web Server for  RNA-Seq Data Analysis and Visualization"),
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    sidebarPanel(
      radioButtons("data_file_type",h3("Use example data or upload your own data"),c("Example Data"="examplecounts","Upload Data"="upload"),selected = "examplecounts"),
      
      # Input: Select a file ----
      conditionalPanel(condition = "input.data_file_type=='upload'",
                       radioButtons("file_type",h3("please choose data type"),c("MATRIX"="matrix",'CSV'="csv"),selected = "matrix"),
                       
                       # Input: Select a file ----
                       conditionalPanel(condition = "input.file_type=='matrix'",
                                        fileInput("file1", "Input matrix Data",
                                                  multiple = TRUE)
                       ),
                       conditionalPanel(condition = "input.file_type=='csv'",
                                        fileInput("file3", "Input csv Data",
                                                  multiple = TRUE,
                                                  accept = c("text/csv","text/comma-separated-values,text/plain",".csv")
                                        )
                       )
      ),
      tags$hr(),
      h3("DEG Analysis"),
      # Input: Select a Group ----
      #choose condition and replicates
      selectInput("Group_factor","Please select experimental conditions",list("1 Control and 1 Experimental"="type1","1 Control and 2 expeimetnal"="type2"),selected = "type1"),
      selectInput("Group_replicate","Please select experimental replicates ",list("Two replicates for each condition "="type1","Three replicates for each condition"="type2"),selected = "type2"),
      
      sliderInput("slider1", "FDR",min = 0.0000000001, max = 0.05, value = 0.05),
      sliderInput("slider2", "log2Foldchange",min = 1.0, max = 5.0, value = 2,step = 0.5),
      
      # Horizontal line ----
      tags$hr(),
      h3("DEG Visualization"),
      h4("Boxplot Figure"), 
      h4("Volcano Figure"), 
      h4("Heatmap Figure"),
      #Z-score Choice
      selectInput("score", label=("Z-score Choice"),list("by matrix","by column")),
      #Distance Choice
      selectInput("distance", label=h5(strong("Distance Choice")),list("Euclidean","Pearson correlation distance","Manhattan")),
      #Heatmap title
      textInput("text", "Figure Title",value = "DEG"),
      # dispaly gene name
      checkboxInput("checkbox1", "Show Gene Name", value = FALSE),
      #dispaly gene cluster
      checkboxInput("checkbox2", "Show Gene Cluster", value = FALSE),
      
      # Horizontal line ----
      tags$hr(),
      h4("PCA Figure"),
      #Heatmap title
      textInput("text2", "Figure Title", 
                value = "PAC"),
      # dispaly PCA legend
      checkboxInput("checkbox4", "Show Legend", value = T),
      
      tags$hr(),
      h3("Download Tables and Figures"),
      #  download DEG Table
      h4("Download Tables"),
      #Table Format Choice
      radioButtons("checkGroup2", "table Format Choice",
                   choices = list("csv" = 1, "txt" = 2),
                   selected = 1),
      downloadButton("downloadCsv", "Download DEG Table"),
      h4("Download Figures"),
      #Figure Format Choice
      radioButtons("checkGroup", "Figure Format Choice",
                   choices = list("png" = 1, "jpeg" = 2,"pdf"=3),
                   selected = 1),
      
      #Download Heatmap Figure
      downloadButton("downloadFigure", "Download Heatmap Figure"),
      
      #Download PCA Figure
      downloadButton("downloadFigure1", "Download PCA Figure"),
      
      #Download Boxplot Figure
      downloadButton("downloadFigure2", "Download Boxplot Figure"),
      
      #Download volcanoplot Figure
      downloadButton("downloadFigure3", "Download Volcano Plot Figure")
    ),
    # Show a plot of the generated distribution
    mainPanel(
      # Output: Tabset w/ plot, summary, and table ----
      tabsetPanel(type = "tabs",
                  tabPanel("lnstruction",
                           br(),
                           em(strong(h4("The RDAV app allows users to visualize differentially expressed genes (DEG) starting with count data"))),
                           h4(em("Explore")," the app's features with the example data set pre-loaded by clicking on the tabs above."),
                           #br(),
                           strong(h4("Features")),
                           img(src = "example.png", height = 300, width = 250),
                           img(src = "example2.png", height = 300, width = 250),
                           #img(src = "example5.jpeg", height = 300, width = 225),
                           img(src = "PCA1.png", height = 270, width = 240),
                           img(src = "volcanoplot.png", height = 300, width = 300),
                           img(src = "Boxplot.png", height = 200, width = 400),
                           ###
                           #em(h4("Visualize your data:")),
                           #h5("clustering (PCA plots, heatmaps)."),
                           #h5("gene-level boxplots of expression values."),
                           #h5("gene-level volcanoplot of expression values."),
                           ###
                           br(),
                           strong(h4("Data Format")),
                           h5("Data must be uploaded as a matrix or CSV file"),
                           h5("File must be the raw counts,not normalized data,e.g.FPKM,TPKM,TPM"),
                           h5("File must have a header row."),
                           h5("First/Left-hand column(s) must be gene identifiers."),
                           h5("First/Left-hand column(s) must be defined as row_name."),
                           br(),
                           em(h4("Example of Data format")),
                           h5("Each row denotes a gene, each column denotes a sample."),
                           img(src = "example3.png", height = 210, width = 600),
                           
                           ####specify each column
                           em(h4("DEG Table")),
                           h5("Column A provide gene name."),
                           h5("Column C and column G provide Fold Changes and FDR,respectively."),
                           h5(" We use both log2FC and FDR to filter DEG data."),
                           img(src = "example4.png", height = 210, width = 600)
                  ),
                  
                  
                  tabPanel("InputData",
                           #tags$hr(),
                           br(),
                           dataTableOutput('countdataDT')),
                  tabPanel("DEG",
                           #tags$hr(),
                           br(),
                           verbatimTextOutput("analysis1"),
                           dataTableOutput('countdataDT2'),
                           dataTableOutput('countdataDT3')),
                  tabPanel("Boxplot", 
                           #tags$hr(),
                           br(),
                           
                           plotOutput("boxplot")),
                  tabPanel("Volcano plot", 
                           
                           #tags$hr(),
                           br(),
                           plotOutput("plot1", height = 600),
                           plotOutput("scatterplot", height = 600,
                                      dblclick = "scatterplot_dblclick",
                                      brush = brushOpts(
                                        id = "scatterplot_brush",
                                        resetOnNew = TRUE
                                      )
                           )
                           
                  ),
                  tabPanel("Heatmap", 
                           #tags$hr(),
                           br(),
                           br(),
                           plotOutput("plot")),
                  
                  tabPanel("PCA", 
                           #tags$hr(),
                           br(),
                           plotOutput("plot2")),
                  
                  tabPanel("Help", 
                           #tags$hr(),
                           br(),
                           # PDF of instructions link
                           downloadLink("instructionspdf",label=h4("Download Instructions (pdf)")),
                           em(strong(h4("The RDAV app allows users to visualize differentially expressed genes (DEG) starting with count data"))),
                           h4(em("Explore")," the app's features with the example data set pre-loaded by clicking on the tabs above."),
                           br(),
                           strong(h3("Instruction")),
                           h5("The app is hosted on the website:https://RDAV/
                              Code can be found on github:https://github.com/starryHK/RDAV/Instruction
                              To run this app locally on your machine,download R or Rstudio and run the following
                              commands once to set up the environment:"
                           ),
                           img(src = "code1.png", height = 110, width = 600),
                           h5("You may now run the shiny app with just one command in R:"),
                           img(src = "code2.png", height = 57, width = 600),
                           br(),
                           strong(h3("Input Data")),
                           h5("You may use this app by",
                              "1.Exploring the pre-loaded example data set.This is a pre-loaded mouse macrophages 
                              RNA-seq example for exploring the app's features.",
                              "2.Upload your own data that is Count data (or log2-expression data) which come from
                              transcriptome sequencing."),
                           
                           strong(h4("Data Format")),
                           h5("Data must be uploaded as a matrix or CSV file"),
                           h5("File must be the raw counts,not normalized data,e.g.FPKM,TPKM,TPM"),
                           h5("File must have a header row."),
                           h5("First/Left-hand column(s) must be gene identifiers."),
                           h5("First/Left-hand column(s) must be defined as row_name."),
                           br(),
                           em(h5("Example of Data format")),
                           h5("Each row denotes a gene, each column denotes a sample."),
                           img(src = "example3.png", height = 210, width = 600),
                           p("Analysis: When raw counts are uploaded, the data is then analyzed by the app. The app uses
                             the voom method from the ‘limma’ Bioconductor package to transform the raw counts into
                             logged and normalized intensity values. These values are then analyzed via linear regression
                             where gene intensity is regressed on the group factor. P-values from all pairwise regression
                             tests for group effect are computed and Benjamini-Hochberg false discovery rate adjusted pvalues are computed for each pairwise comparison"),
                           br(),
                           ####specify each column
                           em(h5("DEG Table")),
                           h5("Column A provide gene name."),
                           h5("Column C and column G provide Fold Changes and FDR,respectively."),
                           h5(" We use both log2FC and FDR to filter DEG data."),
                           img(src = "example4.png", height = 210, width = 600),
                           p("Analyzed data must contain some kind of expression measure for each sample (i.e. counts,
                             normalized intensities, CPMs), and a set of p-values with corresponding fold changes for
                             those p-values. For instance, if you have a p-value for the comparison of control vs exp , you
                             can upload the observed fold change or log2(fold change) between control vs exp. If you
                             have a more complex design and do not have fold changes readily available, you may upload
                             the test statistics or other similar measures of effect size as placeholders. The fold changes
                             are mainly used in the volcano plots. We recommend uploading p-values that are adjusted
                             for multiple comparisons (such as q-values from the qvalue package, or adjusted p-values
                             from p.adjust() function in R)."),
                           
                           br(),
                           strong(h3("Visualization")),
                           img(src = "example.png", height = 300, width = 225),
                           img(src = "example2.png", height = 300, width = 225),
                           #img(src = "example5.jpeg", height = 300, width = 225),
                           img(src = "PCA1.png", height = 250, width = 225),
                           img(src = "volcanoplot.png", height = 300, width = 300),
                           img(src = "Boxplot.png", height = 200, width = 400),
                           br(),
                           em(h5("Visualize your data:")),
                           h5("clustering (PCA plots, heatmaps)."),
                           h5("gene-level boxplots of expression values."),
                           h5("gene-level volcanoplot of expression values."),
                           br(),
                           strong(h4("PCA plots")),
                           h5("This plot uses Principal Component Analysis (PCA) to calculate the principal components of
                              the expression data using data from all genes. Euclidean distances between expression
                              values are used. Samples are projected on the first two principal components (PCs) and the
                              percent variance explained by those PCs are displayed along the x and y axes. Ideally your
                              samples will cluster by group identifier"),
                           img(src = "PCA1.png", height = 250, width = 225),
                           br(),
                           strong(h4("volcanoplots")),
                           h5("This is a scatter plot log fold changes vs –log10(p-values) so that genes with the largest fold
                              changes and smallest p-values are shown on the extreme top left and top right of the plot.
                              Hover over points to see which gene is represented by each point."),
                           
                           img(src = "volcanoplot.png", height = 300, width = 300),
                           br(),
                           strong(h5("Gene Expression Boxplots")),
                           h5("Use the search bar to look up genes in your data set. For selected gene(s) the stripchart
                              (dotplot) and boxplots of the expression values are presented for each group. You may plot
                              one or multiple genes along side each other. Hover over points for more information about
                              the data."),
                           img(src = "Boxplot.png", height = 200, width = 400),
                           
                           br(),
                           strong(h5("Heatmaps")),
                           h5("A heatmap of expression values are shown, with genes and samples arranged by
                              unsupervised clustering. You may filter on test results as well as P-value cutoffs. By default the top genes (with lowest P-values) are shown."),
                           
                           
                           img(src = "example.png", height = 300, width = 225),
                           img(src = "example2.png", height = 300, width = 225)
                           #img(src = "example5.jpeg", height = 300, width = 225)
                           
                           )
                          )
                        )
                      )
                  )



# Define server logic required to draw a histogram
server <- function(input, output) {
  #Use Example file or upload your own data
  dataInput<-reactive({
    print("inputting data")
    validate(
      need((input$data_file_type=="examplecounts")|((!is.null(input$file1))|(!is.null(input$file3))),
           message = "Please select a file")
    )
    if(input$data_file_type=="examplecounts"){
      inFile<-read.delim("EXAMPLE.matrix",
                         header = TRUE,
                         sep = "\t",
                         quote = "\t",dec = ".",
                         fill = TRUE,
                         stringsAsFactors = F,
                         row.names = 1)
    }else {
      if(input$file_type=="matrix"){
        req(input$file1)
        inFile<-read.delim(input$file1$datapath,
                           header = TRUE,
                           sep = "\t",
                           quote = "\t",dec = ".",### when upload matrix data saved by ourself ,but raw matrix,be careful of dec
                           fill = TRUE,
                           stringsAsFactors = F,
                           row.names = 1)
      }else{      
        req(input$file3)
        inFile <- read.csv(input$file3$datapath,
                           header = TRUE,
                           sep = ",",row.names = 1)
      }
    }
    return(inFile)
  })
  #ensure condition by group choice
  datasetInputcondition <- reactive({
    inFile <- dataInput()
    if (is.null(inFile))
      return(NULL)
    if (dim(inFile)[1] != 1) {
      if (input$Group_factor=="type1"){
        if (input$Group_replicate=="type1"){
          condition=c("control","control","exp","exp")
        }else {
          condition=c("control","control","control","exp","exp","exp")
        }
      }else if (input$Group_factor=="type2"){
        if (input$Group_replicate=="type1"){
          condition=c("control","control","exp","exp","exp1","exp1")
        }else{
          condition=c("control","control","control","exp","exp","exp","exp1","exp1","exp1")
        }
      }
      return(condition)
    }
  })
  
  # screen gene and get DESeq result by our defined function (mydds)
  datasetInput <- reactive({
    inFile <- dataInput()
    condition<-datasetInputcondition()
    mydds(inFile,condition)
  })
  
  #screen DEG data by self-define FDR and log2FC,and sort data
  datasetInput3 <- reactive({
    # self-define FDR and log2FC
    selfFDR <- input$slider1
    selflog2FC <- input$slider2
    resSig <- subset(datasetInput(), (datasetInput()$padj < selfFDR & abs(datasetInput()$log2FoldChange) >= selflog2FC))
    #sort the Data
    resSig<-resSig[order(resSig$log2FoldChange,decreasing = FALSE),]
    return(resSig)
  })
  
  
  #Display chosen data(upload,previous or example data) by tables
  output$countdataDT <- renderDataTable({
    alldata <-dataInput()
    if(!is.null(alldata)){
      if(input$data_file_type=="examplecounts") {
        names <- c("Control1","Control2","Control3","Exp1","Exp2","Exp3")
        colnames(alldata) <- names
      }
    }
    alldata<-as.data.frame(alldata,keep.rownames=TRUE)
  })
  
  #Display whole DEG data by table
  output$countdataDT2 <- renderDataTable({
    if (is.null(datasetInput())){
      return(NULL)
    }
    #head(myres)
    resSig<- datasetInput3()
    resSig<-as.data.frame(resSig,keep.rownames=TRUE)
  })
  #Boxplot
  output$boxplot <- renderPlot({
    alldata <- dataInput()
    if(input$data_file_type=="examplecounts") {
      names <- c("Control1","Control2","Control3","Exp1","Exp2","Exp3")
      colnames(alldata) <- names
    }
    #boxplot(alldata,col=rainbow(9),ylim=c(0,150))
    boxplot(log10(alldata+1),col=rainbow(9),pch=20, main="Boxplot", cex=1.0, xlab="Group", ylab="log[10](Counts+1)")
  })
  
  # Single zoomable plot (on left)
  ranges <- reactiveValues(x = NULL, y = NULL)
  output$scatterplot <- renderPlot({
    myres<-datasetInput()
    par(mar=c(5,5,5,5), cex=1.0, cex.main=1.4, cex.axis=1.4, cex.lab=1.4)
    degTotal<-myres
    topT <- as.data.frame(degTotal)
    # self-define FDR and log2FC
    selfFDR <- input$slider1
    selflog2FC <- input$slider2
    #plot(log2FoldChange, -log10(padj), pch=20, main="Volcano plot", cex=1.0, xlab=bquote(~Log[2]~Fold~Change), ylab=bquote(~-log[10]~FDR),xlim=c(-13,13))
    ggplot(topT, aes(log2FoldChange, -log10(padj))) +
      geom_point() +xlim(-13,13)+ylim(-5,100)+geom_vline(xintercept=c(-selflog2FC,0,selflog2FC), linetype="dotted")+ geom_hline(aes(yintercept=-log10(max(topT$pvalue[topT$padj<selfFDR], na.rm=TRUE))),linetype="dashed")+
      coord_cartesian(xlim = ranges$x, ylim = ranges$y, expand = FALSE)
  })
  
  # When a double-click happens, check if there's a brush on the plot.
  # If so, zoom to the brush bounds; if not, reset the zoom.
  observeEvent(input$scatterplot_dblclick, {
    brush <- input$scatterplot_brush
    if (!is.null(brush)) {
      ranges$x <- c(brush$xmin, brush$xmax)
      ranges$y <- c(brush$ymin, brush$ymax)
    }else {
      ranges$x <- NULL
      ranges$y <- NULL
    }
  })
  #volcano plot
  output$plot1 <- renderPlot({
    myres<-datasetInput()
    par(mar=c(5,5,5,5), cex=1.0, cex.main=1.4, cex.axis=1.4, cex.lab=1.4)
    degTotal<-myres
    topT <- as.data.frame(degTotal)
    # self-define FDR and log2FC
    selfFDR <- input$slider1
    selflog2FC <- input$slider2
    #Adjusted P values (FDR Q values)
    with(topT, plot(log2FoldChange, -log10(padj), pch=20, main="Volcano plot", cex=1.0, xlab=bquote(~Log[2]~Fold~Change), ylab=bquote(~-log[10]~FDR),xlim=c(-11,11)))
    with(subset(topT, padj<=selfFDR & log2FoldChange>=selflog2FC), points(log2FoldChange, -log10(padj), pch=20, col="red", cex=0.5))
    with(subset(topT, padj<=selfFDR & log2FoldChange<=(-selflog2FC)), points(log2FoldChange, -log10(padj), pch=20, col="green", cex=0.5))
    #with(subset(topT, padj<0.05 & abs(log2FoldChange)>2), text(log2FoldChange, -log10(padj), labels=subset(rownames(topT), topT$padj<0.05 & abs(topT$log2FoldChange)>2), cex=0.8, pos=3))
    
    #Add lines for absolute FC>2 and P-value cut-off at FDR Q<0.05
    abline(v=0, col="black", lty=3, lwd=1.0)
    abline(v=-selflog2FC, col="black", lty=4, lwd=2.0)
    abline(v=selflog2FC, col="black", lty=4, lwd=2.0)
    abline(h=-log10(max(topT$pvalue[topT$padj<selfFDR], na.rm=TRUE)), col="black", lty=4, lwd=2.0)
    legend("topright", legend=c("Up","Down","Normal"),title = "Significance",col=c("red","green","black"), pch=16, xpd=T, cex=0.5,horiz=F)
    
  })
  #heatmap
  output$plot <- renderPlot({
    
    resSig<-datasetInput3()
    deg<-row.names(resSig)
    sig <- matrix (0,nc=dim(dataInput())[2],nr=length(deg))
    
    #names <- c("NG1","NG2","NG3","RPM1","RPM2","RPM3")
    # new matrix colnames,which match raw data 
    if (input$Group_factor=="type1"){
      if (input$Group_replicate=="type1"){
        names <-c("Control1","Control2","Exp1","Exp2")
      }else {
        names <-c("Control1","Control2","Control3","Exp1","Exp2","Exp3")
      }
    }else if (input$Group_factor=="type2"){
      if (input$Group_replicate=="type1"){
        names <-c("Control1","Control2","Exp1-1","Exp1-2","Exp2-1","Exp2-2")
      }else{
        names <-c("Control1","Control2","Control3","Exp1-1","Exp1-2","Exp1-3","Exp2-1","Exp2-2","Exp2-3")
      }
    }
    sig <- as.data.frame(sig)
    rownames(sig) <- deg
    colnames(sig) <- names
    for(i in 1:length(deg)){
      sig[i,] <- (dataInput())[which(row.names(dataInput()) == deg[i] ), ]
    }
    # log transform then normalize
    sig_matrix <- data.matrix(sig)
    sig_matrix <- log10(sig_matrix+1)
    file <- sig_matrix
    #normalized data either by column or by matrix
    if (input$score == "by column"){
      activity.mean <- apply(file, 2,mean, na.rm=T)
      activity.sd <- apply(file, 2,sd,na.rm=T)
      zscore.mat <- sweep(file, 2, activity.mean, "-")
      zscore.mat <- sweep(zscore.mat, 2, activity.sd, "/")
    } else if (input$score == "by matrix") {
      activity.mean <- mean(file, na.rm=T)
      activity.sd <- sd(file, na.rm=T)
      zscore.mat <- sweep(file, 1, activity.mean, "-")
      zscore.mat <- sweep(zscore.mat, 1, activity.sd, "/")
    }
    ###
    dist.choice <- input$distance
    cluster.dist<-cmcdistance(dist.choice,zscore.mat)
    cluster.clust <- agnes(cluster.dist, diss=T, method="average")
    ordered.mat <- zscore.mat[, cluster.clust$order]
    
    # get data and heatmap order and type
    zscore.mat <- ordered.mat
    rc <- colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(length(deg))
    
    
    if (input$checkbox1==FALSE){
      labRow=NA
    } else{
      labRow=NULL  
    }
    if (input$checkbox2==FALSE){
      Rowv=NA
    }else{
      Rowv=NULL
    }
    heatmap(ordered.mat,main=input$text,col=rc,Rowv = Rowv,RowSideColors = rc,labRow = labRow)
    if (input$checkbox2==F){
    mtext(paste("-3.30  -2.57 -1.96 -1.65  0 1.65 1.96  2.57  3.30 "),side = 2)
    }
  })
  
  #PCA Figure
  output$plot2 <- renderPlot({
    dataT <- t(dataInput())
    dataT2 <-  dataT[, colSums(dataT != 0) > 0.1]
    dataT3 <- log10(dataT2+1)
    dataPCA <- prcomp(dataT3)
    dataPCA3 <- dataPCA
    pc1var <- (summary(dataPCA)$importance[2])*100
    pc1var <-paste0(as.character(pc1var),"%")
    pc1var <-paste("PC1 explains", pc1var,"variance")
    pc2var <- (summary(dataPCA)$importance[5])*100
    pc2var <-paste0(as.character(pc2var),"%")
    pc2var <-paste("PC2 explains", pc2var,"variance")
    if (input$Group_factor=="type1"){
      if (input$Group_replicate=="type1"){
        refClass <- c(1,1,2,2)      
      }else {
        refClass <- c(1,1,1,2,2,2)
      }
    }else if (input$Group_factor=="type2"){
      if (input$Group_replicate=="type1"){
        refClass <- c(1,1,2,2,3,3)      
      }else{
        refClass <- c(1,1,1,2,2,2,3,3,3)      
      }
    }
    refClass <- factor(refClass)
    par(pin=c(3.6,4)) 
    if(input$Group_factor=="type1"){
      plot(dataPCA$x[,1:2],col =refClass,cex=2,pch=16,main=input$text1,cex.main=1.5,xlab=pc1var,ylab=pc2var,cex.lab=1.3)
      if(input$checkbox4==T){
        legend(35,9, legend=c("Control","Exp"),col=c(1,2), pch=16, xpd=T, cex=0.7,horiz=F)
      }
    }else if(input$Group_factor=="type2"){
      plot(dataPCA$x[,1:2],col = refClass,cex=2,pch=16,main=input$text1,cex.main=1.5,xlab=pc1var,ylab=pc2var,cex.lab=1.3)
      if(input$checkbox4==T){
        legend(35,9, legend=c("Control","Exp1","Exp2"),col=c(1,2,3), pch=16, xpd=T, cex=0.7,horiz=F)
      }
    }
  })
  output$instructionspdf <- downloadHandler(filename="Instructions.pdf",
                                            content=function(file){
                                              file.copy("instructions/Instructions.pdf",file)
                                            })
  
  # Downloadable  csv of DEG dataset ----
  output$downloadCsv <- downloadHandler(
    filename = function() {
      if(input$checkGroup2==1){
        paste("DEG", ".csv", sep = "")
      }else{
        paste("DEG", ".txt", sep = "")
      }
    },
    content = function(file) {
      resSig<-datasetInput3()
      #write DEG data to csv
      if(input$checkGroup2==1){
        write.csv(resSig, file)
      }else {
        write.table(resSig, file) 
      }
    })
  
  # Download heatmap Figure ----
  output$downloadFigure <- downloadHandler(
    filename = function(){
      if(input$checkGroup==1){
        paste('Heatmap.png')
      }else if(input$checkGroup==2){
        paste('Heatmap.jpeg')
      }else{
        paste('Heatmap.pdf')
      }
    },
    content = function(file) { 
      if(input$checkGroup==1){
        png(file)
      }else if(input$checkGroup==2){
        jpeg(file,width=3200,height=2100,res=300)
      }else{
        pdf(file)
      }
      resSig<-datasetInput3()
      deg<-row.names(resSig)
      sig <- matrix (0,nc=dim(dataInput())[2],nr=length(deg))
      if (input$Group_factor=="type1"){
        if (input$Group_replicate=="type1"){
          names <-c("Control1","Control2","Exp1","Exp2")
        }else {
          names <-c("Control1","Control2","Control3","Exp1","Exp2","Exp3")
        }
      }else if (input$Group_factor=="type2"){
        if (input$Group_replicate=="type1"){
          names <-c("Control1","Control2","Exp1-1","Exp1-2","Exp2-1","Exp2-2")
        }else{
          names <-c("Control1","Control2","Control3","Exp1-1","Exp1-2","Exp1-3","Exp2-1","Exp2-2","Exp2-3")
        }
      }
      sig <- as.data.frame(sig)
      rownames(sig) <- deg
      colnames(sig) <- names
      for(i in 1:length(deg)){
        sig[i,] <- (dataInput())[which(row.names(dataInput()) == deg[i] ), ]
      }
      # log transform then normalize
      sig_matrix <- data.matrix(sig)
      sig_matrix <- log10(sig_matrix+1)
      file <- sig_matrix
      #normalized data either by column or by matrix
      if (input$score == "by column"){
        activity.mean <- apply(file, 2,mean, na.rm=T)
        activity.sd <- apply(file, 2,sd,na.rm=T)
        zscore.mat <- sweep(file, 2, activity.mean, "-")
        zscore.mat <- sweep(zscore.mat, 2, activity.sd, "/")
      } else if (input$score == "by matrix") {
        activity.mean <- mean(file, na.rm=T)
        activity.sd <- sd(file, na.rm=T)
        zscore.mat <- sweep(file, 1, activity.mean, "-")
        zscore.mat <- sweep(zscore.mat, 1, activity.sd, "/")
      }
      dist.choice <- input$distance
      cluster.dist<-cmcdistance(dist.choice,zscore.mat)
      cluster.clust <- agnes(cluster.dist, diss=T, method="average")
      ordered.mat <- zscore.mat[, cluster.clust$order]
      # get data and heatmap order and type
      zscore.mat <- ordered.mat
      rc <- colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(length(deg))
      
      if (input$checkbox1==FALSE)
      {
        labRow=NA
      }else{
        labRow=NULL
      }
      if (input$checkbox2==FALSE){
        Rowv=NA
      }else{
        Rowv=NULL
      }
      heatmap(ordered.mat,main="DEG",col=rc,Rowv = Rowv,RowSideColors = rc,labRow = labRow)  #labRow control genename
      if (input$checkbox2==F){
         mtext(paste("-3.30   -2.57  -1.96  -1.65  0  1.65   1.96    2.57   3.30 "), line = 1,side = 2, cex = 1,adj=-0.2)
      }
      dev.off()
    })
  
  #download PCA Figure---
  output$downloadFigure1 <- downloadHandler(
    filename = function(){
      if(input$checkGroup==1){
        paste('PCA.png')
      }else if(input$checkGroup==2){
        paste('PCA.jpeg')
      }else{
        paste('PCA.pdf')
      }
    },
    content = function(file) { 
      if(input$checkGroup==1){
        png(file)
      }else if(input$checkGroup==2){
        jpeg(file,width=3200,height=2100,res=300)
      }else{
        pdf(file)
      }
      dataT <- t(dataInput())
      dataT2 <-  dataT[, colSums(dataT != 0) > 0.1]
      dataT3 <- log10(dataT2+1)
      dataPCA <- prcomp(dataT3)
      dataPCA3 <- dataPCA
      pc1var <- (summary(dataPCA)$importance[2])*100
      pc1var <-paste0(as.character(pc1var),"%")
      pc1var <-paste("PC1 explains", pc1var,"variance")
      pc2var <- (summary(dataPCA)$importance[5])*100
      pc2var <-paste0(as.character(pc2var),"%")
      pc2var <-paste("PC2 explains", pc2var,"variance")
      
      if (input$Group_factor=="type1"){
        if (input$Group_replicate=="type1"){
          refClass <- c(1,1,2,2)      
        }else {
          refClass <- c(1,1,1,2,2,2)
        }
      }else if (input$Group_factor=="type2"){
        if (input$Group_replicate=="type1"){
          refClass <- c(1,1,2,2,3,3)      
        }else{
          refClass <- c(1,1,1,2,2,2,3,3,3)      
        }
      }
      refClass <- factor(refClass)
      par(pin=c(4.5,4.2))
      if(input$Group_factor=="type1"){
        plot(dataPCA$x[,1:2],col =refClass,cex=2,pch=16,main=input$text1,cex.main=1.5,xlab=pc1var,ylab=pc2var,cex.lab=1.3)
        if(input$checkbox4==T){
          legend(35,9, legend=c("Control","Exp"),col=c(1,2), pch=16, xpd=T, cex=0.7,horiz=F)
        }
      }else if(input$Group_factor=="type2"){
        plot(dataPCA$x[,1:2],col = refClass,cex=2,pch=16,main=input$text1,cex.main=1.5,xlab=pc1var,ylab=pc2var,cex.lab=1.3)
        if(input$checkbox4==T){
          legend(35,9, legend=c("Control","Exp1","Exp2"),col=c(1,2,3), pch=16, xpd=T, cex=0.7,horiz=F)
        }
      }
      dev.off() 
    })
  #download Boxplot Figure---
  output$downloadFigure2 <- downloadHandler(
    filename = function(){
      if(input$checkGroup==1){
        paste('Boxplot.png')
      }else if(input$checkGroup==2){
        paste('Boxplot.jpeg')
      }else{
        paste('Boxplot.pdf')
      }
    },
    content = function(file) { 
      if(input$checkGroup==1){
        png(file)
      }else if(input$checkGroup==2){
        jpeg(file,width=3200,height=2100,res=300)
      }else{
        pdf(file)
      }
      alldata <- dataInput()
      if(input$data_file_type=="examplecounts") {
        names <- c("Control1","Control2","Control3","Exp1","Exp2","Exp3")
        colnames(alldata) <- names
      }
      #boxplot(alldata,col=rainbow(9),ylim=c(0,150))
      boxplot(log10(alldata+1),col=rainbow(9),pch=20, main="Boxplot", cex=1.0, xlab="Group", ylab="log[10](Counts+1)")
      dev.off() 
    })
  
  #download Volcano Figure---
  output$downloadFigure3 <- downloadHandler(
    filename = function(){
      if(input$checkGroup==1){
        paste('Volcano.png')
      }else if(input$checkGroup==2){
        paste('Volcano.jpeg')
      }else{
        paste('Volcano.pdf')
      }
    },
    content = function(file) { 
      if(input$checkGroup==1){
        png(file)
      }else if(input$checkGroup==2){
        jpeg(file,width=3200,height=2100,res=300)
      }else{
        pdf(file)
      }
      myres<-datasetInput()
      par(mar=c(5,5,5,5), cex=1.0, cex.main=1.4, cex.axis=1.4, cex.lab=1.4)
      degTotal<-myres
      topT <- as.data.frame(degTotal)
      # self-define FDR and log2FC
      selfFDR <- input$slider1
      selflog2FC <- input$slider2
      #Adjusted P values (FDR Q values)
      with(topT, plot(log2FoldChange, -log10(padj), pch=20, main="Volcano plot", cex=1.0, xlab=bquote(~Log[2]~Fold~Change), ylab=bquote(~-log[10]~FDR),xlim=c(-11,11)))
      with(subset(topT, padj<=selfFDR & log2FoldChange>=selflog2FC), points(log2FoldChange, -log10(padj), pch=20, col="red", cex=0.5))
      with(subset(topT, padj<=selfFDR & log2FoldChange<=(-selflog2FC)), points(log2FoldChange, -log10(padj), pch=20, col="green", cex=0.5))
      #with(subset(topT, padj<0.05 & abs(log2FoldChange)>2), text(log2FoldChange, -log10(padj), labels=subset(rownames(topT), topT$padj<0.05 & abs(topT$log2FoldChange)>2), cex=0.8, pos=3))
      #Add lines for absolute FC>2 and P-value cut-off at FDR Q<0.05
      abline(v=0, col="black", lty=3, lwd=1.0)
      abline(v=-selflog2FC, col="black", lty=4, lwd=2.0)
      abline(v=selflog2FC, col="black", lty=4, lwd=2.0)
      abline(h=-log10(max(topT$pvalue[topT$padj<selfFDR], na.rm=TRUE)), col="black", lty=4, lwd=2.0)
      legend("topright", legend=c("Up","Down","Normal"),title = "Significance",col=c("red","green","black"), pch=16, xpd=T, cex=0.5,horiz=F)
      dev.off() 
    })
}

# Run the application 
shinyApp(ui = ui, server = server)

