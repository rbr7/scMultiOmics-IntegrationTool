library(BiocManager)
options(repos = BiocManager::repositories())
library(shiny)
library(shinythemes)
library(dplyr)
library(shinycssloaders)
#library(scATAC.Explorer)
#library(scRNAseq)
library(shinyjs)
library(DT)
library(rsconnect)
library(EnsDb.Hsapiens.v75)
library(Signac)
#library(SeuratData)
library(Seurat)
library(ggplot2)
library(patchwork)
library(BSgenome.Hsapiens.UCSC.hg38)
library(mailr)

######UI
ui <- fluidPage(theme=shinytheme("yeti"),
                navbarPage("scARIA: scATAC-Seq and scRNA-Seq Integration App", 
                           tabPanel("Analysis",
                                    sidebarLayout (
                                      sidebarPanel(
                                        h3(HTML("Upload datasets"), style={'text-align: center;'}),
                                        h5(HTML("Upload scATAC-Seq data"), style={'text-align: center; font-weight: bold;'}),
                                        fileInput("atacfrag", "Choose scATAC-Seq fragments.tsv File*"),
                                        fileInput("atach5", "Choose scATAC-Seq H5 File*"),
                                        fileInput("atacmeta", "Choose scATAC-Seq Metadata File*"),
                                        fileInput("atacindex", "Choose scATAC-Seq Fragments Index File*"),
                                        h5(HTML("Upload scRNA-Seq data"), style={'text-align: center; font-weight: bold;'}),
                                        fileInput("scrnamatrix", "Choose scRNA-Seq RDS/H5AD/Matrix File*"),
                                        selectInput("genome_ver", "Choose the reference genome version:",
                                                    list("hg19" = "hg19", "hg38"="hg38")
                                        ),
                                        textInput("useremail", "Input email address")
                                        actionButton('integrateupload_btn', 'Analyze and Integrate',
                                                     style="color: #fff; background-color: #4CAF50; text-align: center"),
                                        h6(HTML("*=required"), style={'text-align: center; font-weight: bold; color: red;'})
                                      ),
                                    tabsetPanel(
                                      tabPanel("scATAC-Seq", plotOutput("atacplot")),
                                      tabPanel("scRNA-Seq", plotOutput("scrnaplot")),
                                      tabPanel("Refined Clusters", plotOutput("refined_clust")),
                                      tabPanel("Peak to Gene linkages", plotOutput("peakgenelinks")),
                                      downloadLink("downloadrds", "Download Integrated RDS object")
                                    ))),
                           tabPanel("About", span(htmlOutput("about"), style="font-size: 25px;")),
                           tabPanel("GitHub", span(htmlOutput("github"), style="font-size: 25px;"))
                
                
)
)


######Server
server <- function(input, output) {
  #To upload user datasets
  options(shiny.maxRequestSize=5000*1024^2)
  

  
  ##About Section
  output$about<- renderUI({
    HTML(paste("scARIA is a Shiny app to faciliate integration of scATAC-seq and 
  scRNA-seq data from the same biological system. To use scARIA, please upload the following:",
               "1. A fragments.tsv, H5 file and a metadata file of scATAC-seq data",
               "2. An RDS file of scRNA-Seq data",
               "Click on the 'Analyze and Integrate' button to begin the integration. 
  The refined clusters and peak to gene linkages
  plots for the uploaded datasets will be displayed on the same page.",
               "The analysis pipeline is based on the packages Signac and Seurat.", sep="<br/>"))
  })
  
  ##GitHub
  output$github<- renderUI({
    tags$a(href="https://github.gatech.edu/BIOL8803-2022/team4", "https://github.gatech.edu/BIOL8803-2022/team4")
  })
  
  ##Analysis and Integration
  observeEvent(input$integrateupload_btn, {
    #Check if mandatory files have been uploaded
    
    #print(file.exists(paste0(fragment.path[[1]], 'tbi')))
    
    if(is.null(input$atacfrag) | is.null(input$atach5) | is.null(input$scrnamatrix))
    {showModal(modalDialog(
      title = "Missing File(s)!",
      paste0("You have not uploaded one or more required files. Please ensure you have uploaded the mandatory input files"),
      easyClose = TRUE,
      footer = NULL
    ))
    }
    else{ 
      ##ATAC Pipeline
      print("Welcome to scARIA")
      tempfragpath<-"/projects/team4/data"
      file.copy(input$atacfrag$datapath, paste0(tempfragpath, "/fragments.tsv.gz"), overwrite = TRUE)
      file.copy(input$atacindex$datapath, paste0(tempfragpath, "/fragments.tsv.gz.tbi"),overwrite = TRUE)
      counts <- Read10X_h5(filename = input$atach5$datapath)
      metadata <- read.csv(
        file = input$atacmeta$datapath,
        header = TRUE,
        row.names = 1
      )
      
      chrom_assay <- CreateChromatinAssay(
        counts = counts,
        sep = c(":", "-"),
        genome = input$genome_ver,
        fragments = paste0(tempfragpath, "/fragments.tsv.gz"),
        min.cells = 10,
        min.features = 200
      )
      
      pbmc <- CreateSeuratObject(
        counts = chrom_assay,
        assay = "peaks",
        meta.data = metadata
      ) 
      
      
      # extract gene annotations from EnsDb
      annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v75)
      
      # change to UCSC style since the data was mapped to hg19
      seqlevelsStyle(annotations) <- 'UCSC'
      
      # add the gene information to the object
      Annotation(pbmc) <- annotations
      
      # compute nucleosome signal score per cell
      pbmc <- NucleosomeSignal(object = pbmc)
      
      # compute TSS enrichment score per cell
      pbmc <- TSSEnrichment(object = pbmc, fast = TRUE)
      
      # add blacklist ratio and fraction of reads in peaks
      pbmc$pct_reads_in_peaks <- pbmc$peak_region_fragments / pbmc$passed_filters * 100
      pbmc$blacklist_ratio <- pbmc$blacklist_region_fragments / pbmc$peak_region_fragments
      
      #Remove cells that are outliers for these QC metrics
      pbmc <- subset(
        x = pbmc,
        subset = peak_region_fragments > 3000 &
          peak_region_fragments < 20000 &
          pct_reads_in_peaks > 15 &
          blacklist_ratio < 0.05 &
          nucleosome_signal < 4 &
          TSS.enrichment > 2
      )
      
      #Normalization and linear dimensional reduction
      pbmc <- RunTFIDF(pbmc)
      pbmc <- FindTopFeatures(pbmc, min.cutoff = 'q0')
      pbmc <- RunSVD(pbmc)
      
      #Nonlinear dimenison reduction and clustering
      pbmc <- RunUMAP(object = pbmc, reduction = 'lsi', dims = 2:30)
      pbmc <- FindNeighbors(object = pbmc, reduction = 'lsi', dims = 2:30)
      pbmc <- FindClusters(object = pbmc, verbose = FALSE, algorithm = 3)
      output$atacplot<-renderPlot({
        DimPlot(object = pbmc, label = TRUE) + NoLegend()
      })
      
      gene.activities <- GeneActivity(pbmc)
      # add the gene activity matrix to the Seurat object as a new assay and normalize it
      pbmc[['RNA']] <- CreateAssayObject(counts = gene.activities)
      pbmc <- NormalizeData(
        object = pbmc,
        assay = 'RNA',
        normalization.method = 'LogNormalize',
        scale.factor = median(pbmc$nCount_RNA)
      )
      
      
      DefaultAssay(pbmc) <- 'RNA'
      FeaturePlot(
        object = pbmc,
        features = c('MS4A1', 'CD3D', 'LEF1', 'NKG7', 'TREM1', 'LYZ'),
        pt.size = 0.1,
        max.cutoff = 'q95',
        ncol = 3
      )
      
      
      ##RNA-Seq
      pbmc_rna <- readRDS(input$scrnamatrix$datapath)
      
      transfer.anchors <- FindTransferAnchors(
        reference = pbmc_rna,
        query = pbmc,
        reduction = 'cca'
      )
      output$scrnaplot<-renderPlot({
        DimPlot(object = pbmc_rna, group.by = 'celltype',
                label = TRUE,
                repel = TRUE) + NoLegend() + ggtitle('scRNA-seq')
      })
      #DimPlot(object = pbmc, label = TRUE) + NoLegend()
      predicted.labels <- TransferData(
        anchorset = transfer.anchors,
        refdata = pbmc_rna$celltype,
        weight.reduction = pbmc[['lsi']],
        dims = 2:30
      )
      
      
      pbmc <- subset(pbmc, idents = 14, invert = TRUE)
      pbmc <- RenameIdents(
        object = pbmc,
        '0' = 'CD14 Mono',
        '1' = 'CD4 Memory',
        '2' = 'CD8 Effector',
        '3' = 'CD4 Naive',
        '4' = 'CD14 Mono',
        '5' = 'DN T',
        '6' = 'CD8 Naive',
        '7' = 'NK CD56Dim',
        '8' = 'pre-B',
        '9' = 'CD16 Mono',
        '10' = 'pro-B',
        '11' = 'DC',
        '12' = 'NK CD56bright',
        '13' = 'pDC'
      )
      pbmc <- AddMetaData(object = pbmc, metadata = predicted.labels)
            output$refined_clust<-renderPlot({
        DimPlot(object = pbmc, group.by = 'predicted.id',
                label = TRUE,
                repel = TRUE) + NoLegend() + ggtitle('scATAC-seq')
      })
      #Differential peak analysis
      DefaultAssay(pbmc) <- 'peaks'
      da_peaks <- FindMarkers(
        object = pbmc,
        ident.1 = "CD4 Naive",
        ident.2 = "CD14 Mono",
        test.use = 'LR',
        latent.vars = 'peak_region_fragments'
      )
      print(head(da_peaks))
     


      output$peakgenelinks<-renderPlot({
        CoveragePlot(
          object = pbmc,
          region = rownames(da_peaks)[1],
          extend.upstream = 40000,
          extend.downstream = 20000
        )
      })
      
      
    output$downloadrds <- downloadHandler(
        filename = "scARIA_Integrated_scATAC.rds",
        content = function(file) {
        saveRDS(object = pbmc, file = "pbmc_atac_integrated.rds")
        }
      )
    send.mail(from = "scaria8803@gmail.com",
          to = input$useremail,
          subject = "Your scATAC and scRNA-seq integration results",
          body = paste("Please access the integrated RDS object in the attachment."),
          smtp = list(host.name = "smtp.gmail.com", port = 465, user.name = "scaria8803", passwd = "scaria8803!", ssl = TRUE),
          authenticate = TRUE,
          send = TRUE,
          attach.files = #AttachmentPath,
          file.names = "scARIA_Integrated_RDS.rds")
    }
  })

  
  
}

# Run the application 
shinyApp(ui = ui, server = server)
