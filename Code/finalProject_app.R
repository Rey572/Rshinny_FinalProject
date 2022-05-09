#ReyWu

## Author: Rey
## BU BF591
## Final Project

#load library
library(shiny)
library(ggplot2)
library(colourpicker) 
library(dplyr)
library(glue)
library(psych) 
library(DT)
library(tibble)
library(RColorBrewer)
library(stringr) 
library(gridExtra)


ui <- fluidPage(
  #introduction
  titlePanel("BF591 Final Project -- Rey Wu"),
  p("Reference:A nineteen gene‐based risk score classifier predicts prognosis of colorectal cancer patients. PMID:25049118"),
  p("Raw data consists of 54 samples (normal colon, primary CRC, and liver metastases) from 18 individuals."),
  #inputs
  mainPanel(
    tabsetPanel(
      tabPanel("Sample Information", 
               sidebarPanel(fileInput("file1", paste("Sample information matrix in CSV format"))),
               mainPanel(tabsetPanel(
                 tabPanel("Summary Table", tableOutput("sample_summary")),
                 tabPanel("Sortable Table", dataTableOutput("sample_table")),
                 tabPanel("Histogram", plotOutput("sample_plot"))
               ))),
      tabPanel("Counts Matrix", 
               sidebarPanel(
                  fileInput("file2", paste("Normalized counts matrix in CSV format")),
                  sliderInput(inputId="slider1", min=0, max=1, label="Select the percentile of variance: ", value=0, step=0.1),
                  sliderInput(inputId="slider2", min=0, max=54, label="Select the number of non-zero samples: ", value=0, step=1),
                  submitButton("Plot", width = "100%")
                  ),
               mainPanel(tabsetPanel(
                 tabPanel("Summary Table", tableOutput("counts_summary")),
                 tabPanel("Diagnostic Scatter Plots", plotOutput("counts_scatter_D")),
                 tabPanel("Clustered Heatmap", plotOutput("counts_heatmap")),
                 tabPanel("PCA Scatter Plot", 
                          sidebarPanel(
                            fileInput("file4", paste("Metadata in CSV format")),
                            radioButtons('pca_x', "Choose x-axis principal component",
                                         choices = c("1","2","3","4")),
                            radioButtons('pca_y', "Choose y-axis principal component",
                                         choices=c("1","2","3","4")),
                            submitButton("Plot", width = "100%")
                          ),
                          mainPanel( 
                            plotOutput("counts_scatter_PCA")
                          )
                         )
               ))),
      tabPanel("Differential Expression", 
               sidebarPanel(
                 fileInput("file3", paste("Results of a differential expression analysis in CSV format")),
                 p("A volcano plot can be generated with log2 fold-change on the x-axis and p-adjusted on the y-axis."),
                 radioButtons('x_name', "Choose thw column for the x-axis",
                              choices = c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")),
                 radioButtons('y_name', "Choose thw column for the y-axis",
                              choices=c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")),
                 colourInput('color1', "Base point color", value='#22577A'),
                 colourInput('color2', "Highlight point color", value='#FFCF56'),
                 sliderInput(inputId="slider3", min=-23, max=0, label="Select the magnitude of the p adjusted: ",
                             value=0, step=1),
                 submitButton("Plot", width = "100%")
               ),
               mainPanel(
                 tabsetPanel(
                   tabPanel("Sortable Table", dataTableOutput("DE_table")),
                   tabPanel("Volcano Plot", plotOutput("DE_volcano")),
                   tabPanel("Filtered Table", tableOutput("DE_table_2"))
               ))),
      tabPanel("Gene Set Enrichment Analysis",
               sidebarPanel(
                 fileInput("file5", paste("fgsea results from the differential expression data in CSV format")),
               ),
               mainPanel(
                 tabsetPanel(
                   tabPanel("Barplot", 
                            sidebarPanel(
                              sliderInput(inputId="slider4", min=1, max=50, label="Select the number of top pathways: ",
                                          value=0, step=1),
                              submitButton("Plot", width = "100%")
                            ),
                            plotOutput("enrich_bar")
                            ),
                   tabPanel("Filtered Table", 
                            sidebarPanel(
                              sliderInput(inputId="slider5", min=0.01, max=1, label="Select the value of the adjusted p-value: ",
                                          value=0, step=0.01),
                              radioButtons('pathways_type', "Choose NES pathways type:",
                                           choices = c("all", "positive", "negative")),
                              submitButton("Plot", width = "100%")
                            ),
                            dataTableOutput("enrich_table")

                            ),
                   tabPanel("Scatter plot ",
                            sidebarPanel(
                              sliderInput(inputId="slider6", min=0.01, max=1, label="elect the value of the adjusted p-value: ",
                                          value=0, step=0.01),
                              submitButton("Plot", width = "100%")
                            ),
                            plotOutput("enrich_scatter")
                            )
                 )
               )
        
      )
      
    )
  )
)


# Define server logic required to draw a histogram
server <- function(input, output, session) {
  
#------------------Sample Tab------------------------
  #load_Data
  load_data1 <- reactive({
    data <- read.csv(input$file1$datapath, header=T)  %>%
    return(data)
  })
  
  sample_summary <- function(dataf) {
    summary_table <- describe(dataf, fast=TRUE)
    return(summary_table)
  }
  
  sample_table <- function(dataf) {
    dataf <- DT::datatable(dataf, extensions = 'Buttons', class = "display")
    return(dataf)
  }
  
  sample_plot <- function(dataf) {
    count <- colSums(Filter(is.numeric, dataf))
    bar_table <- tibble('sample' = colnames(dataf[,-c(1,2)]),
                         'count' = count)
    p <- ggplot(bar_table, aes(x=sample, y=count)) +
          geom_bar(stat="identity") +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    return(p)
  }
  
  #outputs
  output$sample_summary <- renderTable({
    sample_summary(load_data1())
  })
  output$sample_table <- DT::renderDataTable({
    sample_table(load_data1())
  }) 
  output$sample_plot <- renderPlot({
    sample_plot(load_data1())
  })
#---------------------------------------------------
  
  
  
  
#------------------Counts Tab------------------------
  load_data2 <- reactive({
    data <- read.csv(input$file2$datapath, header=T) %>%
      select(!"X.1") %>%
      dplyr::rename('gene'='X') %>%
      column_to_rownames(var="gene")
    return(data)
    
  })
  load_data4 <- reactive({
    data <- read.csv(input$file4$datapath, header=T) %>%
    return(data)
  })
  
  counts_summary <- function(dataf, slider1, slider2){
    sample_num <- ncol(dataf)
    gene_num <- nrow(dataf)
    
    dataf <- dataf %>% 
      mutate(non_zero = rowSums(dataf != 0),
             variance=apply(dataf, 1, var)) 
    passed <- dataf %>%
      filter(non_zero>=slider2)%>%
      filter(variance>=quantile(variance, slider1))
    
    passed_num <- nrow(passed)
    passed_per <- (passed_num/gene_num)*100
    
    table <- tibble(Samples = sample_num,
                    Genes = gene_num,
                    num_passed = passed_num,
                    num_passed_percentage = passed_per,
                    num_failed = gene_num-passed_num,
                    num_failed_percentage = 100-passed_per)
   
    return(table)
  }
  
  counts_scatter_D <- function(dataf, slider1, slider2){
    #---median vs. variance
    plot_data <- tibble(median=apply(dataf, 1, median), 
                        variance=apply(dataf, 1, var),
                        rank=rank(median),
                        non_zero = rowSums(dataf != 0))
    passed <- filter(plot_data, quantile(variance, slider1)<variance)
    passed_min <- min(passed$variance)
    plot_data$pass_var=TRUE
    plot_data$pass_var[plot_data$variance>passed_min] = TRUE
    plot_data$pass_var[plot_data$variance<passed_min] = FALSE
    
    plot1 <- ggplot2::ggplot(plot_data, aes(x=rank, y=variance, col=pass_var)) +
      ggplot2::geom_point() +
      ggplot2::xlab("Rank(Median)") +
      ggplot2::ylab("Variance") 
    plot1 <- plot1 + ggplot2::scale_y_log10()
    
    #---median vs. non-zero
    plot_data$pass_zero=TRUE
    plot_data$pass_zero[plot_data$non_zero>slider2] = TRUE
    plot_data$pass_zero[plot_data$non_zero<slider2] = FALSE
    plot2 <- ggplot2::ggplot(plot_data, aes(x=rank, y=non_zero, col=pass_zero)) +
      ggplot2::geom_point() +
      ggplot2::xlab("Rank(Median)") +
      ggplot2::ylab("Non_zero") +
      ggplot2::scale_x_log10(oob = scales::squish_infinite)
    return(grid.arrange(plot1, plot2, ncol=2))
  }
  
  counts_heatmap <- function(dataf, slider1, slider2){
    dataf <- dataf %>% 
      mutate(non_zero = rowSums(dataf != 0),
             variance=apply(dataf, 1, var)) 
    passed <- dataf %>%
      filter(non_zero>=slider2)%>%
      filter(variance>=quantile(variance, slider1))
    # passed <- filter(dataf, rowSums(across(where(is.numeric)))!=0) 
    # variance=apply(passed, 1, var)
    # passed <- filter(passed[,1:54], quantile(variance, slider1)<variance)
    
    passed <- as.matrix(passed[,1:54])
    col.pal <- RColorBrewer::brewer.pal(10, 'RdBu')
    heatmap <- heatmap(passed,col = colorRampPalette(brewer.pal(8,"Blues"))(3))
    
    legend(x = "right", legend = c("normal", "primary tumor", "liver metastases"),
           cex = 0.8, fill = colorRampPalette(brewer.pal(8, "Blues"))(3))
    return(heatmap)
  }
  
  counts_scatter_PCA <- function(dataf, meta_info, pca_x, pca_y){
    pca <- prcomp(t(dataf))
    plot_data <- meta_info
    plot_data$PCx <- pca$x[ , pca_x]
    plot_data$PCy <- pca$x[ , pca_y]
    percent_var <- pca$sdev^2 / sum(pca$sdev^2 )
    pca_plot <- ggplot2::ggplot(plot_data, ggplot2::aes(x=PCx, y=PCy, col=tissueType)) +
      ggplot2::geom_point() +
      ggplot2::xlab(paste0(glue("PC{pca_x}: ",round(percent_var[1] * 100),"% variance"))) +
      ggplot2::ylab(paste0(glue("PC{pca_y}: ",round(percent_var[2] * 100),"% variance"))) +
      ggplot2::theme(legend.position="top")
    
    return(pca_plot)
  }
  
  output
  output$counts_summary <- renderTable({
    counts_summary(load_data2(), as.numeric(input$slider1), as.numeric(input$slider2))
  })
  output$counts_scatter_D <- renderPlot({
    counts_scatter_D(load_data2(), as.numeric(input$slider1), as.numeric(input$slider2))
  })
  output$counts_heatmap <- renderPlot({
    counts_heatmap(load_data2(), as.numeric(input$slider1), as.numeric(input$slider2))
  })
  output$counts_scatter_PCA <- renderPlot({
    counts_scatter_PCA(load_data2(), load_data4(),as.numeric(input$pca_x), as.numeric(input$pca_y))
  })
#---------------------------------------------------
 
  
   
  
#------------------DE Tab------------------------
  load_data3 <- reactive({
    data <- read.csv(input$file3$datapath, header=T) %>%
      dplyr::rename(gene=X)
    return(data)
  })

  DE_table <- function(dataf) {
    dataf <- DT::datatable(dataf, extensions = 'Buttons', class = "display")
    return(dataf)
  }

  DE_volcano <- function(dataf, x_name, y_name, slider, color1, color2) {
      plot <- ggplot(dataf, aes(!!sym(x_name), -log10(!!sym(y_name)), colour=!!sym(y_name)<10^slider)) +
        geom_point() +
        theme_bw() +
        scale_color_manual(values=c(color1, color2)) +
        xlab(x_name) +
        guides(colour = guide_legend(title = glue("pvalue<10^{slider}"))) +
        theme(legend.position = "bottom")
      return(plot)
    }

  DE_table_2 <- function(dataf, slider) {
    table <- filter(dataf, padj< 10^slider) %>%
      apply(2, function(x) formatC(x))
    return(table)
  }
  #output
  output$DE_table <- DT::renderDataTable({
    DE_table(load_data3())
  })
  output$DE_volcano <- renderPlot({
    DE_volcano(load_data3(), input$x_name, input$y_name, as.numeric(input$slider3), input$color1, input$color2)
  })
  output$DE_table_2 <- renderTable({
    DE_table_2(load_data3(), as.numeric(input$slider3))
  })
#----------------------------------------
  
  
  
#-------------enrichment analysis-------------------
  load_data5 <- reactive({
    data <- read.csv(input$file5$datapath, header=T)
    return(data)
  })
  
  enrich_bar <- function(fgsea_results, num_paths){
    top_pos <- fgsea_results%>% filter(padj < .25 & NES > 0) %>% slice_max(NES, n=num_paths)
    top_neg <- fgsea_results%>% filter(padj < .25 & NES < 0) %>% slice_min(NES, n=num_paths)
    
    subset <- fgsea_results %>% 
      filter(pathway %in% c(top_pos$pathway, top_neg$pathway)) %>%
      mutate(pathway = factor(pathway)) %>%
      mutate(plot_name = str_replace_all(pathway, '_', ' '))
    
    barplot <- subset %>% 
      mutate(plot_name = forcats::fct_reorder(factor(plot_name), NES)) %>%
      ggplot() +
      geom_bar(aes(x=plot_name, y=NES, fill = NES > 0), stat='identity', show.legend = FALSE) +
      scale_fill_manual(values = c('TRUE' = 'red', 'FALSE' = 'blue')) + 
      theme_minimal(base_size = 8) +
      ggtitle('fgsea results for Hallmark MSigDB gene sets') +
      ylab('Normalized Enrichment Score (NES)') +
      xlab('') +
      scale_x_discrete(labels = function(x) str_wrap(x, width = 80)) +
      coord_flip()
    return(barplot)
  }
  
  enrich_table <- function(dataf, slider5, pathways_type){
    table <- filter(dataf, padj<slider5)
    if(pathways_type=="negative"){
      subset <- filter(table, NES<0)
    }
    if(pathways_type=="positive"){
      subset <- filter(table, NES>0)
    }
    if(pathways_type=="all"){
      subset <- table
    }
    subset <- DT::datatable(subset, extensions = 'Buttons', class = "display")
    return(subset)
  }
  
  enrich_scatter <- function(dataf, slider6){
    sca_plot <- ggplot2::ggplot(dataf, aes(x=NES, y=-log10(padj), color=padj<slider6)) +
      ggplot2::geom_point() +
      guides(colour = guide_legend(title = glue("padj<{slider6}"))) 
    
    return(sca_plot)
  }
  
  #output
  output$enrich_bar <- renderPlot({
    enrich_bar(load_data5(),as.numeric(input$slider4))
  })
  output$enrich_table <- DT::renderDataTable({
    enrich_table(load_data5(), as.numeric(input$slider5), input$pathways_type)
  })
  output$enrich_scatter <- renderPlot({
    enrich_scatter(load_data5(), as.numeric(input$slider6))
  })
}

# Run the application
shinyApp(ui = ui, server = server)
