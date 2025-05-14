library(shiny)
library(bslib)
library(ggplot2)
library(tidyverse)
library(forcats)
library(dplyr)
library(ReactomePA)
library(clusterProfiler)
library(colourpicker) 
library(DT) 
library(patchwork) 
library(biomaRt)
library(gplots) # For heatmap.2()
library(purrr) 
library(org.Hs.eg.db)
library(rtracklayer)
library(fgsea)

# Define UI for application:
ui <- fluidPage(
  
  # Font for the entire webpage
  tags$head(
    tags$style(
      HTML("* { font-family: 'Times New Roman', serif; }"), 
      HTML(".nav-tabs>li>a {color: #2C3539}") 
    )
  ),
  
  # Page headings
  h1("BF591 Final Project - Vassanth M"),
  h4("A simple bioinformatic analyses of the dataset - 'mRNA-Seq Expression profiling of human post-mortem BA9 brain tissue for Huntington’s Disease and neurologically normal individuals'"),
  
  # Separate main panels for plot and table
  tabsetPanel(id = "tabset",
              
              # Sample info Tab panel
              tabPanel("Samples",
                       sidebarLayout(
                         sidebarPanel(
                           tags$head(tags$style(".btn-file {background-color:#2980B9; border-color: #2980B9;}.progress-bar{color:black; background-color:#3498DB;}")),
                           fileInput("sample_info_upload", paste0("Load sample information"), accept = c(".csv"))
                         ), 
                         mainPanel(
                           tabsetPanel(
                             tabPanel("Summary", DTOutput("sample_info_summary_table")),
                             tabPanel("Table", DTOutput("sample_data")),
                             tabPanel("Plots", plotOutput("sample_hist"))
                           )
                         )
                       )
              ), 
              
              # Counts info Tab panel
              tabPanel("Counts",
                       sidebarLayout(
                         sidebarPanel(
                           tags$head(tags$style(".btn-file {background-color:#3498DB; border-color: #3498DB;}.progress-bar{color:black; background-color:#2ECC71;}")),
                           fileInput("count_info_upload", paste0("Load counts data file"), accept = c(".csv", ".tsv")),
                           tags$style(HTML(".js-irs-0 .irs-single, .js-irs-0 .irs-bar-edge, .js-irs-0 .irs-bar {background: #66cc99; color: black}")),
                           sliderInput("var_range", "Select the minimum percentile of variance:",
                                       min = 0, max = 100, value = 50),
                           sliderInput("count_range", "Select the non-zero sample threshold for genes:",
                                       min = 0, max = 69, value = 0),
                           h5("Genes Passing Filters"),
                           colourInput("count_passed_filter_color", NULL, "#C0C7CC"), # base point color
                           h5("Genes Not Passing Filters"),
                           colourInput("count_failed_filter_color", NULL, "#F22F2C"), # highlight point color
                           actionButton("do", "Filter", width = "100%", icon = icon("filter"), 
                                        style = "color: black; background-color: #66cc99; border-color: #66cc99")
                         ), 
                         mainPanel(
                           tabsetPanel(
                             tabPanel("Summary", DTOutput("count_info_summary_table")),
                             tabPanel("Scatter Plots", plotOutput("count_summary_scatters")),
                             tabPanel("Heatmap", plotOutput("count_heat")),
                             tabPanel("PCA", 
                                      fluidRow(
                                        sidebarPanel(
                                          radioButtons("pca_x_axis", "Pick the PCA to plot on the x-axis",
                                                       c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10"), selected = "PC1"),
                                          radioButtons("pca_y_axis", "Pick the PCA to plot on the y-axis",
                                                       c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10"), selected = "PC2"),
                                          actionButton("pca_plot_do", "Plot", width = "100%", icon = icon("brush"), 
                                                       style = "color: black; background-color: #6FCFD1; border-color: #6FCFD1")
                                        ),
                                        column(8, plotOutput("count_pca"))
                                      )
                             )
                           )
                         )
                       )
              ),
              
              # Differential Expression info Tab panel
              tabPanel("DE",
                       sidebarLayout(
                         sidebarPanel(
                           tags$head(tags$style(".btn-file {background-color:#e86866; border-color: #e86866;}.progress-bar{color:black; background-color:#66cc99;}")),
                           fileInput("deseq_upload", paste0("Load differential expression results"), accept = c(".csv", ".tsv")),
                           tags$style(HTML(".js-irs-0 .irs-single, .js-irs-0 .irs-bar-edge, .js-irs-0 .irs-bar {background: #66cc99; color: black}")),
                           sliderInput("padj_range", "Select the Adjusted P-Value Filtering Threshold:",
                                       min = 0, max = 1, step = 0.01, value = 0.1),
                           actionButton("deseq_do", "Filter", width = "100%", icon = icon("filter"), 
                                        style = "color: black; background-color: #66cc99; border-color: #66cc99")
                         ),  
                         mainPanel(
                           tabsetPanel(
                             tabPanel("DE Results", DTOutput("deseq_tab")),
                             tabPanel("Plots",
                                      tabsetPanel(
                                        tabPanel("Raw P-values", plotOutput("pval_hist")),
                                        tabPanel("Histogram of Log2 Fold Changes", plotOutput("logfc_hist")),
                                        tabPanel("Volcano Plot of Differential Expression Results",
                                                 fluidRow(
                                                   sidebarPanel(
                                                     h5("Base point color"),
                                                     colourInput("volc_col", NULL, "#EDE60C"), # base point color
                                                     h5("Highlight point color"),
                                                     colourInput("volc_col2", NULL, "#15AD0A"), # highlight point color
                                                     sliderInput("volc_plot_range", "Select -log10(adjusted p-value) to filter DESeq2 results:",
                                                                 min = 0, max = 35, value = 15),
                                                     actionButton("volc_do", "Filter", width = "100%", icon = icon("filter"), 
                                                                  style = "color: black; background-color: #66cc99; border-color: #66cc99")
                                                   ), 
                                                   column(8, plotOutput("volc_plot"))
                                                 )
                                        )
                                      )
                             )
                           )
                         )
                       )
              ),
              
              # Reactome PA info Tab panel
              tabPanel("Reactome Pathway Enrichment Analysis",
                       sidebarLayout(
                         sidebarPanel(
                           tags$head(tags$style(".btn-file {background-color:#e86866; border-color: #e86866;}.progress-bar{color:black; background-color:#66cc99;}")),
                           fileInput("reactome_upload", paste0("Load DE results as csv"), accept = c(".csv"))
                         ),
                         mainPanel(
                           tabsetPanel(
                             tabPanel("Reactome PA", actionButton("run_analysis", "Run Analysis", class = "btn btn-primary")),
                             tabPanel("Pathway Table", DTOutput("pathway_table")),
                             tabPanel("Barplot", plotOutput("barplot")),
                             tabPanel("Dotplot", plotOutput("dotplot")),
                             tabPanel("Network Plot", plotOutput("cnetplot"))
                           )
                         )
                       )
              )
  ) # Close tabsetPanel
) # Close fluidPage

  

# Define server logic
server <- function(input, output) 
  options(shiny.maxRequestSize = 30 * 1024^2) #set file size limit high so you don't exceed limit
  
  ##SAMPLE INFO
    # sample_info_data
    #  loads sample information dataframe
  
  sample_info_data <- reactive({   
    req(input$sample_info_upload)
    df <- read.csv(input$sample_info_upload$datapath, stringsAsFactors = FALSE)    
    return(df)
    
  })
  
  
    # sample_info_summary
    #  takes sample information data, creates special summary table with means of continuous variables
    #  sum_data (df): a data frame with sample information
  
  
  sample_info_summary <- function(sum_data) {
    
    # Retrieve and process sample info data
    summary_df <- sample_info_data() %>%
      mutate(Condition = factor(Condition, levels = c("HD", "Control")))
    
    # Create a new tibble summarizing column info
    summary_tib <- tibble(
      "Column Name" = names(summary_df),
      "Type" = sapply(summary_df, class),
      "Mean (sd) or Distinct Values" = lapply(seq_along(summary_df), function(i) {
        x <- summary_df[[i]]  # Extract column data
        col_name <- names(summary_df)[i]  # Get column name
        
        if (is.numeric(x) && !grepl("^mRNA_Seq_reads_", col_name)) {
          # Calculate mean and standard deviation
          mean_val <- round(mean(x, na.rm = TRUE), 2)
          sd_val <- round(sd(x, na.rm = TRUE), 2)
          paste0(mean_val, " (+/- ", sd_val, ")")
        } else {
          # Concatenate distinct values for categorical or non-matching data
          paste(unique(x), collapse = ", ")
        }
      })
    )
    
    # Return the summary as a datatable
    datatable(summary_tib, rownames = FALSE, options = list(pageLength = nrow(summary_df)))
  }
  
  
  #sample info summary tab
  output$sample_info_summary_table <- renderDT({
    sample_info_summary()
  })
  
    # sample_data_table
    #  takes sample information data, displays it as a sortable table
    #  sum_data (df): a data frame with sample information
  
  sample_data_table <- function(sum_data){  
    sample_data_df <- sample_info_data()
    datatable(sample_data_df, rownames = FALSE, options = list(pageLength = nrow(sample_data_df)))
  }
  
  #sample data table
  output$sample_data <- renderDT({
    sample_data_table()
  })
  
    # continuous_sample_histograms
    #  creates histogram plots of continuous variables in the input sample data frame
    #  sum_data (df): a data frame with sample information
  continuous_sample_histograms <- function(sum_data){
    sample_data_df <- sample_info_data()
    
    # Extract continuous data and variables
    continuous_data_df <- sample_data_df %>%
      select_if(is.numeric)%>%
      gather()
    
    ggplot(gather(continuous_data_df), aes(value, fill = key))+
      geom_histogram(bins = 10)+
      theme_bw()+
      facet_wrap(~key, scales = "free_x")
  }
  
  #sample summary histograms
  output$sample_hist <- renderPlot({
    continuous_sample_histograms()
  })
  
  ##COUNTS INFORMATION
    # counts_data
    #  loads RNA counts data frame
  counts_data <- reactive({
    req(input$count_info_upload)
    count_df <- read.csv(input$count_info_upload$datapath, stringsAsFactors = FALSE, row.names = 1)
    return(count_df)
  })
  
    # counts_summary_table
    #  creates table that summarizes effects of counts filtering, including:
    #number of samples, total number of genes, number and % of genes passing current filter,
    #number and % of genes not passing current filter
    #  count_info_data(df): a data frame with counts information
    #  var_filter (int): a variance percentile value for which to filter the counts
    #  count_filter (int): a number of zero counts per gene to use for filtering
  
  counts_summary_table <- function(count_info_data, var_filter, count_filter){
    #filtering step
    #calculate row/gene variances
    gene_variances <- apply(count_info_data, 1, var)
    
    #convert slider value to probability value between 0 and 1
    var_filter_prob <- var_filter/100 
    
    #calculate variance threshold
    var_threshold <- quantile(gene_variances, var_filter_prob)
    filtered_counts <- count_info_data %>%
  
      filter(
        rowSums(. != 0) >= count_filter &
          apply(., 1, function(row) var(row) >= var_threshold)      
      )
    
    #create new tibble with all the info we want
    count_summary_tib <- tibble(
      "Number of Samples" = ncol(count_info_data),
      "Total Number of Genes" = nrow(count_info_data),
      "Number and % of Genes Passing Filters" = paste0(nrow(filtered_counts), " ( ", round((nrow(filtered_counts)/nrow(count_info_data))*100, 2), "%)"),
      "Number and % of Genes Not Passing Filters" =  paste0(nrow(count_info_data) - nrow(filtered_counts), " ( ", round(((nrow(count_info_data) - nrow(filtered_counts))/nrow(count_info_data))*100, 2), "%)")
    )
    
    datatable(count_summary_tib, rownames = FALSE)
  }
  
  
  #counts info summary tab
  output$count_info_summary_table <- renderDT({
    input$do #For slider input
    counts_summary_table(count_info_data = counts_data(), var_filter = isolate(input$var_range), count_filter = isolate(input$count_range))
    
  })
  
    # counts_diagnostic_scatterplots
    #  creates two scatter plots that show genes ranked by median count vs. log10(variance) and
    # genes ranked by median count vs. number of zeros for each gene that passes the filters. Genes
    # Passing filters are in darker colors, filtered genes in lighter color
    #  count_info_data(df): a data frame with counts information
    #  var_filter (int): a variance percentile value for which to filter the counts
    #  count_filter (int): a number of zero counts per gene to use for filtering
    #  color_1: a color to indicate which counts passed the filters
    #'  color_2: a color to indicate which counts did not the filters
  
  counts_diagnostic_scatterplots <- function(count_info_data, var_filter, count_filter, color_1, color_2){
    
    #creating a tibble storing the gene medians, variances, the sum of zero 
    #counts for each gene, and rankings  
    #filtering steps
    
    #calculate row/gene variances
    gene_variances <- apply(count_info_data, 1, var)
    
    #convert slider value to probability value between 0 and 1
    var_filter_prob <- var_filter/100 
    
    #calculate variance threshold
    var_threshold <- quantile(gene_variances, var_filter_prob)
    
    
    result_tib <- count_info_data %>%
      tibble("gene_count_medians" = apply(dplyr::select(.,-1), 1, median),
             "gene_variance" = apply(dplyr::select(.,-1), 1, var),
             "sum_zero_counts" = apply(dplyr::select(.,-1) == 0, 1, sum),
             "ranked_medians" = rank(gene_count_medians)
      )
    
    
    med_var_scatter <- ggplot(result_tib, aes(x = ranked_medians, y = log10(gene_variance))) +
      geom_point(aes(color = ifelse(rowSums(count_info_data != 0) >= count_filter &
                                      apply(count_info_data, 1, function(row) var(row) >= var_threshold),
                                    "Passed", "Filtered Out")), position = "jitter") +
      scale_color_manual(values = c("Filtered Out" = color_2, "Passed" = color_1)) +
      theme_bw() +
      labs(title = "Gene Variance vs. Ranked Median Counts",
           color = "Genes Passing Variance and Counts Filters")
    
    med_zeros_scatter <- ggplot(result_tib, aes(x = ranked_medians, y = sum_zero_counts)) +
      geom_point(aes(color = ifelse(rowSums(count_info_data != 0) >= count_filter &
                                      apply(count_info_data, 1, function(row) var(row) >= var_threshold),
                                    "Passed", "Filtered Out")), position = "jitter") +
      scale_color_manual(values = c("Filtered Out" = color_2, "Passed" = color_1)) +
      theme_bw() +
      labs(title = "Sum of Zero Count Samples vs. Ranked Median Counts",
           color = "Genes Passing Variance and Counts Filters")
    
    
    med_var_scatter / med_zeros_scatter #create figure with both scatter plots
    
  }
  
  output$count_summary_scatters <- renderPlot({
    
    input$do #connects slider input to action button
    counts_diagnostic_scatterplots(count_info_data = counts_data(),
                                   var_filter = isolate(input$var_range), 
                                   count_filter = isolate(input$count_range),
                                   color_1 = isolate(input$count_passed_filter_color),
                                   color_2 = isolate(input$count_failed_filter_color))
    
  })
  
    # counts_filtered_heatmap
    #  creates a clustered heatmap of counts remaining after filtering. Counts
    #are log10-transformed before plotting
    #  count_info_data(df): a data frame with counts information
    #  var_filter (int): a variance percentile value for which to filter the counts
    #  count_filter (int): a number of zero counts per gene to use for filtering
  
  
  counts_filtered_heatmap <- function(count_info_data, var_filter, count_filter) {
    
    #  filtering step, convert to matrix, do log10 transformation
    #filtering step
    #calculate row/gene variances
    gene_variances <- apply(count_info_data, 1, var)
    
    #convert slider value to probability value between 0 and 1
    var_filter_prob <- var_filter/100 
    
    #calculate variance threshold
    var_threshold <- quantile(gene_variances, var_filter_prob)
    
    filtered_counts <- count_info_data %>%
      filter(
        rowSums(. != 0) >= count_filter &
          apply(., 1, function(row) var(row) >= var_threshold)
      )
    
    
    # Extract the counts matrix without row names and column names
    counts_matrix <- as.matrix(filtered_counts[, -1])
    
    # Log10 transformation
    log10_filtered_matrix <- log10(1 + counts_matrix)
    
    # Add back row names and column names to the log10-transformed matrix
    rownames(log10_filtered_matrix) <- rownames(filtered_counts)  
    colnames(log10_filtered_matrix) <- colnames(counts_matrix)
    
    #heatmap
    heatmap.2(log10_filtered_matrix,
              key = TRUE,
              density.info = "none",
              trace = "none",
              
              # Viridis color palette
              col = colorRampPalette(c("#440154", "#3B528B", "#21908C", "#5DC963", "#FDE725"))(256),
              
              sepwidth = c(0.2, 0.2),         # Set width for column separations
              margins = c(10, 5),             # Adjust top and bottom margins
              main = "Log10-Transformed Counts",     # Main title
              cex.main = 3,                   # Font size
              key.xlab = "Samples",
              key.ylab = "Genes"
    )
  }
  
  output$count_heat <- renderPlot({
    input$do #connects slider input to action button
    counts_filtered_heatmap(count_info_data = counts_data(),
                            var_filter = isolate(input$var_range),
                            count_filter = isolate(input$count_range))  
  })
  
    # counts_pca_plot
    #  conducts principal components analysis (PCA) and generates a scatter plot
    # of PCA results for filtered counts
    #   data (tib): a (n x _S_) data set
    #   x_axis (string): which principal component to plot on the x-axis
    #   y_axis (string): which principal component to plot on the y-axis
    #   ggplot: scatter plot showing each sample in the first two PCs.
    # `plot_pca(data, meta, "Raw Count PCA")`
  counts_pca_plot <- function(count_info_data, x_axis, y_axis) {
    
    
    #do pca on transposed counts matrix (samples as rows, columns as genes)
    pca <- prcomp(
      t(count_info_data),
      center = TRUE,
      scale = FALSE
    )
    
    
    #Merging PCA output and metadata together by sample names
    pca_output <- data.frame(pca$x[,1:10])    
    
    
    #Include variance_explained in x and y axis labels:
    var_explained <- (pca$sdev)^2/sum((pca$sdev)^2) * 100
    
    #extract selected principal components for the x and y axes
    x_index <- as.numeric(sub("PC","",x_axis))
    y_index <- as.numeric(sub("PC","",y_axis))
    
    x_component <- var_explained[x_index]
    y_component <- var_explained[y_index]
    
    #make plot
    pca_gg <- ggplot(pca_output, aes(x=!!sym(x_axis),y=!!sym(y_axis))) +
      geom_point(color = "#0AADA8")+
      theme_bw()+
      theme(panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), axis.line = element_line(colour = "#050505"))+
      
      labs(x = paste0(x_axis,': ', round(x_component,1), '% variance'),
           y = paste0(y_axis, ': ', round(y_component,1), '% variance'),
           title = "Counts Principal Components Analysis Results")
    
    return(pca_gg)
  }
  
  
  output$count_pca <- renderPlot({
    input$pca_plot_do #connects slider input to action button
    counts_pca_plot(count_info_data = counts_data(), 
                    x_axis = isolate(input$pca_x_axis),
                    y_axis = isolate(input$pca_y_axis))
    
  })

  ##DIFFERENTIAL EXPRESSION ANALYSIS RESULTS FUNCTIONS
    # load_deseq_output
    #  loads the output statistics from a DESEQ2 differential expression analysis.
  
  load_deseq_output <- reactive({
    
    req(input$deseq_upload)
    deseq_df <- read.csv(input$deseq_upload$datapath, stringsAsFactors = FALSE, row.names = 1)
    
    return(deseq_df)
    
  })

    # deseq_results_table
    #  takes DESEQ2 differential expression analysis results, adds an additional column
    #indicating upregulation/downregulation/NS, and outputs a stortable table with gene search functionality
    #  deseq_mat (df): a data frame containing deseq2 results
    #  padj_threshold (num): an adjusted p-value number for which to filter the deseq results 
  
  
  deseq_results_table <- function(deseq_mat, padj_threshold){
    
    #load results
    #applying padj filter,adding additional up/down/ns column   
    annotated_deseq_res <- deseq_mat %>%
      filter(padj <= padj_threshold) %>% 
      mutate(Expression_Status = case_when(padj < padj_threshold & log2FoldChange > 0 ~ 'UP',
                                           padj < padj_threshold & log2FoldChange < 0 ~ 'DOWN',
                                           TRUE ~ 'NS')) 
    
    
    datatable(annotated_deseq_res, rownames = TRUE, options = list(pageLength = nrow(annotated_deseq_res)))
  }
  
  #sample data table
  output$deseq_tab <- renderDT({  
    input$deseq_do #link action button to adjusted pvalue slider
    deseq_results_table(deseq_mat = load_deseq_output(),
                        padj_threshold = isolate(input$padj_range))
    
  })

    # Function to plot the unadjusted p-values as a histogram
    #
    #   deseq_results (tib): Tibble with DESeq2 results
    #
    #   ggplot: a histogram of the raw p-values from the DESeq2 results
    #    
    #
    #    pval_plot <- plot_pvals(labeled_results)
  deseq_pval_histogram <- function(deseq_results) {
    
    pval_gg <- ggplot(deseq_results, aes(x = pvalue)) +
      geom_histogram(bins = 60, color='black', fill='lightblue2')+
      theme_bw()+
      labs(x = 'P-Value', y = 'Count', title = 'Histogram of Raw P-Values Obtained from DE Analysis (HD vs. Control)')
    pval_gg
  }
  
  output$pval_hist<- renderPlot({
    deseq_pval_histogram(deseq_results = load_deseq_output())
    
  })
  

    # Function to plot the log2FoldChange from DESeq2 results in a histogram
    #   labeled_results (tibble): Tibble with DESeq2 results 
    #   padj_threshold (float): threshold for considering significance (padj)
    #   ggplot: a histogram of log2FC values from genes significant at padj
    # threshold
  
  plot_log2fc <- function(deseq2_results, padj_threshold) {
    
    #filter input to only include results within defined padj_threshold
    
    filtered_results <- subset(deseq2_results, padj <= padj_threshold)
    
    #plotting filtered results
    gg2 <- ggplot(filtered_results, aes(x = log2FoldChange)) +
      geom_histogram(bins = 100, color='black', fill='cyan3')+
      theme_bw()+
      labs(x = 'log2FoldChange', y = 'Counts', title = 'Histogram of Log2FoldChanges for DE Genes (HD vs. Control)')
    gg2
    
  }
  
  output$logfc_hist<- renderPlot({
    input$deseq_do #links action button to pvalue slider
    plot_log2fc(deseq2_results = load_deseq_output(), 
                padj_threshold = isolate(input$padj_range))
    
  })
  
    # Volcano plot
    #   dataf The loaded data frame.
    #   slider A negative integer value representing the magnitude of
    # p-adjusted values to color. Most of our data will be between -1 and -300.
    #   color1 One of the colors for the points.
    #   color2 The other colors for the points. 
    #   A ggplot object of a volcano plot
    #   Returns a volcano plot of -log10(adjusted-p values) and log2Fold Change
    # of DESEQ2 results
    #    volcano_plot(df, "blue", "taupe")
  volcano_plot <- function(deseq2_results, slider, color1, color2) {
    gg <- ggplot(deseq2_results, aes(x = log2FoldChange, y = -log10(padj))) +
      geom_point(aes(color = ifelse(-log10(padj) < slider, "FALSE", "TRUE")), size = 0.5)+
      scale_color_manual(values = c("TRUE" = color2, "FALSE"= color1))+
      theme_bw()+
      geom_hline(yintercept = 0, linetype = "dashed") +
      
      labs(x = "log2FoldChange", y = "-log10(padj)", 
           title = 'Volcano plot of DESeq2 differential expression results',
           color = paste("-log10(padj) >=", slider))+
      theme(legend.position = "bottom")
    return(gg)
  }
  
  output$volc_plot<- renderPlot({
    
    input$volc_do #links action button to pvalue slider
    
    volcano_plot(deseq2_results = load_deseq_output(),
                 slider = isolate(input$volc_plot_range),
                 color1 = isolate(input$volc_col),
                 color2 = isolate(input$volc_col2)
    )
    
  })

  ##ReactomePA Reactome Pathway Analysis
  #Load the original data (CSV file)
  reactome_data <- reactive({
    req(input$reactome_upload)
    read.csv(input$reactome_upload$datapath)
  })
  
  #Extract relevant columns using org.Hs.eg.db
  library(org.Hs.eg.db)  # Load the annotation database
  
  #Process gene list to convert Ensembl IDs to Entrez IDs
  genelist_process <- reactive({
    df <- reactome_data()
    
    #Process gene list: Remove version numbers from Ensembl IDs
    genelist <- data.frame(
      GeneID = gsub("\\..*$", "", df[, 1]),  # Remove version numbers
      FC = df[, 6],  # Fold Change
      padj = df[, 10]  # Adjusted p-value
    )
    
    #Use org.Hs.eg.db to map Ensembl IDs to Entrez IDs
    genelist$entrezgene_id <- mapIds(
      org.Hs.eg.db,
      keys = genelist$GeneID,
      column = "ENTREZID",
      keytype = "ENSEMBL",
      multiVals = "first"
    )
    
    #Filter the gene list for adjusted p-value and fold change criteria
    final <- genelist %>%
      filter(!is.na(entrezgene_id)) %>%
      filter(padj < 0.05 & abs(FC) > 1.5)
    
    final$entrezgene_id  #Return the Entrez IDs for pathway analysis
  })
  
  #Perform Reactome pathway enrichment
  pathway_results <- eventReactive(input$run_analysis, {
    enrichPathway(
      gene = genelist_process(),
      organism = "human",
      pvalueCutoff = 0.05,
      readable = TRUE
    )
  })
  
  # Render pathway table
  output$pathway_table <- renderDT({
    req(pathway_results())
    results_df <- as.data.frame(pathway_results()) %>%
      arrange(p.adjust)  # Sort by p.adjust ascending
    datatable(results_df, options = list(scrollX = TRUE))
  })
  
  output$barplot <- renderPlot({
    req(pathway_results())
    data <- as.data.frame(pathway_results()) %>%
      arrange(p.adjust) %>%  # Sort by p.adjust ascending
      mutate(Description = fct_reorder(Description, p.adjust))
    
    ggplot(data, aes(x = Description, y = p.adjust)) +
      geom_bar(stat = "identity", aes(fill = p.adjust), show.legend = FALSE) +
      scale_fill_gradient(low = "#5DC863", high = "#440154") +
      coord_flip() +
      labs(
        title = "Top Enriched Pathways",
        x = "Pathways",
        y = "Adjusted P-value"
      ) +
      theme_minimal()
  })
    output$dotplot <- renderPlot({
      req(pathway_results())
      data <- as.data.frame(pathway_results()) %>%
        arrange(p.adjust) %>%  #Sort by p.adjust ascending
        mutate(Description = fct_reorder(Description, p.adjust))
      
      ggplot(data, aes(x = Description, y = p.adjust, size = Count, color = p.adjust)) +
        geom_point(alpha = 0.8) +
        scale_color_gradient(low = "#5DC863", high = "#440154") +
        coord_flip() +
        labs(
          title = "Pathway Enrichment Dotplot",
          x = "Pathways",
          y = "Adjusted P-value"
        ) +
        theme_minimal()
    })
  
    output$cnetplot <- renderPlot({
      req(pathway_results())
      data <- as.data.frame(pathway_results())
      
      ggplot(data, aes(x = Description, y = p.adjust, 
                       size = as.numeric(sub("/.*", "", GeneRatio)), 
                       color = p.adjust)) +
        geom_point(alpha = 0.7, shape = 9) +
        scale_color_gradientn(
          colors = c("#440154", "#3B528B", "#21908D", "#5DC863", "#FDE725"),
          name = "Adjusted P-value",
          guide = guide_colorbar(reverse = TRUE)
        ) +
        scale_size_continuous(
          name = "Gene Ratio",
          range = c(3, 10),
          breaks = c(5, 10, 15, 20, 25)
        ) +
        coord_flip() +
        labs(
          title = "Pathway Enrichment Analysis",
          subtitle = "Pathways Ranked by Adjusted P-value",
          x = "Pathways",
          y = "Adjusted P-value"
        ) +
        theme_minimal(base_family = "sans") +
        theme(
          plot.title = element_text(face = "bold", size = 16, color = "#2C3E50"),
          plot.subtitle = element_text(color = "#7F8C8D"),
          axis.text.y = element_text(size = 10, color = "#34495E"),
          axis.title = element_text(face = "bold", color = "#2C3E50"),
          legend.position = "right",
          legend.title = element_text(face = "bold"),
          panel.grid.major.y = element_line(linetype = "dotted", color = "#BDC3C7"),
          plot.background = element_rect(fill = "#F4F6F7"),
          panel.background = element_rect(fill = "#F4F6F7")
        ) +
        guides(
          color = guide_colorbar(order = 1),
          size = guide_legend(order = 2))
    })
    
   
}
# Run the application 
 shinyApp(ui = ui, server = server)
