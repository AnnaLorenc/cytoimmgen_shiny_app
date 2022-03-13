library(shiny)
library(pheatmap)
library(tidyverse)
library(RColorBrewer)

#Data preprocessing
#see retrace_Blagojes_processing.R

main_dir="data"
#averaged expression values -files
files_with_matrices <- list(naive="data/matrices_for_heatmaps/naive_mean_expression_per_window_full.txt",
                        memory="data/matrices_for_heatmaps/memory_mean_expression_per_window_full.txt")

#eQTL genes-files
files_with_lists_of_genes <- list(naive=
                                            "data/genes_to_use_as_starting_lists/naive_genelist_dynamic_clean.csv",
                            memory="data/genes_to_use_as_starting_lists/memory_genelist_dynamic_clean.csv")

window_plot <- read_csv("data/matrices_for_heatmaps/data_withpseud.csv")


#averaged expression values -data
data_with_matrices <- sapply(files_with_matrices, function(x) {
    (read_table(x, col_names = TRUE))%>%
        select(index, paste0("win", 1:10))}, simplify = FALSE, USE.NAMES = TRUE)

#eQTL genes -data
data_with_list_of_genes <- sapply(files_with_lists_of_genes, function(x) {
    (read_csv(x))}, simplify = FALSE, USE.NAMES = TRUE)

#here we assume we want only genes which were found to be eQTLs
final_genes_for_big_heatmap <- sapply(names(data_with_list_of_genes), function(for_heatmap) intersect(data_with_list_of_genes[[for_heatmap]]$gene, data_with_matrices[[for_heatmap]]$index), simplify = FALSE, USE.NAMES = TRUE)



#here we assume we want only genes which were found to be eQTLs
#genes_for_dropdown_menu <-sort(unique(unlist(final_genes_for_big_heatmap)))
genes_for_dropdown_menu <- sort(unique(c(data_with_matrices[[1]]$index,data_with_matrices[[1]]$index)))

main_text="Change of gene expression in CD4 cells after activation.
We performed scRNAseq of naive and memory CD4 cells from 119 individuals.
We analysed these two types of CD4 cells before (0h) and at several timepoints (16h, 40h, 5days) after activation.
We located each cell on a pseudotime trajectory, from resting to highly activated cells. Here, we present average expression values for each gene (across individuals and cells) in relation to pseudotime - each of 10 pseudotime windows comprises 10% of all naive or memory cells. To relate pseudotime windows to activation timeline, underneath the heatmaps we show number of cells originating from each experimental timepoint."

upper_heatmap_text="For each gene, expression is normalised across pseudotime windows (expression values are not comparable between genes, only temporal trends)."
  
lower_heatmap_text="Expression values might be compared between genes. Normalised expression values in log10 scale. Download  the underlying data with the button on the left."  

ui = fluidPage(
    
        titlePanel("Gene expression in CD4 T cells after activation"),
    
        fluidRow(column(7,
    
          p( main_text, style = "font-size:18px;"),
          div("Detailed methods description is included  in ",
            a("preprint", href="https://doi.org/10.1101/2021.12.06.470953", target="_blank",style = "font-size:18px;"),
            style = "font-size:18px;"),
          p("Please resize the window if the plots do not show",style = "font-size:18px;color:red;")),
        ),
        
    sidebarLayout(
        
               sidebarPanel(width=2,
                   selectInput("assayNames",
                               p("Select cell type",style = "font-size:15px"),
                               choices = names(files_with_matrices)
                   ),
                   selectizeInput("geneName",
                               p("Select at least 2 genes manually. Start typing to see the options.",style = "font-size:15px"),
                               choices = NULL,
                               selected = final_genes_for_big_heatmap[[1]][1:3],
                               multiple = TRUE),
                   
                   fileInput("file_with_genenames",
                             p("...or provide a file with one gene symbol per line.",style = "font-size:15px")),
                   
                   downloadButton("downloadData", "Download") 
                   
  
     ),

     mainPanel(
               h3("Expression - normalised counts"),
               p(lower_heatmap_text, style = "font-size:15px"),
               plotOutput("naiveheatmap"),
               
               h3("Expression normalised within each genes"),
               p(upper_heatmap_text, style = "font-size:15px"),
               
               plotOutput("subsetmap"),
               plotOutput("celltypes_across_windows")
                 )

    )                    
)

server = function(input, output, session) {
    
    updateSelectizeInput(session, 'geneName', choices = genes_for_dropdown_menu, server = TRUE)
    
  data_for_window_plot <- reactive({window_plot%>%filter(cell_type==input$assayNames)})
    
    data_for_small_heatmap <-
        reactive({data <- data_with_matrices[[input$assayNames]]
        matrix_data <- data %>%select(-c(index))%>% data.matrix()
        rownames(matrix_data) <- data$index
        return(matrix_data)
            })
    
    
    gene_name_from_a_menu <- reactive({
        if(is.null(input$geneName)| length(input$geneName)<2){
            data_with_matrices[[1]][1:2,]$index
        }else{input$geneName}
    })
    
    gene_names_from_a_file <- reactive({
        inFile <- input$file_with_genenames

        if (is.null(inFile)){
            return(NULL)
        }else{

            return(
                intersect(unlist(read.table(inFile$datapath)),
                                 rownames(data_for_small_heatmap()))
                )}
    })

        
    
    small_heatmap <-reactive({
        if(is.null(gene_names_from_a_file())|length(gene_names_from_a_file())<2){
       #     gene_names_possible_to_plot = gene_name_from_a_menu()
        gene_names_possible_to_plot <-  gene_name_from_a_menu()
         }else{
             gene_names_possible_to_plot=gene_names_from_a_file()
       }
        return(data_for_small_heatmap()[gene_names_possible_to_plot,])
        
    })
    

    output$subsetmap = renderPlot( 
        pheatmap(small_heatmap() ,
                 cluster_cols=FALSE, cluster_rows=FALSE, show_rownames=TRUE, scale="row", show_colnames=TRUE,,fontsize_col=20, fontsize_row = 13,angle_col = 0))#, height=300)
    
   output$naiveheatmap = renderPlot( 
          pheatmap(log10(small_heatmap() +0.000001) ,
                   cluster_cols=FALSE, cluster_rows=FALSE, show_rownames=TRUE, scale="none", show_colnames=TRUE,fontsize_col=20, fontsize_row = 13,angle_col = 0,col= colorRampPalette(brewer.pal(8, "BuPu"))(25)))#, height=300)

    output$celltypes_across_windows <- renderPlot(   ggplot(data_for_window_plot()) +
                                                      geom_area(aes(x=factor(pseud_wind,levels=paste0("win", 1:10)),
                                                                                y=cells, fill=Time_point, group=Time_point))+
                                                      theme_classic()+guides(fill="none")+xlab("")+scale_y_continuous(position = "right")+
                                                      facet_wrap(~Time_point, ncol=1))
   
   # Downloadable csv of selected dataset ----
   output$downloadData <- downloadHandler(
     filename = function() {
       paste(input$assayNames, ".csv", sep = "")
     },
     content = function(file) {
       write.csv(small_heatmap(), file, row.names = TRUE)
     }
   )
}

shinyApp(ui, server)