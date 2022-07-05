library(shiny)
library(stringi)

description <-
  "The first step is to upload a FASTA file of protein sequences.
These sequences will be reverse translated into DNA sequence
and then codon optimized with ThermoFisher GeneArt codon optimization
tool to the species you select. The algorithm will framgent the DNA
sequences accoring to the criteria you select, using defined cut sites
and length constraints. The algorithm then finds the set of cuts that
created the minimal total number of fragments using only one of the
defined cut sequences per protein. The output is a codon optimized
sequence for each proten, fragment DNA sequences and sets of fragments
that reconstruct the DNA sequence to code for each input protein.
The computational time is highly dependent on the length of the
protein and the minimal fragment length. If you submit long sequences,
>1400 amino acids, it is recommended to increase the minimal fragment length to
800-900."

ui <- fluidPage(

  # Application title
  titlePanel("Fragment design for Golden Gate cloning"),
  h4("This tool creates DNA fragments that can be used for
      Golden Gate cloning from a set of protein sequences."),

  # Input data is set and submitted in this sidebar
  sidebarLayout(
    sidebarPanel(
      p(description),

      # Fasta file upload field
      fileInput("protein_fasta", "Upload a protein fasta file"),

      # Selection of species to perform codon optimization for
      selectInput("species",
                  "Species to codon optimize for",
                  c("Arabidopsis thaliana",
                    "Bacillus subtilis",
                    "Bos taurus",
                    "Brassica napus",
                    "Caenorhabditis elegans",
                    "Cricetulus griseus",
                    "Drosophila melanogaster",
                    "Escherichia coli",
                    "Glycine max",
                    "Homo sapiens",
                    "Hordeum vulgare",
                    "Lycopersicon esculentum",
                    "Mus musculus",
                    "Nicotiana benthamiana",
                    "Nicotiana tabacum",
                    "Oryctolagus cuniculus",
                    "Oryza sativa",
                    "Pichia pastoris",
                    "Saccharomyces cerevisiae",
                    "Schizosaccharomyces pombe",
                    "Spodoptera frugiperda",
                    "Synechococcus elongatus",
                    "Vaccinia virus",
                    "Zea mays"
                    ),
                  selected = "Homo sapiens"
                  ),

      # Fragment length input
      sliderInput("fragment_range",
                  "Fragment length range",
                  min = 300,
                  max = 1200,
                  value = c(300, 1000),
                  step = 10),
      # Sequence to exclude in codon optimization
      textInput("prevent_seq",
                "Sequence to prevent after codon optimization",
                "GGTCTCA"),

      # Overhangs to use
      textInput("overhangs",
                "Allowed fragment overhangs",
                value = "CCGA,GGGC,AGAA,CTTA,TCAA"),

      # Overhangs to use
      textInput("api_token",
                "Thermo Fisher API token",
                value = ""),


      # Submit/Go Button
      actionButton("go", "Go"),

      # Progress
      htmlOutput("progress")
    ),

    # Results are displayed in the main panel
    mainPanel(
      # This section contains 3 tabs, one to display uploaded sequences in the
      # fasta file, one for codon optimized sequences and one for the
      # final design
      tabsetPanel(
        tabPanel(
          "Uploaded fasta file",
          h4("Uploaded fasta file"),
          uiOutput("input_fasta")
        ),
        tabPanel(
          "Codon optimized sequence",
          h4("Codon optimized sequence"),
          p("Codon optimized sequence will appear below
             when calculations are complete"),
          uiOutput("download_codon_opt"),
          htmlOutput("python_stdout")
        ),
        tabPanel(
          "Fragment design",
          h4("Fragment design"),
          p("Fragment design will appear below when calculations are complete"),
          uiOutput("download_fragment_design"),
          uiOutput("fragment_design")
        )
      )
    )
  )
)


server <- function(input, output) {
  # Input fasta file
  output$species <- renderText({
    input$species
  })

  output$input_fasta <- renderUI({
    in_file <- input$protein_fasta
    if (is.null(in_file)) {
      return(NULL)
    }
    raw_text <- readLines(in_file$datapath) # get raw text
    split_text <- stri_split(str = raw_text, regex = "\\n")
    replaced_text <- lapply(split_text, p)

    return(replaced_text)
  })

  # Codon optimization
  python_stdout <- eventReactive(input$go, {

    in_file <- input$protein_fasta
    min_length <- round(input$fragment_range[1] / 3)
    max_length <- round(input$fragment_range[2] / 3)
    species <- input$species
    prevent_seq <- input$prevent_seq
    api_token <- input$api_token

    if (is.null(in_file)) {
      return(NULL)
    }

    system(paste0("dos2unix ", in_file$datapath))
    system(paste0("python3 ./src/codon_optimize.py ",
                  in_file$datapath, " ",
                  min_length, " ",
                  max_length, " ",
                  paste0("'", species, "'"),
                  " ", prevent_seq,
                  " ",
                  api_token,
                  " 2>&1"),
           intern = TRUE,
           ignore.stderr = TRUE)
  })

  output$python_stdout <- renderText({
    python_stdout()
  })

  output$download_output <- downloadHandler(
    filename = "CodonOptimizedSeq.fa",
    content = function(file) {
      writeLines(python_stdout(), file)
    })

  output$download_codon_opt <- renderUI({
    if (!is.null(python_stdout())) {
      downloadButton("download_output", "Download Output File")
    }
  })

  # Fragmentation
  fragment_design <- reactive({
    if (is.null(python_stdout()) |
          is.null(input$overhangs) |
          input$overhangs == "") {
      return(NULL)
    }

    opt_seq_file <- tempfile()
    writeLines(python_stdout(), opt_seq_file)
    system(paste0("python3 ./src/find_optimal_fragments.py ",
                  opt_seq_file,
                  " ",
                  input$overhangs,
                  " ",
                  input$fragment_range[1],
                  " ",
                  input$fragment_range[2],
                  " ",
                  input$prevent_seq),
           intern = TRUE)
  })

  output$fragment_design <- renderUI({
    req(fragment_design())
    split_text <- stringi::stri_split(str = fragment_design(), regex = "\\n")
    replaced_text <- lapply(split_text, p)
    return(replaced_text)
  })

  output$download_fragments <- downloadHandler(
    filename = "FragmentDesign.txt",
    content = function(file) {
      writeLines(fragment_design(), file)
    })

  output$download_fragment_design <- renderUI({
    if (!is.null(python_stdout())) {
      downloadButton("download_fragments", "Download Output File")
    }
  })

  # Progress messages
  output$progress <- renderText({
    if (!is.null(input$protein_fasta) &
        input$go &
        !is.null(python_stdout()) &
        !is.null(fragment_design())) {
      return("Upload complete<br>
              Codon optimization complete<br>
              Fragment design complete")

    } else if (!is.null(input$protein_fasta) &
               input$go &
               !is.null(python_stdout())) {
      return("Upload complete<br>Codon optimization complete")
    } else if (!is.null(input$protein_fasta) & input$go) {
      return("Upload complete")
    }
  })
}

# Run the application
shinyApp(ui = ui, server = server)
