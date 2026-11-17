# Lab 1 - MT Alignment

library(shiny)
library(bslib)
library(Biostrings)
library(msa)

base_url <- "https://raw.githubusercontent.com/nzupp/Evolution-Labs/main/Lab1/"

aln_file <- tempfile(fileext = ".fasta")
download.file(paste0(base_url, "mammal_mt_aligned.fasta"), aln_file, quiet = TRUE)
aligned_seqs <- readDNAStringSet(aln_file)
aligned_chars <- as.character(aligned_seqs)
all_seqs <- DNAStringSet(gsub("-", "", aligned_chars))
names(all_seqs) <- names(aligned_seqs)

species_info <- data.frame(
  scientific = c("Homo sapiens", "Macaca mulatta",
                 "Canis lupus", "Panthera onca", "Ursus maritimus",
                 "Balaenoptera musculus", "Tursiops truncatus",
                 "Pteropus vampyrus", "Mus musculus",
                 "Loxodonta africana",
                 "Macropus giganteus", "Ornithorhynchus anatinus"),
  common = c("Human", "Rhesus Macaque",
             "Wolf", "Jaguar", "Polar Bear",
             "Blue Whale", "Bottlenose Dolphin",
             "Flying Fox", "Mouse",
             "African Elephant",
             "Eastern Grey Kangaroo", "Platypus"),
  stringsAsFactors = FALSE
)

img_urls <- setNames(
  paste0(base_url, "mammal_images/", gsub(" ", "_", species_info$scientific), ".png"),
  species_info$scientific
)

nt_colors <- c(A = "#64f73f", C = "#3c88ee", G = "#eb413c", T = "#ffb340", U = "#ffb340")

render_alignment <- function(chars) {
  rows <- lapply(names(chars), function(nm) {
    res <- strsplit(chars[[nm]], "")[[1]]
    col <- nt_colors[toupper(res)]
    col[is.na(col)] <- "transparent"
    spans <- paste0('<span style="background:', col, '">', res, "</span>", collapse = "")
    div(
      class = "d-flex text-nowrap",
      style = "line-height: 24px;",
      div(
        class = "d-flex align-items-center gap-2 bg-white",
        style = "position: sticky; left: 0; width: 250px; flex: none; z-index: 1;",
        if (!is.na(img_urls[nm])) img(src = img_urls[nm], height = 20),
        em(nm)
      ),
      div(style = "font-family: monospace;", HTML(spans))
    )
  })
  div(class = "border rounded p-2 overflow-auto", rows)
}

full_alignment_ui <- render_alignment(aligned_chars)

ui <- page_fluid(
  theme = bs_theme(primary = "#566db0"),
  title = "Mammal Mitochondrial Genome Alignment",
  h1("Mammal Mitochondrial Genome Alignment", class = "my-4"),
  h2("Introduction"),
  p("Good pedagogy here"),
  p("Here we will compare mitochondrial genome sequences from twelve species spanning class Mammalia."),
  h2("Full Alignment"),
  full_alignment_ui,
  h2("Select Subset", class = "mt-4"),
  p("Optionally select species to re-align a subset."),
  uiOutput("species_grid"),
  input_task_button("align", "Re-Align Subset", label_busy = "Aligning...", class = "my-3"),
  uiOutput("subset_alignment")
)

server <- function(input, output, session) {
  selected <- reactiveVal(character(0))
  subset_aln <- reactiveVal(NULL)
  
  observeEvent(input$card_click, {
    sp <- input$card_click
    cur <- selected()
    selected(if (sp %in% cur) setdiff(cur, sp) else c(cur, sp))
  })
  
  output$species_grid <- renderUI({
    cards <- lapply(seq_len(nrow(species_info)), function(i) {
      sp <- species_info$scientific[i]
      is_sel <- sp %in% selected()
      card(
        class = if (is_sel) "border-2 border-primary bg-primary-subtle" else "border-2 opacity-50",
        style = "cursor: pointer; user-select: none;",
        onclick = sprintf("Shiny.setInputValue('card_click', '%s', {priority: 'event'})", sp),
        card_body(
          class = "d-flex align-items-center gap-3",
          img(src = img_urls[[sp]], width = 40, height = 40, style = "object-fit: contain;"),
          div(
            div(class = "fw-bold fst-italic", sp),
            div(class = "small text-muted", species_info$common[i])
          )
        )
      )
    })
    do.call(layout_column_wrap, c(list(width = 1/3), cards))
  })
  
  observeEvent(input$align, {
    if (length(selected()) < 2) {
      showNotification("Select at least two species.", type = "warning")
      return()
    }
    seqs <- all_seqs[names(all_seqs) %in% selected()]
    aln <- msa(seqs, method = "ClustalOmega")
    subset_aln(as.character(unmasked(as(aln, "DNAMultipleAlignment"))))
  })
  
  output$subset_alignment <- renderUI({
    req(subset_aln())
    render_alignment(subset_aln())
  })
}

shinyApp(ui, server)


