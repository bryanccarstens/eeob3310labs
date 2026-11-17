# Lab 1 - BRCA1 Alignment

library(shiny)
library(bslib)
library(Biostrings)
library(msa)

base_url <- "https://github.com/bryanccarstens/eeob3310labs/tree/main/Apps/Lab1/"

aln_file <- tempfile(fileext = ".fasta")
download.file(paste0(base_url, "primate_BRCA1_aligned.fasta"), aln_file, quiet = TRUE)
aligned_seqs <- readAAStringSet(aln_file)
aligned_chars <- as.character(aligned_seqs)
all_seqs <- AAStringSet(gsub("-", "", aligned_chars))
names(all_seqs) <- names(aligned_seqs)

species_info <- data.frame(
  scientific = c("Homo sapiens", "Pan troglodytes", "Gorilla gorilla",
                 "Pongo pygmaeus", "Hylobates lar", "Macaca mulatta",
                 "Papio anubis", "Saimiri sciureus", "Callithrix jacchus"),
  common = c("Human", "Chimpanzee", "Western Gorilla",
             "Bornean Orangutan", "White-handed Gibbon", "Rhesus Macaque",
             "Olive Baboon", "Common Squirrel Monkey", "Common Marmoset"),
  group = c("Great Ape", "Great Ape", "Great Ape",
            "Great Ape", "Lesser Ape", "Old World Monkey",
            "Old World Monkey", "New World Monkey", "New World Monkey"),
  stringsAsFactors = FALSE
)

img_urls <- setNames(
  paste0(base_url, "primate_images/", gsub(" ", "_", species_info$scientific), ".png"),
  species_info$scientific
)

aa_colors <- c(
  A = "#80a0f0", I = "#80a0f0", L = "#80a0f0", M = "#80a0f0", F = "#80a0f0",
  W = "#80a0f0", V = "#80a0f0", C = "#f08080",
  K = "#f01505", R = "#f01505", E = "#c048c0", D = "#c048c0",
  N = "#15c015", Q = "#15c015", S = "#15c015", T = "#15c015",
  G = "#f09048", P = "#c0c000", H = "#15a4a4", Y = "#15a4a4"
)

render_alignment <- function(chars) {
  rows <- lapply(names(chars), function(nm) {
    res <- strsplit(chars[[nm]], "")[[1]]
    col <- aa_colors[res]
    col[is.na(col)] <- "transparent"
    spans <- paste0('<span style="background:', col, '">', res, "</span>", collapse = "")
    div(
      class = "d-flex text-nowrap",
      style = "line-height: 24px;",
      div(
        class = "d-flex align-items-center gap-2 bg-white",
        style = "position: sticky; left: 0; width: 210px; flex: none; z-index: 1;",
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
  title = "Primate BRCA1 Protein Alignment",
  h1("Primate BRCA1 Protein Alignment", class = "my-4"),
  h2("Introduction"),
  p("Good pedagogy here"),
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
            div(class = "small text-muted",
                paste0(species_info$common[i], " (", species_info$group[i], ")"))
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
    subset_aln(as.character(unmasked(as(aln, "AAMultipleAlignment"))))
  })
  
  output$subset_alignment <- renderUI({
    req(subset_aln())
    render_alignment(subset_aln())
  })
}

shinyApp(ui, server)
