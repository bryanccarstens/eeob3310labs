# Lab 3 - Genetic Diversity in Hoary Bats

library(shiny)
library(bslib)

base_url <- "https://github.com/bryanccarstens/eeob3310labs/tree/main/Data/Lab3"
data_url <- paste0(base_url, "HO_Data/")
img_url  <- paste0(base_url, "hoary.jpg")

snp  <- read.delim(paste0(data_url, "SNPMatrix_100.55.txt"),
                   check.names = FALSE, stringsAsFactors = FALSE)
mono <- read.table(paste0(data_url, "Monomorphics_100.txt"),
                   header = FALSE, stringsAsFactors = FALSE)

ids  <- snp[, 1]
snp  <- as.matrix(snp[, -1])
a1   <- snp[seq(1, nrow(snp), 2), ]
a2   <- snp[seq(2, nrow(snp), 2), ]
bats <- ids[seq(1, nrow(snp), 2)]
rownames(a1) <- rownames(a2) <- bats

n_poly      <- length(unique(colnames(snp)))
n_mono      <- nrow(mono)
L           <- nchar(mono[1, 1])
total_sites <- (n_poly + n_mono) * L

ho_per_locus <- function(a1, a2) {
  sapply(seq_len(ncol(a1)), function(j) {
    ok <- a1[, j] != "-" & a2[, j] != "-"
    if (sum(ok) == 0) return(NA)
    sum(a1[ok, j] != a2[ok, j]) / sum(ok)
  })
}

pi_per_locus <- function(a1, a2) {
  sapply(seq_len(ncol(a1)), function(j) {
    x <- c(a1[, j], a2[, j]); x <- x[x != "-"]; k <- length(x)
    if (k <= 1) return(0)
    f <- table(x)
    k / (k - 1) * (1 - sum((f / k)^2))
  })
}

ho_j <- ho_per_locus(a1, a2)
pi_j <- pi_per_locus(a1, a2)
Ho   <- mean(ho_j, na.rm = TRUE)
Pi   <- sum(pi_j) / total_sites

af <- do.call(rbind, lapply(seq_len(ncol(a1)), function(j) {
  x <- c(a1[, j], a2[, j]); x <- x[x != "-"]
  f <- sort(table(x) / length(x), decreasing = TRUE)
  data.frame(locus = colnames(a1)[j], snp = j,
             major = names(f)[1], major_freq = round(as.numeric(f[1]), 3),
             minor = paste(names(f)[-1], collapse = "/"),
             minor_freq = paste(round(f[-1], 3), collapse = "/"),
             stringsAsFactors = FALSE)
}))

sample_sizes <- c(5, 10, 15, 20)
sample_sizes <- sample_sizes[sample_sizes <= length(bats)]

compact_css <- "
.bslib-value-box .value-box-area { padding: 0.5rem 0.9rem !important; }
.bslib-value-box .value-box-title { font-size: 0.85rem; margin-bottom: 0.1rem; }
.bslib-value-box .value-box-value { font-size: 1.3rem; }
"

ui <- page_fluid(
  theme = bs_theme(primary = "#566db0"),
  tags$head(tags$style(compact_css)),
  title = "Genetic Diversity in Hoary Bats",
  h1("Genetic Diversity in Hoary Bats", class = "my-4"),
  img(src = img_url, class = "mb-4", style = "width: 40%; min-width: 220px;"),
  
  h2("Illumina Sequencing"),
  p("stuff about Illumina here"),
  
  h2("Stuff about variation here"),
  p("Variation stuff"),
  
  h2("Diversity Estimates", class = "mt-4"),
  layout_column_wrap(
    width = "240px", fixed_width = TRUE, fill = FALSE,
    value_box(title = "Observed Heterozygosity (Ho)",
              value = round(Ho, 5), theme = "primary", height = "80px"),
    value_box(title = "Nucleotide Diversity (π)",
              value = round(Pi, 5), theme = "primary", height = "80px")
  ),
  p(class = "small text-muted mt-2",
    sprintf("%d individuals, %s total sites.", length(bats),
            format(total_sites, big.mark = ","))),
  
  h2("Effect of Sample Size", class = "mt-4"),
  p("Randomly subsample bats and recompute diversity. Resample to see how much the estimates change."),
  actionButton("resample", "Resample", class = "btn-primary mb-3"),
  tableOutput("subsample_table"),
  
  h2("Allele Frequencies", class = "mt-4"),
  p("A random sample of 10 SNPs and their allele frequencies."),
  actionButton("more_af", "Show another 10", class = "btn-primary mb-3"),
  tableOutput("af_table"),
  
  h2("Hardy-Weinberg Equilibrium", class = "mt-4"),
  p("Pedagogy here"),
  
  h2("Genetic Drift", class = "mt-4"),
  tags$a(
    href = "https://media.hhmi.org/biointeractive/click/population-genetics-explorer/replicated",
    target = "_blank",
    class = "btn btn-primary",
    "Open Population Genetics Explorer"
  ),
  div(style = "height: 80px;")
)

server <- function(input, output, session) {
  
  subsample <- eventReactive(input$resample, {
    if (input$resample == 0) set.seed(1)  # first draw matches the original Rmd
    rows <- lapply(sample_sizes, function(n) {
      s <- sample(bats, n)
      h <- ho_per_locus(a1[s, ], a2[s, ])
      p <- pi_per_locus(a1[s, ], a2[s, ])
      data.frame(n = n,
                 Ho = mean(h, na.rm = TRUE), Ho_var = var(h, na.rm = TRUE),
                 pi = sum(p) / total_sites,  pi_var = var(p))
    })
    do.call(rbind, rows)
  }, ignoreNULL = FALSE)
  
  output$subsample_table <- renderTable(subsample(), digits = 6)
  
  af_sample <- eventReactive(input$more_af, {
    af[sample(nrow(af), min(10, nrow(af))), ]
  }, ignoreNULL = FALSE)
  
  output$af_table <- renderTable(af_sample(), digits = 3)
}

shinyApp(ui, server)
