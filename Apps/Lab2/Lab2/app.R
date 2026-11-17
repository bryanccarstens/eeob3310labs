# Lab 2 - Radiometric Dating with IsoplotR

library(shiny)
library(bslib)
library(IsoplotR)

base_url <- "https://raw.githubusercontent.com/nzupp/Evolution-Labs/main/Lab2/"
img_url  <- paste0(base_url, "tb.jpg")

# Melchor et al. (2023): Terror Bird Sample C893
C893 <- rbind(
  c(0.01059453,0.00123637, 0.00116294,0.00006392, 0.06491293,0.00797675),
  c(0.00783476,0.00065408, 0.00118553,0.00003199, 0.04676350,0.00359908),
  c(0.01404437,0.00246090, 0.00119645,0.00009993, 0.08379147,0.01600426),
  c(0.00757531,0.00122958, 0.00119979,0.00005128, 0.04571039,0.00747712),
  c(0.01054545,0.00130929, 0.00120152,0.00005469, 0.06385782,0.00732565),
  c(0.00838480,0.00087255, 0.00120603,0.00003637, 0.04873122,0.00464638),
  c(0.00910423,0.00159258, 0.00120890,0.00009885, 0.05723694,0.01061144),
  c(0.00871621,0.00114461, 0.00122649,0.00006270, 0.05171649,0.00708740),
  c(0.00981905,0.00135750, 0.00123154,0.00006710, 0.05935341,0.01079807),
  c(0.01112948,0.00127776, 0.00123943,0.00006076, 0.06562737,0.00793142),
  c(0.00690803,0.00089396, 0.00124082,0.00006078, 0.03979429,0.00529944),
  c(0.00705550,0.00102591, 0.00124382,0.00006095, 0.04122998,0.00620584),
  c(0.01421014,0.00155341, 0.00125949,0.00006264, 0.08242323,0.00941422),
  c(0.00918893,0.00117303, 0.00126369,0.00007560, 0.05291751,0.00705122),
  c(0.01411747,0.00198310, 0.00127483,0.00006782, 0.08260222,0.01259772),
  c(0.00811268,0.00107844, 0.00130899,0.00007362, 0.04525619,0.00623449),
  c(0.01022219,0.00135343, 0.00133066,0.00009408, 0.05604869,0.00777958),
  c(0.01536299,0.00187293, 0.00133471,0.00006605, 0.08143233,0.00961059),
  c(0.00853645,0.00083095, 0.00135121,0.00004638, 0.04620153,0.00458553),
  c(0.02134814,0.00222678, 0.00135567,0.00009763, 0.11263227,0.01288018),
  c(0.00830456,0.00130819, 0.00136317,0.00007506, 0.04455572,0.00730852),
  c(0.01937069,0.00279554, 0.00144525,0.00010747, 0.09991406,0.01587200),
  c(0.01707393,0.00243746, 0.00145404,0.00006341, 0.08767800,0.01420042),
  c(0.01208045,0.00147190, 0.00148645,0.00007225, 0.05949306,0.00759592),
  c(0.00615633,0.00131914, 0.00151504,0.00009291, 0.02939965,0.00648556),
  c(0.06215198,0.00724045, 0.00171724,0.00008557, 0.26241958,0.03082849),
  c(0.06946469,0.00472452, 0.00173371,0.00010485, 0.29090176,0.02684098),
  c(0.01624965,0.00390405, 0.00186456,0.00035777, 0.07464943,0.01941652),
  c(0.26780920,0.02578520, 0.00341201,0.00034392, 0.55086127,0.05263143),
  c(0.09834125,0.01222332, 0.01554128,0.00065949, 0.04913779,0.00604558)
)
colnames(C893) <- c("Pb207U235", "errPb207U235", "Pb206U238",
                    "errPb206U238", "Pb207Pb206", "errPb207Pb206")

grain_ids <- seq_len(nrow(C893))

make_upb <- function(rows) {
  d <- C893[rows, , drop = FALSE]
  as.UPb(cbind(d[, 1], d[, 2], d[, 3], d[, 4], rep(0.5, nrow(d))), format = 1)
}

ui <- page_fluid(
  theme = bs_theme(primary = "#566db0"),
  title = "Radiometric Dating with IsoplotR",
  h1("Radiometric Dating with IsoplotR", class = "my-4"),
  div(
    class = "d-flex align-items-center gap-4 mb-4",
    img(src = img_url, style = "width: 25%; min-width: 180px;"),
    div(
      h2("Introduction"),
      p("Good pedagogy here"),
      p("Here we will analyze U-Pb zircon data from sample C893 (Melchor et al., 2023), ",
        "associated with a terror bird fossil, to estimate its age.")
    )
  ),
  
  h2("Data"),
  p("Each row is one zircon grain. Ratios are given with their 1σ errors."),
  div(class = "border rounded p-2 overflow-auto mb-4", style = "max-height: 300px;",
      tableOutput("data_table")),
  
  h2("Exclude Grains"),
  p("Optionally select grains to exclude from the analysis. Both plots update automatically."),
  selectizeInput("exclude", NULL, choices = grain_ids, multiple = TRUE,
                 options = list(placeholder = "No grains excluded")),
  
  h2("Concordia Plot", class = "mt-4"),
  plotOutput("concordia_plot", height = "750px"),
  
  h2("Radial Plot", class = "mt-4"),
  plotOutput("radial_plot", height = "750px")
)

server <- function(input, output, session) {
  keep <- reactive(setdiff(grain_ids, as.integer(input$exclude)))
  
  upb <- reactive({
    validate(need(length(keep()) >= 2, "Keep at least two grains."))
    make_upb(keep())
  })
  
  output$data_table <- renderTable(
    cbind(Grain = grain_ids, C893),
    digits = 8
  )
  
  output$concordia_plot <- renderPlot({
    concordia(upb(), type = 2, show.age = 2)
  })
  
  output$radial_plot <- renderPlot({
    radialplot(upb(), k = "auto")
  })
}

shinyApp(ui, server)