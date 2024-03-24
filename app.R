### Least Shiny Path v0.23 ###
### 2024-03-24

library(shiny)
library(shinyjs)
library(sf)
library(terra)
library(leastcostpath)
library(shinyscreenshot)

create_conductance_matrix <- function(x) {
  cells <- which(!is.na(terra::values(x)))
  na_cells <- which(is.na(terra::values(x)))
  adj <- terra::adjacent(x = x, cells = cells, pairs = TRUE)
  adj <- adj[!adj[, 2] %in% na_cells, ]
  
  spatvals <- terra::values(x)[, 1]
  spatvals <- spatvals[adj[, 2]]
  
  ncells <- length(cells) + length(na_cells)
  cs_matrix <- Matrix::Matrix(data = 0, nrow = ncells, ncol = ncells)
  cs_matrix[adj] <- spatvals
  
  cs <- list(conductanceMatrix = cs_matrix, 
             costFunction = NA, 
             "maxSlope" = NA,
             exaggeration = FALSE,
             criticalSlope = NA,
             neighbours = NA, 
             nrow = terra::nrow(x), 
             ncol = terra::ncol(x), 
             "resolution" = terra::res(x), 
             "extent" = as.vector(terra::ext(x)),  
             crs = terra::crs(x, proj = TRUE))
  
  class(cs) <- "conductanceMatrix"
  
  return(cs)
}


ui <- fluidPage(
  tags$head(tags$style(HTML("
    .shiny-notification {
      font-size: 20px;  /* Adjust the font size */
      width: 300px;     /* Adjust the width */
      height: auto;     /* Adjust the height */
    }
    #shiny-notification-panel {
      top: 50%;         /* Position from the top */
      left: 50%;        /* Position from the left */
      transform: translate(-50%, -50%);  /* Center the notifications */
    }
  "))),
  sidebarLayout(
    sidebarPanel(
      tabsetPanel(
        tabPanel("Cost Surface",
          div( actionButton("about", "About")),
          div(
            h4("Cost Surfaces")
          ),
          div(
          actionButton("toggle1", "Use Your Own Cost Surface"),  
          conditionalPanel(
            condition = "input.toggle1 % 2 == 1",
            h3("Use Your Own Cost Surface"),
            fileInput("userCostSurfaceFile", "Upload User Cost Surface (GeoTIFF)"),
            radioButtons("resistOrConduct", "Cost surface type:",
                         c("Conductance" = "conduct",
                           "Resistance" = "resist"))
          )),
          div(
          actionButton("toggle2", "Generate New Cost Surface"),
          conditionalPanel(
            condition = "input.toggle2 % 2 == 1",
            h3("Calculate Cost Surface"),
            fileInput("demFile", "Upload DEM file (GeoTIFF)"),
            div(
              actionButton("useDemoDEM", "Use Demo DEM")
            ),
            selectInput("costFunction", "Select Cost Function", choices = c(
              #"basic slope", "enerscape log conductivity",
              "tobler", "tobler offpath", "davey", "rees", "irmischer-clarke male", 
              "irmischer-clarke offpath male", "irmischer-clarke female", "irmischer-clarke offpath female",
              "modified tobler", "garmy", "kondo-saino", "wheeled transport", "herzog", 
              "llobera-sluckin", "naismith", "minetti", "campbell", "campbell 2019 1",
              "campbell 2019 5", "campbell 2019 10", "campbell 2019 15", "campbell 2019 20", 
              "campbell 2019 25", "campbell 2019 30", "campbell 2019 35", "campbell 2019 40", 
              "campbell 2019 45", "campbell 2019 50", "campbell 2019 55", "campbell 2019 60", 
              "campbell 2019 65", "campbell 2019 70", "campbell 2019 75", "campbell 2019 80", 
              "campbell 2019 85", "campbell 2019 90", "campbell 2019 95", "campbell 2019 99", 
              "sullivan 167", "sullivan 5", "sullivan 833")),
            selectInput("neighbours", "Select Neighbours", choices = c(4, 8, 16, 32)),
            numericInput("critSlope", "Critical Slope", value = 12),
            numericInput("maxSlope", "Max Slope", value = NULL),
            checkboxInput("exaggeration", "Exaggeration", value = FALSE),
            div(
              actionButton("calculateCostSurface", "Calculate Cost Surface")
            ),
          )),
          
          div(
          h4("Waterways"),
          actionButton("addWaterways", "Load Waterways"),
          checkboxInput("showWaterways", "Show Waterways", value = FALSE)
          ),

          div(
            h4("Cost Surface Stochasticity"),
            actionButton("addStochasticity", "Add Global Stochasticity")
          ),
          div(
            h4("Download Raster"),
            downloadButton("downloadRaster2", "Download Displayed Raster"),
            actionButton("screenshot", "Save map image")
          )
        ), #close tabPanel
        tabPanel("Least cost path",
           div( actionButton("about2", "About")),
           div(
           fileInput("pointFile", "Upload Point File (CSV format)")
           ),
           div(
             actionButton("useDemoPts", "Use Demo Points")
           ),
           uiOutput("idMenu"),
           uiOutput("xMenu"),
           uiOutput("yMenu"),
           actionButton("setCoords", "Set Coordinates"),
           actionButton("showData", "Show Points Table"),
           actionButton("plotPoints", "Plot Points"),
           div(
           selectInput("lcpMethod", "Least cost path method", 
                       choices = c("Sequential", "From one to everywhere", "From everywhere to everywhere")
                       ),
           conditionalPanel(
             condition = "input.lcpMethod == 'Sequential'",
             uiOutput("sequenceFieldMenu")
           ),
           conditionalPanel(
             condition = "input.lcpMethod == 'From one to everywhere'",
             uiOutput("idValueMenu")
           ),
           selectInput("costSurface", "Select Cost Surface", 
                       choices = NULL),
           actionButton("calculatePaths", "Calculate Paths"),
           actionButton("plotPaths", "Plot Paths"),
           downloadButton("downloadPaths", "Download Paths")
           ),
           div(
            h3("Movement Corridors"),
            actionButton("calculateCorridors", "Calculate Path Corridors")
           ),
           div(
             h3("Path Density"),
             actionButton("createDensity", "Calculate Path Density"),
             downloadButton("downloadDensity", "Download Path Density Map")
           ),
           div(
             h3("Download Raster"),
             downloadButton("downloadRaster", "Download Displayed Raster"),
             actionButton("screenshot2", "Save map image")
           )
        ), #close tabPanel
        selectInput("rasterData", "Display Raster Map:", choices = NULL)
      ) #close tabsetPanel

    ), #close sidebarPanel
    
    mainPanel(
      fluidRow(
        id = "main_panel",
          plotOutput("rasterPlot", height = "100vh", width = "100%"),
          outputId = "rasterMap",
            height = "100vh" ,
            width = "100%"
      ) #close fluidRow
    ) #close mainPanel
  ) #close sidebarLayout
) #close UI fluidPage


server <- function(input, output, session) {
  options(shiny.maxRequestSize = 30*1024^2)
  # Initialize an empty list to store the raster data
  raster_data <- reactiveValues(data = list())
  values <- reactiveValues()
  cc <- list()
  
  
  observeEvent({
    c( input$about,  input$about2)
    }, {
    showModal(modalDialog(
      title = "About",
      tabsetPanel(
        tabPanel("About the App", 
                 p(HTML("<b>Welcome to the leastcostpath interface app</b>"),tags$br(),
                 "Nick Waber 2024",tags$br(),tags$br(),
                 "The functions in this app either use or are based primarily around Joseph Lewis's leastcostpath library (v2.0.12)", # and Emilio Berti's enerscape library (v1.1.0).
                   tags$br(),tags$br(),
                   "Citation:",
                   tags$br(),
                   "Lewis, J. (2023) leastcostpath: Modelling Pathways and Movement Potential Within a Landscape (version 2.0.12). 
Available at: https://github.com/josephlewis/leastcostpath"),tags$br(),
                 a("Link to leastcostpath", href="https://github.com/josephlewis/leastcostpath"),
                    tags$br(),tags$br(),
  #                 "Berti E (2024). _enerscape: Compute Energy Landscapes_. R package version 1.1.0,
  # <https://CRAN.R-project.org/package=enerscape>."
  ),
        
        tabPanel("Instructions", 
                 p(HTML("<b>IMPORTANT: all data must be in a projected CRS (i.e. a UTM grid).</b>"),
                   tags$br(),tags$br(),
                   "1. Get a cost surface.", tags$br(),
                   "In leastcostpath, cost surfaces are conductance matrices which quantify the ease of movement, rather than the difficulty, between grid squares. Higher values = easier travel.", tags$br(),
                   "If you have a cost surface raster, load it using the Upload User Cost Surface interface and procede to step 1.5 or 2.",
                   tags$br(),
                   "If you wish to create a cost surface, upload a DEM raster file using the Upload DEM interface:",
                   tags$br(),
                   "Select your preferred settings for cost surface creation.",tags$br(),
                   "Cost Function determines the calculation for the surface; different methods weight slope, elevation, distance, etc. differently.",tags$br(),
                   "Neighbours determines the search directions/distance for slope calculation for each DEM cell.",tags$br(),
                   "Critical Slope - from Lewis 2023: is  the transition where switchbacks become more effective than direct uphill or downhill paths . 
                   Cost of climbing the critical slope is twice as high as those for moving on flat terrain and is used for estimating the cost of using 
                   wheeled vehicles. Default value is 12, which is the postulated maximum gradient traversable by ancient transport (Verhagen and Jeneson, 
                   2012). Critical slope only used in  wheeled transport  cost function", tags$br(),
                   "Maximum Slope denotes the value above which a slope becomes an absolute obstacle.", tags$br(),
                   "Exaggeration (boolean) indicates whether uphill and downhill movement should be calculated seperately to make downhill less costly.", tags$br(),
                   tags$br(),
                   "1a Optional: Add Waterway Preference", tags$br(),
                   "Add a binary raster for navigable waterways (1 = water; 0 = non-water).  Weight indicates preference factor to be applied to waterways.", tags$br(),
                   tags$br(),
                   "1b Optional: Add Global Stochasticity", tags$br(),
                   "From Lewis:",tags$br(),
                   "The add_global_stochasticity to a conductanceMatrix is based on the method proposed by Pinto and Keitt (2009). 
                   Rather than using a static neighbourhood (for example as supplied in the neighbours function in the 
                   create_slope_cs), the neighbourhood is redefined such that the adjacency is non-deterministic and is instead 
                   determined randomly based on the threshold value", tags$br(),
                   "A lower percent quantile means less stochasticity will be incorporated.",tags$br(),
                   tags$br(),
                   "2. Add your points",tags$br(),
                   "Upload a CSV, selected the relevant fields, set the coordinates, and plot the points.",tags$br(),
                   tags$br(),
                   "3. Make least cost paths",tags$br(),
                   "Select a least cost path method.",tags$br(),
                   "Sequential follows your points according to the order in the ID field.",tags$br(),
                   "From One To Everywhere (FOTE) starts at a single point and then plots the path from that point directly to each other point.",tags$br(),
                   "From Everywhere To Everywhere (FETE/Pairwise) calculates the path between every pair of points",tags$br(),
                   "Paths can be downloaded as GPKG vector files using the interface.",tags$br(),
                   tags$br(),
                   "3.5a: Movement corridors",tags$br(),
                   "This creates rasters with the cumulative cell cost between the pairs of points. 
                   It offers a more nuanced understanding of terrain than a distinct line.",tags$br(),
                   "3.5b: Path density",tags$br(),
                   "This creates a raster file that quantifies the overlap between pathways.  
                   It does not display well on this app because the least cost path vector is being plotted 
                   directly over the path density rasters. Download the raster and load it in a GIS to see the density."
                   
                 )),

              tabPanel("Cost Surface Methods", 
                       p(HTML("<b>References for cost surface methods (from Lewis 2023):</b>"),
                         tags$br(),tags$br(),
                         
                         "Tobler, W. 1993. Three Presentations on Geographical Analysis and Modeling. Technical Report 93-1 (Santa Barbara, CA)",
                           tags$br(),tags$br(),
                           "Davey, R.C., M. Hayes and J.M. Norman 1994. “Running Uphill: An Experimental Result and Its Applications,” The Journal of the Operational Research Society 45, 25",
                           tags$br(),tags$br(),
                           "Rees, W.G. 2004. “Least-cost paths in mountainous terrain,” Computers &amp; Geosciences 30, 203–09",
                           tags$br(),tags$br(),
                           "Irmischer, I.J. and K.C. Clarke 2018. “Measuring and modeling the speed of human navigation,” Cartography and Geographic Information Science 45, 177–86",
                           tags$br(),tags$br(),
                           "Márquez-Pérez, J., I. Vallejo-Villalta and J.I. Álvarez-Francoso 2017. “Estimated travel time for walking trails in natural areas,” Geografisk Tidsskrift-Danish Journal of Geography 117, 53–62",
                           tags$br(),tags$br(),
                           "Garmy, P. et al. 2005. “Logiques spatiales et ‘systèmes de villes’ en Lodévois de l’Antiquité à la période moderne,” Temps et espaces de l’homme en société, analyses et modèles spatiaux en archéologie 335–46",
                           tags$br(),tags$br(),
                           "Kondo, Y. and Y. Seino 2010. “GPS-aided walking experiments and data-driven travel cost modeling on the historical road of Nakasendō-Kisoji (Central Highland Japan),” Making History Interactive (Proceedings of the 37th International Conference, Williamsburg, Virginia, United States of America) 158–65",
                           tags$br(),tags$br(),
                           "Herzog, I. 2013. “The potential and limits of Optimal Path Analysis,” in Bevan, A. and M. Lake (edd.), Computational approaches to archaeological spaces (Publications of the Institute of Archaeology, University College London) 179–211",
                           tags$br(),tags$br(),
                           "Llobera, M. and T.J. Sluckin 2007. “Zigzagging: Theoretical insights on climbing strategies,” Journal of Theoretical Biology 249, 206–17",
                           tags$br(),tags$br(),
                           "Naismith, W. 1892. “Excursions: Cruach Ardran, Stobinian, and Ben More,” Scottish Mountaineering club journal 2, 136",
                           tags$br(),tags$br(),
                           "Minetti, A.E. et al. 2002. “Energy cost of walking and running at extreme uphill and downhill slopes,” Journal of Applied Physiology 93, 1039–46",
                           tags$br(),tags$br(),
                           "Campbell, M.J., P.E. Dennison and B.W. Butler 2017. “A LiDAR-based analysis of the effects of slope, vegetation density, and ground surface roughness on travel rates for wildland firefighter escape route mapping,” Int. J. Wildland Fire 26, 884",
                           tags$br(),tags$br(),
                           "Campbell, M.J. et al. 2019. “Using crowdsourced fitness tracker data to model the relationship between slope and travel rates,” Applied Geography 106, 93–107",
                           tags$br(),tags$br(),
                           "Sullivan, P.R. et al. 2020. “Modeling Wildland Firefighter Travel Rates by Terrain Slope: Results from GPS-Tracking of Type 1 Crew Movement,” Fire 3, 52"
                       )
                      )
        

    ),      
    size = "l",  # Large modal window
    easyClose = T)
    )
  })

  observeEvent(input$toggle1, {
    shinyjs::toggle(id = "section1")  # Toggle visibility of section 1
  })
  
  observeEvent(input$toggle2, {
    shinyjs::toggle(id = "section2")  # Toggle visibility of section 2
  })
  
  observeEvent(input$useDemoDEM, {
    raster_data$data$DEM <- terra::rast("PNW_clipped_sm.tif")
    updateSelectInput(session, "rasterData", choices = names(raster_data$data))
  })
  
  user_cost_surface <- reactive({
    req(input$userCostSurfaceFile)
    inFile <- input$userCostSurfaceFile
    if (is.null(inFile)) return(NULL)

    # Load the raster file
    ucs <- terra::rast(inFile$datapath)

    # Check if the cost surface type is "resist"
    if (input$resistOrConduct == "resist") {
      # Convert the resistance map to a conductance map
      ucs <- 1 / ucs
    }

    return(ucs)
  })

  observe({
    ucs <- user_cost_surface()
    if (!is.null(ucs)) {
      raster_data$data[["User Cost Surface"]] <- ucs
      updateSelectInput(session, "rasterData", choices = names(raster_data$data))
    }
  })
  
  uploaded_dem <- reactive({
    req(input$demFile)
    inFile <- input$demFile
    if (is.null(inFile)) return(NULL)
    
    # Load the DEM
    dem <- terra::rast(inFile$datapath)
    
    # Generate a hillshade from the DEM
    slope <- terra::terrain(dem, "slope", unit = "radians")
    aspect <- terra::terrain(dem, "aspect", unit = "radians")
    hillshade <- terra::shade(slope, aspect, angle = 35, direction = 315)
    
    return(list(dem = dem, hillshade = hillshade))
  })
  
  # When a new DEM file is uploaded, add it to the list and update the selectInput choices
  observe({
    dem_data <- uploaded_dem()
    if (!is.null(dem_data)) {
      raster_data$data$DEM <- dem_data$dem
      raster_data$data$Hillshade <- dem_data$hillshade
      updateSelectInput(session, "rasterData", choices = names(raster_data$data))
    }
  })
  
  waterways <- reactive({
    req(input$waterwaysFile)
    inFile <- input$waterwaysFile
    if (is.null(inFile)) return(NULL)
    
    # Load the waterways layer
    waterways <- terra::rast(inFile$datapath)
    return(waterways)
  })
  
  
  # Reactive variable for points data
  points_data <- reactiveVal()
  
  # Load demo data when the button is pressed
  observeEvent(input$useDemoPts, {
    points_data(read.csv("Trade Centres DEMO CONIC.csv", stringsAsFactors = F))
  })
  
  # Load data from uploaded file
  observeEvent(input$pointFile, {
    inFile <- input$pointFile
    if (!is.null(inFile)) {
      points_data(read.csv(inFile$datapath, stringsAsFactors = FALSE))
    }
  })
  
  output$idMenu <- renderUI({
    req(points_data())
    selectInput("idField", "ID", choices = names(points_data()))
  })
  
  output$xMenu <- renderUI({
    req(points_data())
    selectInput("xField", "X Coordinate Field", choices = names(points_data()))
  })
  
  output$yMenu <- renderUI({
    req(points_data())
    selectInput("yField", "Y Coordinate Field", choices = names(points_data()))
  })
  
  observeEvent(input$setCoords, {
    req(input$idField, input$xField, input$yField, points_data())
    
    points <- points_data()
    
    points_sf <- sf::st_sf(
      geometry = sf::st_sfc(
        lapply(1:nrow(points), function(i) {
          sf::st_point(c(points[i, input$xField], points[i, input$yField]))
        })
      ),
      crs = terra::crs(raster_data$data[[input$rasterData]])
    )
    
    # Store the spatial points in a reactiveValues object for use elsewhere in the app
    values$points_sf <- points_sf
  })
  
  
  observeEvent(input$showData, {
    showModal(modalDialog(
      title = "Points Data",
      tableOutput("pointsTable")
    ))
  })
  output$pointsTable <- renderTable({
    points_data()
  })
  
  
  

  
  selected_cost_function <- reactive({
    cost_function_name <- input$costFunction
    cost_function(cost_function_name)
  })
  
  output$rasterPlot <- renderPlot({
    dem <- uploaded_dem()
    if (!is.null(dem)) {
      terra::plot(dem, main = "Uploaded DEM")
      if (input$showWaterways && !is.null(waterways())) {
        cols <- c("0" = "transparent", "1" = "white")
        terra::plot(waterways(), add = TRUE, col = cols)
      }
    }
  })
  
  
  observeEvent(input$calculateCostSurface, {
    cost_function <- input$costFunction
    neighbours <- as.numeric(input$neighbours)
    crit_slope <- input$critSlope
    max_slope <- input$maxSlope
    exaggeration <- input$exaggeration
    
    
    if(input$costFunction == "basic slope"){
      slp <- terrain(raster_data$data$DEM, "slope", unit = "degrees")
      slp <- tan( pi/180 * slp)*100
      max_val <- terra::minmax(slp)
      slp_inverted <- max_val[2] - slp
      slope_cs <- create_cs(
        x = slp_inverted
        # x = raster_data$data$DEM,
        # neighbours = neighbours,
        # exaggeration = exaggeration,
        # max_slope = max_slope
        )
      slope_cs_rast <- leastcostpath::rasterise(slope_cs)
    }
    if(input$costFunction == "enerscape log conductivity"){
      en <- enerscape(raster_data$data$DEM, 150/2.2, neigh = neighbours, unit = "kcal")
      en_log <- log(en)
      max_val <- terra::minmax(en_log)
      en_log_inverted <- max_val[2] - en_log
      slope_cs <- create_cs(
        x = en_log_inverted,
        # x = raster_data$data$DEM,
        # neighbours = neighbours,
        # exaggeration = exaggeration,
        # max_slope = max_slope
      )
      slope_cs_rast <- leastcostpath::rasterise(slope_cs)
    }
    else {
      # Calculate the cost surface
      slope_cs <- create_slope_cs(
        x = raster_data$data$DEM, 
        cost_function = cost_function, 
        neighbours = neighbours,
        crit_slope = crit_slope,
        max_slope = max_slope,
        exaggeration = exaggeration)
      slope_cs_rast <- leastcostpath::rasterise(slope_cs)
    }
    
    observe({
      if (!is.null(slope_cs)) {
        showNotification("Cost surface calculated successfully.")
      }
    })
    
    if (!is.null(slope_cs_rast)) {
      raster_data$data[["Cost Surface"]] <- slope_cs_rast
      updateSelectInput(session, "rasterData", choices = names(raster_data$data))
    }
  })

  
  # Observe the raster_data$data list for changes
  observe({
    # Get the names of the rasters
    raster_names <- names(raster_data$data)
    
    # Filter the names to include only those that contain the word "Cost"
    cost_raster_names <- raster_names[grep("Cost", raster_names)]
    
    # Update the selectInput choices
    updateSelectInput(session, "costSurface", choices = cost_raster_names)
    updateSelectInput(session, "costSurfaceSource", choices = cost_raster_names)
  })
  
  observeEvent(input$addWaterways, {
    showModal(modalDialog(
      title = "Add Waterways",
      fileInput("waterwaysFile", "Upload Waterways Layer (GeoTIFF)"),
      uiOutput("costSurfaceSource"),
      #      selectInput("costSurfaceSource", "Apply waterway preference to surface:", choices = NULL),
      checkboxInput("preferWaterways", "Prefer Waterways", value = FALSE),
      numericInput("waterWeight", "Preference Weight", value = 2.5),
      numericInput("bufferDistance", "Waterway preference buffer distance (m)", value = 5000),
      actionButton("adjustWaterways", "Adjust for Waterways"),
      actionButton("adjustWaterwaysBuffer", "Adjust for waterways with buffer")
    ))
  })
  
  
  output$costSurfaceSource <- renderUI({
    selectInput("costSurfaceSource", "Apply waterway preference to surface:", choices = names(raster_data$data))
  })

  observeEvent(input$adjustWaterways, {
    req(input$waterwaysFile, input$costSurfaceSource)

    # Load the waterways layer
    waterways <- terra::rast(input$waterwaysFile$datapath)

    # Determine which cost surface to use
    if (input$costSurfaceSource == "Cost Surface") {
      cs <- raster_data$data[["Cost Surface"]]
    } else if (input$costSurfaceSource == "User Cost Surface") {
      cs <- raster_data$data[["User Cost Surface"]]
    }

    # Define a function to adjust the cost surface based on the waterway preference
    adjust_func <- function(x, y) {
      # If the value is greater than 0 (indicating a waterway), multiply by the weight
      # If the value is 0 (indicating non-waterway), multiply by 1
      return(terra::ifel(y > 0, x * input$waterWeight, x))
    }

    # Use the adjust_func function to adjust the cost surface
    cs_rast_waterways <- adjust_func(cs, waterways)

    # Add the adjusted cost surface to the raster_data$data list
    raster_data$data[[paste(input$costSurfaceSource, "with Waterways")]] <- cs_rast_waterways

    # Update the selectInput choices
    updateSelectInput(session, "rasterData", choices = names(raster_data$data))
    showNotification("Waterway cost adjustment calculated successfully.")
  })
  
  observeEvent(input$adjustWaterwaysBuffer, {
    req(input$waterwaysFile, input$costSurfaceSource, input$bufferDistance)
    
    # Load the waterways layer
    waterways <- terra::rast(input$waterwaysFile$datapath)
    waterways[waterways == 0] <- NA
    
    # Calculate the proximity grid
    proximity <- terra::distance(waterways, unit = "m")
    
    # Adjust the preference values based on the proximity
    preference <- terra::ifel(proximity <= input$bufferDistance, 
                              input$waterWeight * (1 - proximity / input$bufferDistance) + 1, 
                              1)
    
    # Determine which cost surface to use
    if (input$costSurfaceSource == "Cost Surface") {
      cs <- raster_data$data[["Cost Surface"]]
    } else if (input$costSurfaceSource == "User Cost Surface") {
      cs <- raster_data$data[["User Cost Surface"]]
    }
    # Multiply the cost surface by the preference weight grid
    cs_rast_waterways_buffer <- cs * preference
    
    
    # Add the adjusted cost surface to the raster_data$data list
    raster_data$data[[paste(input$costSurfaceSource, "with Waterways Buffer")]] <- cs_rast_waterways_buffer
    
    # Update the selectInput choices
    updateSelectInput(session, "rasterData", choices = names(raster_data$data))
    
    showNotification("Waterway buffer cost adjustment calculated successfully.")
  })
  
  
  observeEvent(input$addStochasticity, {
    showModal(modalDialog(
      title = "Add Global Stochasticity",
      numericInput("percentQuantile", "Percent Quantile", value = 0.5, min = 0.00, max = 1.00),
      selectInput("costSurfaceSource", "Select Cost Surface Source", choices = c("Cost Surface", "User Cost Surface")),
      actionButton("goButton", "Go")
    ))
  })
  observeEvent(input$goButton, {
    req(input$percentQuantile)
    
    # Determine which cost surface to use
    if (input$costSurfaceSource == "Cost Surface") {
      cs <- create_conductance_matrix(raster_data$data[["Cost Surface"]])
    } else if (input$costSurfaceSource == "User Cost Surface") {
      cs <- create_conductance_matrix(raster_data$data[["User Cost Surface"]])
      cs$crs <- terra::crs(raster_data$data[["User Cost Surface"]])
    }
    
    stochastic_cs <- leastcostpath::add_global_stochasticity(x = cs, percent_quantile = input$percentQuantile)
    raster_data$data[["Stochastic Cost Surface"]] <- rasterise(stochastic_cs)
    #updateSelectInput(session, "rasterData", choices = names(raster_data$data))
    
    removeModal()
  })
  
  
  output$rasterPlot <- renderPlot({
    selected_data <- input$rasterData
    if (!is.null(selected_data)) {
      if (selected_data == "Hillshade") {
        terra::plot(raster_data$data[[selected_data]], 
                    col = grey(0:100/100), 
                    main = selected_data,
                    legend = FALSE)
      } else {
        terra::plot(raster_data$data[[selected_data]], main = selected_data)
      }
      if (input$showWaterways && !is.null(waterways())) {
        cols <- c("0" = "transparent", "1" = "white")
        terra::plot(waterways(), add = TRUE, col = cols)
      }
    }
  })

  
  observeEvent(input$plotPoints, {
    output$rasterPlot <- renderPlot({
      selected_data <- input$rasterData
      if (!is.null(selected_data)) {
        if (selected_data == "Hillshade") {
          terra::plot(raster_data$data[[selected_data]], 
                      col = grey(0:100/100), 
                      main = selected_data,
                      legend = FALSE)
        } else {
          terra::plot(raster_data$data[[selected_data]], main = selected_data)
        }
        if (input$showWaterways && !is.null(waterways())) {
          cols <- c("0" = "transparent", "1" = "white")
          terra::plot(waterways(), add = TRUE, col = cols)
        }
        if (!is.null(values$points_sf)) {
          plot(values$points_sf,
               add = TRUE,
               pch = 16)
          points <- points_data()
          coords <- sf::st_coordinates(values$points_sf)
          text(coords[,1], coords[,2], labels = points[[input$idField]], pos = 4, cex = 0.8)
        }

      }
    })
  })
  
  observeEvent(input$submitEPSG, {
    req(input$epsgCode)
    selected_data <- input$rasterData
    rast <- raster_data$data[[selected_data]]
    terra::crs(rast) <- as.character(input$epsgCode)
    removeModal()
  })
  
  
  output$sequenceFieldMenu <- renderUI({
    req(points_data())
    selectInput("sequenceField", "Sequence field", choices = names(points_data()))
  })
  
  output$idValueMenu <- renderUI({
    req(points_data(), input$idField)
    selectInput("idValue", "Origin", choices = unique(points_data()[[input$idField]]))
  })
  
  observeEvent(input$calculatePaths, {
    req(values$points_sf, input$lcpMethod)
    
    if (input$costSurface == "Cost Surface") {
      cs <- create_conductance_matrix(raster_data$data[["Cost Surface"]])
    } else if (input$costSurface == "User Cost Surface") {
      cs <- create_conductance_matrix(raster_data$data[["User Cost Surface"]])
      cs$crs <- terra::crs(raster_data$data[["User Cost Surface"]])
    } else if (input$costSurface == "Stochastic Cost Surface") {
      cs <- create_conductance_matrix(raster_data$data[["Stochastic Cost Surface"]])
      cs$crs <- terra::crs(raster_data$data[["Stochastic Cost Surface"]])
    } else if (input$costSurface == "Cost Surface with Waterways") {
      cs <- create_conductance_matrix(raster_data$data[["Cost Surface with Waterways"]])
      cs$crs <- terra::crs(raster_data$data[["Cost Surface with Waterways"]])
    } else if (input$costSurface == "User Cost Surface with Waterways") {
      cs <- create_conductance_matrix(raster_data$data[["User Cost Surface with Waterways"]])
      cs$crs <- terra::crs(raster_data$data[["User Cost Surface with Waterways"]])
    } else if (input$costSurface == "Cost Surface with Waterways Buffer") {
      cs <- create_conductance_matrix(raster_data$data[["Cost Surface with Waterways Buffer"]])
      cs$crs <- terra::crs(raster_data$data[["Cost Surface with Waterways Buffer"]])
    } else if (input$costSurface == "User Cost Surface with Waterways Buffer") {
      cs <- create_conductance_matrix(raster_data$data[["User Cost Surface with Waterways Buffer"]])
      cs$crs <- terra::crs(raster_data$data[["User Cost Surface with Waterways Buffer"]])
    }
    
    if (input$lcpMethod == "Sequential") {
      points <- points_data()[order(points_data()[[input$sequenceField]]), ]
      #points_data()[[input$idField]] <- as.numeric(points_data()[[input$idField]])
      values$points_sf <- values$points_sf[order(points_data()[[input$sequenceField]]), ]
      
      lcps <- lapply(1:(nrow(points) - 1), function(i) {
        create_lcp(x = cs, origin = values$points_sf[i,], destination = values$points_sf[i + 1,])
      })
      
      lcps <- do.call(rbind, lcps)       # Combine all the least cost paths into one object
      
    } else if (input$lcpMethod == "From one to everywhere") {
      origin <- values$points_sf[which(points_data()[[input$idField]] == input$idValue), ]
      lcps <- create_lcp(x = cs, origin = origin, destination = values$points_sf)
      
    } else if (input$lcpMethod == "From everywhere to everywhere") {
      lcps <- create_FETE_lcps(x = cs, locations = values$points_sf)
    }
    
    values$lcps <- lcps
    
    observe(
      if(!is.null(lcps)){
        showNotification("Least cost paths calculated successfully.")
      })
  })
  
  observeEvent(input$plotPaths, {
    output$rasterPlot <- renderPlot({
      selected_data <- input$rasterData
      if (!is.null(selected_data)) {
        if (selected_data == "Hillshade") {
          terra::plot(raster_data$data[[selected_data]], 
                      col = grey(0:100/100), 
                      main = selected_data,
                      legend = FALSE)
        } else {
          terra::plot(raster_data$data[[selected_data]], main = selected_data)
        }
        if (input$showWaterways && !is.null(waterways())) {
          cols <- c("0" = "transparent", "1" = "white")
          terra::plot(waterways(), add = TRUE, col = cols)
        }
        if (!is.null(values$points_sf)) {
          plot(values$points_sf, 
               add = TRUE,
               pch = 16)
          coords <- sf::st_coordinates(values$points_sf)
          text(coords[,1], coords[,2], labels = points_data()[[input$idField]], pos = 4, cex = 0.8)
        }
        if (!is.null(values$lcps)) {
          plot(terra::vect(values$lcps), add = TRUE)
        }
      }
    })
  })
  
  output$downloadPaths <- downloadHandler(
    filename = function() {
      paste("paths_", input$lcpMethod, "_", input$costFunction, "_", Sys.Date(), ".gpkg", sep="")
    },
    content = function(file) {
      sf::st_write(values$lcps, file)
    }
  )
  
  observeEvent(input$createDensity, {
    req(values$lcps)
    
    # Calculate the path density
    selected_data <- input$rasterData
    density <- leastcostpath::create_lcp_density(raster_data$data[[selected_data]], values$lcps)
    raster_data$data[["Path density map"]] <- density
    updateSelectInput(session, "rasterData", choices = names(raster_data$data))
    
    observe(
      if(!is.null(density)){
        showNotification("Path density calculated successfully.")
      })
  })
  
  observeEvent(input$calculateCorridors, {
    
    # Check the selected least cost path method
    if (input$lcpMethod == "Sequential") {
      
      # Loop over the points sequentially
      for(i in 1:(nrow(values$points_sf) - 1)) {
        # Create a cost corridor for each pair of points and store it in the list
        cc[[i]] <- create_cost_corridor(slope_cs, values$points_sf[i,], values$points_sf[i + 1,])
      }
      
    } else if (input$lcpMethod == "From one to everywhere") {
      
      # Find the origin point
      origin <- values$points_sf[which(points_data()[[input$idField]] == input$idValue), ]
      
      # Loop over the other points
      for(i in 1:nrow(values$points_sf)) {
        if (i != which(points_data()[[input$idField]] == input$idValue)) {
          # Create a cost corridor from the origin to each other point and store it in the list
          cc[[i]] <- create_cost_corridor(slope_cs, origin, values$points_sf[i,])
        }
      }
      
    } else if (input$lcpMethod == "From everywhere to everywhere") {
      
      # Initialize a counter
      counter <- 1
      
      # Loop over all pairs of points
      for(i in 1:(nrow(values$points_sf) - 1)) {
        for(j in (i + 1):nrow(values$points_sf)) {
          # Create a cost corridor for each pair of points and store it in the list
          cc[[counter]] <- create_cost_corridor(slope_cs, values$points_sf[i,], values$points_sf[j,])
          counter <- counter + 1
        }
      }
      
    }
    
    # Store each cost corridor as a separate raster in your raster_data$data list and update the selectInput choices
    for(i in seq_along(cc)) {
      raster_data$data[[paste0("Cost Corridor ", i)]] <- cc[[i]]
    }
    updateSelectInput(session, "rasterData", choices = names(raster_data$data))

    observe(
      if(!is.null(cc)){
        showNotification("Cost corridors calculated successfully.")
      })
  })
  

  output$downloadRaster <- downloadHandler(
    filename = function() {
      paste(input$rasterData,"_", gsub(" ", "-", input$costFunction), "_", Sys.Date(), ".tif", sep="")
    },
    content = function(file) {
      # Get the currently selected raster data
      selected_data <- input$rasterData
        terra::writeRaster(
          raster_data$data[[selected_data]],
          file)
    }
  )

  output$downloadRaster <- downloadHandler(
    filename = function() {
      if (input$rasterData == "Cost Surface") {
        paste(input$rasterData, "_", gsub(" ", "-", input$costFunction), "_", Sys.Date(), ".tif", sep="")
      } else {
        paste(input$rasterData, "_", Sys.Date(), ".tif", sep="")
      }
    },
    content = function(file) {
      # Get the currently selected raster data
      selected_data <- input$rasterData
      terra::writeRaster(
        raster_data$data[[selected_data]],
        file)
    }
  )
  
  output$downloadRaster2 <- downloadHandler(
    filename = function() {
      if (input$rasterData == "Cost Surface") {
        paste(input$rasterData, "_", gsub(" ", "-", input$costFunction), "_", Sys.Date(), ".tif", sep="")
      } else {
        paste(input$rasterData, "_", Sys.Date(), ".tif", sep="")
      }
    },
    content = function(file) {
      # Get the currently selected raster data
      selected_data <- input$rasterData
      terra::writeRaster(
        raster_data$data[[selected_data]],
        file)
    }
  ) 
  
  observeEvent(input$screenshot, {
    filename <- if (input$rasterData == "Cost Surface") {
      paste(input$rasterData, "_", gsub(" ", "-", input$costFunction), "_", Sys.Date(), sep="")
    } else {
      paste(input$rasterData, "_", Sys.Date(), sep="")
    }
    screenshot(selector = "#main_panel", filename = filename)
    
  })
  
  observeEvent(input$screenshot2, {
    filename <- if (input$rasterData == "Cost Surface") {
        paste(input$rasterData, "_", gsub(" ", "-", input$costFunction), "_", Sys.Date(), sep="")
      } else {
        paste(input$rasterData, "_", Sys.Date(), sep="")
      }
    screenshot(selector = "#main_panel", filename = filename)
    
  })
  
    # ##debugger
  # observe({
  #   # Get the names of the rasters
  #   raster_names <- names(raster_data$data)
  #   
  #   # If there are any rasters
  #   if (length(raster_names) > 0) {
  #     # Print the names and classes of the rasters
  #     print("Names and classes of uploaded rasters:")
  #     for (raster_name in raster_names) {
  #       rast <- raster_data$data[[raster_name]]
  #       print(paste0(raster_name, " (", class(rast), ")"))
  #     }
  #   }
  #   
  #   # Print file upload status
  #   if (!is.null(input$userCostSurfaceFile)) {
  #     print(paste("User Cost Surface file uploaded: ", input$userCostSurfaceFile$name))
  #   }
  #   
  #   if (!is.null(input$demFile)) {
  #     print(paste("DEM file uploaded: ", input$demFile$name))
  #   }
  #   
  #   if (!is.null(input$pointFile)) {
  #     print(paste("Point file uploaded: ", input$pointFile$name))
  #   }
  # })
  
} # close server

shinyApp(ui, server)
