###### script to create a shiny app with the database of Arteaga et al. (2016) ######

## Load library
library(shiny)
library(ggplot2)
library(dplyr)
library(raster)
library(shinydashboard)
library(plotly)
library(leaflet)
library(SNPRelate)
library(tidyr)


####### the user can download the database of eigenvalues, coordinates and 19 WorldClim variables #########
# load the gds file
genofile <- snpgdsOpen("../data/maicesArtegaetal2015.gds", allow.duplicate = TRUE)

# get the names of the samples from the gds file
sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))

##### Metadata #####
# load the database with maize information
fullmat <- read.delim(file= "../meta/maizteocintle_SNP50k_meta_extended.txt")

#get samples names
samples_names <- read.delim("../data/maicesArtegaetal2015.fam",
                            header=FALSE, sep=" ") %>% select(., V2) %>% rename(., Indiv=V2) 

### run PCA with the gds file ###
pca <- snpgdsPCA(genofile, num.thread=2)

# calculate the percent variation in the components
pc.percent <- pca$varprop*100

# round number to two digits
x<-round(pc.percent, 2)

# select the first six components and convert them into a dataframe
eigen_corn <- as.data.frame(pca$eigenvect[,1:6]) 

# change column names
colnames(eigen_corn) <- c("EV1", "EV2", "EV3", "EV4", "EV5", "EV6")

# select the variables latitude and longitude  
coord_corn <- dplyr::select(fullmat, Longitud, Latitud)

# get worldclim variables
r <- getData("worldclim",var="bio",res=2.5, path = "../data/")

# create object from coordinates
points_corn <- SpatialPoints(coord_corn, proj4string = r@crs)

# extract the bioclim values for each species
values_corn <- raster::extract(r, points_corn)

# convert to data frame
wc_corn <- as.data.frame(values_corn)

# Loop to divide temperature data between 10 
x <- c("bio1", "bio5", "bio6", "bio7", "bio8", "bio9", "bio10", "bio11")
for (var in x) {
  wc_corn[var] <- wc_corn[var]/10
}

# bind datafremes by columthe `app.R` script generates the shiny dashboard apps
tab_corn <- cbind(samples_names, fullmat, eigen_corn, wc_corn) 

#### Loops for shiny, create the labels of the x and y axes of the PCA ####
eigen_pca <- vector("character")
eigen_names <- vector("character")
for (i in 1:6){
  eigen_pca[i] <- paste("Eigenvector", i, "explaining", round(pc.percent, 2)[i], "%")
  eigen_names[i] <- paste0("EV", i)
}
attr(eigen_names, "names") <- eigen_pca

#### Admixture ###
## get Q files 
x <- list.files(paste0("../data/admixture/")) 
Qfile<-x[grep("*.Q", x) ]

# Loop to create a list with every K
Qval <- list()
for (j in 1:length(Qfile)){
  # read Q file
  Qval[[j]] <-read.table(paste0("../data/admixture/", Qfile[j]))
  names(Qval[[j]])<-paste0("K", 1:ncol(Qval[[j]]))
}

# loop to extract the elements of the list as dataframe
for (i in 1:length(Qval)) {
  assign(paste0("Qval", i), as.data.frame(Qval[[i]]))
}

# combine elements of data frames to create a df for each K
Qval1 <- cbind(INDIV=samples_names$Indiv, bio1 = tab_corn$bio1, bio12 = tab_corn$bio12, bio3 = tab_corn$bio3, Altitud = tab_corn$Altitud, Qval1)
Qval2 <- cbind(INDIV=samples_names$Indiv, bio1 = tab_corn$bio1, bio12 = tab_corn$bio12, bio3 = tab_corn$bio3, Altitud = tab_corn$Altitud, Qval2)
Qval3 <- cbind(INDIV=samples_names$Indiv, bio1 = tab_corn$bio1, bio12 = tab_corn$bio12, bio3 = tab_corn$bio3, Altitud = tab_corn$Altitud, Qval3)
Qval4 <- cbind(INDIV=samples_names$Indiv, bio1 = tab_corn$bio1, bio12 = tab_corn$bio12, bio3 = tab_corn$bio3, Altitud = tab_corn$Altitud, Qval4)
Qval5 <- cbind(INDIV=samples_names$Indiv, bio1 = tab_corn$bio1, bio12 = tab_corn$bio12, bio3 = tab_corn$bio3, Altitud = tab_corn$Altitud, Qval5)

# transform to long format  
Qval1_long<- gather(Qval1, key=Kgroup, value=Qadmixture, 6:ncol(Qval1))
Qval2_long<- gather(Qval2, key=Kgroup, value=Qadmixture, 6:ncol(Qval2))
Qval3_long<- gather(Qval3, key=Kgroup, value=Qadmixture, 6:ncol(Qval3))
Qval4_long<- gather(Qval4, key=Kgroup, value=Qadmixture, 6:ncol(Qval4))
Qval5_long<- gather(Qval5, key=Kgroup, value=Qadmixture, 6:ncol(Qval5))


#### Shiny App #####
## The user interface (ui) has a title panel and a sidebar layout, which includes a sidebar panel and the main panel

header <- dashboardHeader(title = "Maize Analysis") # title of shiny dashboard

sidebar <- dashboardSidebar(
  sidebarMenu(  # PCA section, chart title = "PCA on SNP genotypes
    menuItem("PCA on SNP genotypes", icon = icon("bar-chart-o"), # bar-chart-o is an icon used to represented chart
             # Input directly under menuItem                      
             selectInput("x_axis", "Choose the eigenvector for the x-axis:", # Create a select list of input control for the x-axis
                         choices = c("Eigenvector 1" = "EV1",                
                                     "Eigenvector 2" = "EV2",
                                     "Eigenvector 3" = "EV3",
                                     "Eigenvector 4" = "EV4",
                                     "Eigenvector 5" = "EV5",
                                     "Eigenvector 6" = "EV6"), 
                         selected = "EV2"),
             
             selectInput("y_axis", "Choose the eigenvector for the y-axis:", # Create a select list of input control for the y-axis
                         choices = c("Eigenvector 1" = "EV1",
                                     "Eigenvector 2" = "EV2",
                                     "Eigenvector 3" = "EV3",
                                     "Eigenvector 4" = "EV4",
                                     "Eigenvector 5" = "EV5",
                                     "Eigenvector 6" = "EV6"),
                         selected = "EV1"), 
             
             selectInput("coloring", "Define the point colours:", ## Create a select list of variables to color the points on the PCA
                         choices = c("Annual Mean Temperature" = "bio1",
                                     "Annual Precipitation" = "bio12",
                                     "Isothermality" = "bio3", 
                                     "Altitude" = "Altitud"),
                         selected = "Altitud") 
    ),
    
    menuItem("Maize distribution map", icon = icon("map")), ## Map section
    
    menuItem("Metadata", icon = icon("table"),             ## Table section
             selectInput("dataset", "Choose a dataset:", 
                         choices = c("Eigenvalues", "Coordinates", "Worlclime")), ## Create a select list of variables of dataframes to download 
             radioButtons("filetype", "File type:",
                          choices = c("csv", "tsv")),    ## Create options of format to download dataframe
             downloadButton('downloadData', 'Download')), ## Create the button download
  
      menuItem("Admixture", icon = icon("bar-chart-o"),  ## Admixture secction
             selectInput("Kvalue", "Choose the value for K:", ## select a K 
                         choices = c("K1" = "Qval1_long", 
                                     "K2" = "Qval2_long",
                                     "K3" = "Qval3_long",
                                     "K4" = "Qval4_long",
                                     "K5" = "Qval5_long"),
                         selected = "Qval3_long"))
    
  )
)

## The main body contains the PCA, Map, and Table. This is the part to output all the sections.

body <- dashboardBody(
  fluidRow(
    box(title = "PCA", status = "primary", solidHeader = TRUE, plotlyOutput("pcaplot")),
    box(title = "MAP", status = "primary", solidHeader = TRUE, leafletOutput("mymap"),
        absolutePanel(top = 10, right = 10)),
        
    box(title = "Table", status = "primary", solidHeader = TRUE, tableOutput('table'), style = "height:400px; overflow-y: scroll;overflow-x: scroll;"),
    box(title = "Structure", status = "primary", solidHeader = TRUE, plotOutput("structure"))
  )
  
)

## We will create the Server part where the program and logic behind valueBoxOutput()
## and plotOutput() are added with renderValueBox() and renderPlot() respectively.
## These are enclosed inside a server function , with input and output as its parameters.

server = function(input, output, session) { 
  
  ##### PCA #####
  
  ## create a render plot, this is that the graph will change each time you select the variables defined above
  
  output$pcaplot <- renderPlotly({
    
    # build graph with ggplot syntax
    p <- ggplot(tab_corn, aes_string(x = input$x_axis, y = input$y_axis, color = input$coloring)) +  # Defined variables above in PCA secction
      geom_point() + labs(x = names(eigen_names[which(eigen_names == input$x_axis)]),   # input$x_axi = PCA eigenvalue selected for x axis 
                          y = names(eigen_names[which(eigen_names == input$y_axis)]),
                          colour = "") +  # input$y_axis = PCA eigenvalue selected for y axis
      scale_color_gradient(low = "#35b7c4", high = "#4455a0") # input$coloring = color points according to variable selected
    
    ggplotly(p)   
    
    
  })
  
  ####### MAP #######
  
  ## Color palette defined according Altitude variable
  ## This part has to change according to variable selected to color the points in PCA
  
  output$mymap <- renderLeaflet({
    # Use leaflet() here, and only include aspects of the map that
    # won't need to change dynamically (at least, not unless the
    # entire map is being torn down and recreated).
    leaflet(tab_corn) %>% addTiles() %>%
      fitBounds(~min(Longitud), ~min(Latitud), ~max(Longitud), ~max(Latitud))
  })
  
  # Incremental changes to the map (in this case, replacing the circles when a new color is chosen) should be performed in an observeEvent.
  # can be used to trigger a piece of code when a certain event occurs.
  observeEvent (input$coloring, {  # Act on input$coloring
    if(input$coloring == "bio1"){  # If bio1 is selected in input $ coloring the following is executed
      pal <- colorNumeric(
        palette = colorRampPalette(c("#35b7c4", "#4455a0"))(length(tab_corn$bio1)), 
        domain = tab_corn$bio1) # The Color palette will be created from bio1
      leafletProxy("mymap", data = tab_corn) %>%
        clearShapes() %>%
        addCircleMarkers(lng = ~Longitud, lat = ~Latitud, color = ~pal(bio1),  fillColor = ~pal(bio1), # The colors of the circular markers 
                         fillOpacity = 1, popup = ~paste(bio1), radius = 4)  # are defined from the palette created above and bio1
    }
    else if(input$coloring == "bio3"){  # to color from another variable we use else if and we repeat the same process as above for each variable
      pal <- colorNumeric(
        palette = colorRampPalette(c("#35b7c4", "#4455a0"))(length(tab_corn$bio3)), 
        domain = tab_corn$bio3)
      leafletProxy("mymap", data = tab_corn) %>%
        clearShapes() %>%
        addCircleMarkers(lng = ~Longitud, lat = ~Latitud, color = ~pal(bio3),  fillColor = ~pal(bio3), 
                         fillOpacity = 1, popup = ~paste(bio3), radius = 4)
    }
    else if(input$coloring == "bio12"){
      pal <- colorNumeric(
        palette = colorRampPalette(c("#35b7c4", "#4455a0"))(length(tab_corn$bio12)), 
        domain = tab_corn$bio12)
      leafletProxy("mymap", data = tab_corn) %>%
        clearShapes() %>%
        addCircleMarkers(lng = ~Longitud, lat = ~Latitud, color = ~pal(bio12),  fillColor = ~pal(bio12), 
                         fillOpacity = 1, popup = ~paste(bio12), radius = 4)
    }
    else if(input$coloring == "Altitud"){
      pal <- colorNumeric(
        palette = colorRampPalette(c("#35b7c4", "#4455a0"))(length(tab_corn$Altitud)), 
        domain = tab_corn$Altitud)
      leafletProxy("mymap", data = tab_corn) %>%
        clearShapes() %>%
        addCircleMarkers(lng = ~Longitud, lat = ~Latitud, color = ~pal(Altitud),  fillColor = ~pal(Altitud), 
                         fillOpacity = 1, popup = ~paste(Altitud), radius = 4)
    } # if there were more variables we would repeat until they were completed
  })
  
  
  ###### TABLE ######
  
  datasetInput <- reactive({
    # Fetch the appropriate data object, depending on the value of input$dataset.
    switch(input$dataset,
           "Eigenvalues" = select(tab_corn, Indiv, EV1, EV2, EV3, EV4, EV5, EV6), #Selected eigenvalues of PCA
           "Coordinates" = select(tab_corn, Indiv, Longitud, Latitud), #Selected Coordinates 
           "Worlclime" = select(tab_corn, Indiv, bio1, bio2, bio3, bio4, bio5, bio6, bio7, bio8, bio9, bio10, 
                                bio11, bio12, bio13, bio14, bio15, bio16, bio17, bio18, bio19) #Selected worlclime variables
           
    )
  })
  
  output$table <- renderTable({
    datasetInput()
  })
  
  # downloadHandler() takes two arguments, both functions.
  # The content function is passed a filename as an argument, and
  # it should write out data to that filename.
  output$downloadData <- downloadHandler(
    
    # This function returns a string which tells the client
    # browser what name to use when saving the file.
    filename = function() {
      paste(input$dataset, input$filetype, sep = ".")
    },
    
    # This function should write data to a file given to it by
    # the argument 'file'.
    content = function(file) {
      sep <- switch(input$filetype, "csv" = ",", "tsv" = "\t")
      
      # Write to a file specified by the 'file' argument
      write.table(datasetInput(), file, sep = sep, 
                  row.names = FALSE)
    }
  )
  
  ###### Admixture ######
  
  data_admix <- reactive({get(input$Kvalue)}) # Select the dataframe from each K
  
  output$structure <- renderPlot({
    data_admix1 <- data_admix() 
    # Order levels of the column that ggplot uses in x so
    # that they are in the desired order
    data_admix1$INDIV <- reorder(data_admix1$INDIV, -data_admix1 [,input$coloring]) 
    # build graph with ggplot syntax
    ggplot(data_admix1 , aes_string(x="INDIV", y="Qadmixture", fill= "Kgroup")) + geom_col() +
      theme(axis.text.x= element_blank())
    
  })
  
}

## So far, we have defined both essential parts of a Shiny app - UI and Server.
## Finally, we have to call / run the Shiny, with UI and Server as its parameters.

shinyApp(
  ui = dashboardPage(header, sidebar, body), 
  server
)
