## cgt_profile
## Jamie Spaulding

#' Criminal Geographic Targeting Model for Geographic Profiling (Rossmo Formula)
#' @description XXXXX
#' @param lat XXXX
#' @param lon XXXX
#' @param b XXXX
#' @param f XXXX
#' @param g XXXX
#' @return XXXX
#' @author Jamie Spaulding, Keith Morris
#' @keywords spatial methods
#' @examples
#' \donttest{
#' #Using provided dataset from Chicago Data Portal:
#' crimes <- data.frame(rcrimeanalysis:::crimes)
#' nr_data <- head(crimes,n=5000) #truncate dataset for near repeat analysis
#' library("rgdal")
#' out <- near_repeat_analysis(data=nr_data,tz="America/Chicago",epsg="32616")
#' path <- paste0(getwd(),"/netout") #path for iGraph networks out
#' name <- 1
#' # Save Image of Each igraph Network to Netpath Directory
#' library(igraph)
#' for(i in out){
#'     png(file = paste(path, "/series", name, ".png",sep = ""))
#'     plot(i, layout=layout_with_lgl, edge.color="orange",
#'     vertex.color="orange", vertex.frame.color="#ffffff",
#'     vertex.label.color="black")
#'     dev.off()
#'     name <- name+1
#' }}
#'
#' @importFrom sp SpatialPoints
#' @importFrom RANN.L1 nn2
#' @importFrom sp spTransform
#' @importFrom igraph graph_from_adjacency_matrix
#' @importFrom igraph components
#' @importFrom stats complete.cases
#' @importFrom stats dist
#' @export
cgt_profile <- function(lat, lon, buffer = NULL, f = NULL, g = NULL){
  # Set Defaults -----
  if (is.null(buffer)) {
    # Calculate Incident Buffer Zone -----
    dat_nn <- cbind(lat,lon) # Extract only lat and lon columns
    nn_list <- RANN.L1::nn2(dat_nn, dat_nn, k=2) # Find nearest neighbors using Manhattan distance
    nn <- nn_list$nn.idx # Extract NN pairs
    
    # Calculate Manhattan Distances Between NN Pairs -----
    nn_md <- NULL
    jj <- 1
    for(i in 1:nrow(nn)){
      incid1 <- dat_nn[nn[i,1],]
      incid2 <- dat_nn[nn[i,2],]
      dx <- geosphere::distHaversine(p1 = c(incid1[2], incid1[1]), p2 = c(incid2[2], incid1[1]), r = 3958) # hold y (lat) constant
      dy <- geosphere::distHaversine(p1 = c(incid1[2], incid1[1]), p2 = c(incid1[2], incid2[1]), r = 3958) # hold x (lon) constant
      nn_md[jj] <- dx + dy
      jj <- jj+1
    }
    buffer <- (mean(nn_md)) / 2} #default: 1/2 Mean NN Manhattan Distance
  c <- length(lat)
  if (is.null(f)) {f <- 1.9} #default: Rossmo (1995)
  if (is.null(g)) {g <- 1.9} #default: Rossmo (1995)
  
  # Computation of Map Boundaries/ Hunting Area -----
  # Rossmo, D. Kim. Geographic profiling. CRC press, 1999. pg 199
  lat_max <- max(lat) + ((max(lat) - min(lat)) / (2 * (length(lat) - 1)))
  lat_min <- min(lat) - ((max(lat) - min(lat)) / (2 * (length(lat) - 1)))
  lon_max <- max(lon) + ((max(lon) - min(lon)) / (2 * (length(lon) - 1)))
  lon_min <- min(lon) - ((max(lon) - min(lon)) / (2 * (length(lon) - 1)))
  
  # Calculate Range of Bounding Box -----
  lat_range <- lat_max - lat_min
  lon_range <- lon_max - lon_min
  
  # Determine Sequence of Lat and Lon Gridlines -----
  lats <- seq(lat_min,lat_max, length.out = 200)
  lons <- seq(lon_min,lon_max, length.out = 200)
  
  # Create a Run Sequence for Each Incident of Grid Points -----
  run_seq <- expand.grid(lats,lons)
  names(run_seq) <- c("lats", "lons")
  
  # CGT Distance Decay Function -----
  jj <- 1
  phi <- NULL
  output <- data.frame()
  # PROGRESS BAR FOR LOOP
  pb = txtProgressBar(min = 0, max = length(lat)*40000, style = 3)
  tick <- 0
  
  for(i in 1:length(lat)){
    for(j in 1:nrow(run_seq)){
      tick <- tick + 1
      setTxtProgressBar(pb, tick)
      xn <- lon[i]
      yn <- lat[i]
      xi <- run_seq$lons[j]
      yi <- run_seq$lats[j]
      dx <- geosphere::distHaversine(p1 = c(xn, yi), p2 = c(xi, yi), r = 3958) #hold y (lat) constant
      dy <- geosphere::distHaversine(p1 = c(xi, yn), p2 = c(xi, yi), r = 3958) #hold y (lat) constant
      if(dx + dy > 1){phi <- 1} else {phi <- 0}
      output[jj,i] <- (phi / ((dx + dy) ^ f)) + (((1 - phi) * (buffer ^ (g - f))) / (((2 * buffer) - (dx - dy)) ^ g)) #Rossmo Formula
      jj <- jj + 1
    }
    jj <- 1
  }
  
  ## Summation of Values for Each Grid Point
  sums <- rowSums(output, na.rm = TRUE)
  dat <- cbind(sums, run_seq)
  return(dat)
}


data <- read.csv("G:/Workspace/PhD_Research/Prediction/new/prediction_alt_methods/cases/desalvo_incidents.csv")
data <- head(data,n=13)
lat <- data$lat
lon <- data$lon
buffer=NULL
f=NULL
g=NULL


g_map = sp::SpatialPixelsDataFrame(points = dat[c("lons", "lats")], data = dat)
g_map <- raster::raster(g_map)
raster::crs(g_map) <- sp::CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs") # assign a Coordinate Reference System

# Parula color pallete for jeopardy surface
library(leaflet)
pal <- colorNumeric(pals::parula(200), raster::values(g_map),
                    na.color = "transparent")

leaflet() %>%
  addTiles() %>%
  addProviderTiles("Esri.WorldTopoMap", group = "Topo") %>%
  addAwesomeMarkers(lng = -71.07357, lat = 42.41322, icon = awesomeIcons(icon = 'home', markerColor = 'green'), popup =
                      "Residence") %>%
  addRasterImage(g_map, colors = pal, opacity = 0.6) %>%
  addLegend(pal = pal, values = raster::values(g_map), title = "Score") %>% #adds legend when using likelihood (k = 1)
  addCircleMarkers(lng = data$lon, lat = data$lat, radius = 4, opacity=1, fill = "black", stroke=TRUE,
                   fillOpacity = 0.75, weight=2, fillColor = "red")
