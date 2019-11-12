## cmd_pred
## Jamie Spaulding

#' Calculation of Center of Minimum Distance for Geographic Profiling
#' @description A calculation of the center of minimum distance for serial crime
#'     analysis. This function is among the centrographic methods which have
#'     been used for geographic profiling. The model assumes that the serial
#'     perpetrator's home base is relatively central among the crime incidents.
#' @param lat a vector of latitudes for the crime incident series
#' @param lon a vector of latitudes for the crime incident series
#' @return A latitude and longitude point of the center of minimum distance of
#'     the incidents. This mean can be used to prioritize the area which contains
#'     the offender's anchor point.
#' @author Jamie Spaulding, Keith Morris
#' @keywords spatial methods
#' @examples
#' #Using provided dataset for the Boston Strangler Incidents:
#' data(desalvo)
#' cmd_pred(desalvo$lat, desalvo$lon)
#' @importFrom aspace distances
#' @importFrom grDevices chull
#' @importFrom splancs gridpts
#' @export
cmd_pred <- function(lat, lon) {
  # Create and Initialize Objects/Counters
  i <- 1
  x <- c() #X-coord
  x[i] <- 0
  y <- c() #Y-coord
  y[i] <- 0
  d <- c() #Distance
  dist <- 0.00111 #Meters (111km in 1 lat)
  d[i] <- dist
  n <- c() #Iteration
  n[i] <- 0
  cells <- c() #Cells
  cells[i] <- 0
  dx <- 10 #Grid Spacing
  dy <- 10 #Grid Spacing

  # Min Convex Polygon -----
  points <- data.frame(lat, lon)
  hpts <- grDevices::chull(points)
  MCP <- cbind(points[hpts,1], points[hpts,2])

  while (d[i] >= dist) {
    grid <- splancs::gridpts(MCP, dx, dy)
    M.CMD <- matrix(0, nrow = nrow(grid), ncol = 3)
    for(j in 1:nrow(grid)) {
      coord.CMD <- grid[j, ]
      sumdist.CMD <- sum(aspace::distances(centre.xy = coord.CMD, points))
      M.CMD[j,1] <- sumdist.CMD
      M.CMD[j,2] <- coord.CMD[1]
      M.CMD[j,3] <- coord.CMD[2]
    }

    if (i >= 1) {
      order.CMD <- M.CMD[order(M.CMD[ ,1]), ]
      CMD <- order.CMD[1, ]
    }else (CMD <- M.CMD[1,])

    # CMD for each Iteration -----
    x[i + 1] <- CMD[2]
    y[i + 1] <- CMD[3]

    #Estimate Dist Between Current and Previous CMD -----
    d[i + 1] <- sqrt((x[i + 1] - x[i]) ^ 2 + (y[i + 1] - y[i]) ^ 2)
    n[i + 1] <- dx
    cells[i + 1] <- nrow(grid)
    rm(grid)
    i <- i + 1
    dx <- dx + 1
    dy <- dy + 1
  }

  # Result: CMD Holds the Min Dist Coord -----
  result <- cbind(n, round(x,2), round(y,2), round(d,2), cells)
  result <- as.data.frame(result[3:nrow(result), ])
  result[,1] <- seq(1, nrow(result), 1)
  colnames(result) <- c("n", "X", "Y", "Dist", "Cells")

  # CMD returned from simulation -----
  CMD <- cbind(result[nrow(result), 2:5])
  return(data.frame(lat = CMD[,1], lon = CMD[,2]))
}

