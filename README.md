# rgeoprofile: Geographic Profiling Methods for Serial Crime Analysis

## Overview

**rgeoprofile** is a package containing various functions for the analysis of serial crime incidents. The package implements algorithms for the geographical profiling of serial incidents in attempt to prioritize the area in which the anchor point or home base of the perpetrator is located. These functions include centrographic methods (center of minimum distance, geometric and harmonic means) and distance decay functions in the following forms: Rossmo's formula, linear, negative exponential (CrimeStat and Dragnet, with/out plateau and buffer), normal, lognormal, and truncated negative exponential.
	
**Keywords**: Crime Analysis, Geographic Profiling, Distance Decay, Open Source

## Getting Started

These instructions will allow any user to install the R package on your local machine for development and testing purposes. 

### Installing

Begin by installing the package from the GitHub repository:

```
devtools::install_github('JSSpaulding/rgeoprofile')
```

After installation, load and attach the package

```
library(rgeoprofile)
```

The following will call example serial crime incident data from within the package for the Boston Strangler (Albert DeSalvo) case. This data structure includes the following columns: victim name, age, date of incident, location (street address), latitude, and longitude.

```
desalvo <- data.frame(rgeoprofile:::boston_strangler)
```

## Distance Decay Model Usage

The following examples illustrate the usage of the various distance decay functions within the **rgeoprofile** package.

### cgt_profile

The *cgt_profile* function applies the criminal geographic targeting (CGT) model for serial crime analysis developed by DK Rossmo for geographic profiling and prediction of perpetrator anchor point/home base. The CGT distance decay function assumes a buffer zone around the incidents where the likelihood of perpetrator anchor point increases to a peak likelihood before decaying as distance further increases. This model also utilizes Manhattan distance for the likelihood calculation. The formula for the distance decay function is:

\begin{equation}
p_{i,j} = \sum_{n=1}^{n} \frac{\phi_{i,j}}{(\mid X_{i} - x_{n} \mid + \mid Y_{i} - y_{n} \mid) ^ f} + \frac{(1 - \phi_{i,j}) (B^{g-f})}{(2B - \mid X_{i} - x_{n} \mid - \mid Y_{i} - y_{n} \mid) ^ g}
\end{equation}

When plotted, the likelihood of the perpetrator's anchor point is illustrated by the following:

![](../master/images/Rossmo_ddf.png)

```
desalvo <- data.frame(rgeoprofile:::boston_strangler)
test <- cgt_profile(desalvo$lat, desalvo$lon)
```

#### Mapping the Resultant Geographic Profile
Using the functions in the *rgeoprofile* package, the resultant geographic profile or jeopardy surface (three-dimensional probability surfaces which indicates the most probable area of offender anchor point) are returned for the input case and can then be mapped. The following will allow for the development of a *leaflet* map which plots the incidents and output from above:
```
g_map = sp::SpatialPixelsDataFrame(points = test[c("lons", "lats")], data = test)
g_map <- raster::raster(g_map)
# Assign a Coordinate Reference System for the Raster
raster::crs(g_map) <- sp::CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
# Define a Parula Color Pallete for Resultant Jeopardy Surface
library(leaflet) #for mapping the geographic profile
pal <- colorNumeric(pals::parula(200), raster::values(g_map),
    na.color = "transparent")
leaflet() %>%
    addTiles() %>%
    addProviderTiles('Esri.WorldTopoMap', group = 'Topo') %>%
    addAwesomeMarkers(lng = -71.07357, lat = 42.41322, icon =
        awesomeIcons(icon = 'home', markerColor = 'green'), popup = 'Residence') %>%
    addRasterImage(g_map, colors = pal, opacity = 0.6) %>%
    addLegend(pal = pal, values = raster::values(g_map), title = 'Score') %>%
    addCircleMarkers(lng = data$lon, lat = data$lat, radius = 4, opacity = 1,
        fill = 'black', stroke = TRUE, fillOpacity = 0.75, weight = 2,
        fillColor = "red")
```

![](../master/images/rossmo_geoprofile.png)

### neg_exp_profile
The *neg_exp_profile* function applies variations of the negative exponential decay model for serial crime analysis. In this model, the decline is at a constant rate, therefore the likelihood of the perpetrator's home base drops quickly from the incident locations until it approaches zero likelihood. The user can select different variants including the CrimeStat base model, the Dragnet model, or whether a buffer and plateau is present at the start of the decay function. This model assumes that the likelihood of the serial perpetrator's home base decreases in a exponential fashion as the distance increases from the crime incidents. This model utilizes Euclidean distance for the likelihood calculation. The general formula for the distance decay function is:

\begin{equation}
p_{i,j} = ae^{-bd_{i,j}}
\end{equation}

An illustration of the difference between the base negative exponential function in comparison to the inclusion of a plateau and buffer is shown in the following:

![](../rgeoprofile/images/neg_exp_ddf.png)

```
desalvo <- data.frame(rgeoprofile:::boston_strangler)
test <- neg_exp_profile(desalvo$lat, desalvo$lon, method = "CrimeStat")
## See Mapping the Resultant Geographic Profile above for plotting the results
```

### Other Distance Decay Models
The following image provides a comparison of the different distance decay functions. Descriptions of the different functions will occur with usage examples thereafter. 

![](../master/images/crimestat_ddf.png)

#### linear_profile

The *linear_profile* function applies the linear decay model for serial crime analysis within CrimeStat. This model assumes that the likelihood of the serial perpetrator's home base decreases in a linear fashion as the distance increases from the crime incidents. This model utilizes Euclidean distance for the likelihood calculation. The formula for the distance decay function is:

\begin{equation}
p_{i,j} = a + bd_{i,j}
\end{equation}

```
desalvo <- data.frame(rgeoprofile:::boston_strangler)
test <- linear_profile(desalvo$lat, desalvo$lon)
## See Mapping the Resultant Geographic Profile above for plotting the results
```

#### lognorm_profile

The *lognorm_profile* function applies the lognormal decay model for serial crime analysis within CrimeStat. This model is very similar to the normal model except with more skew to either side. If there is reason to believe that the perpetrator's residence is closer to the incidents, this function can take the form of a very rapid increase near incident with a gradual decline from the peak likelihood. This model utilizes Euclidean distance for the likelihood calculation. The formula for the distance decay function is:

\begin{equation}
p_{i,j} = a \frac{1}{d_{i,j}^2 S \sqrt{2\pi}}^\frac{[\ln{(d_{i,j}^2 - \overline{d})]^2}}{2S^2}
\end{equation}

```
desalvo <- data.frame(rgeoprofile:::boston_strangler)
test <- lognorm_profile(desalvo$lat, desalvo$lon)
## See Mapping the Resultant Geographic Profile above for plotting the results
```

#### norm_profile

The *norm_profile* function applies the normal decay model for serial crime analysis within CrimeStat. This model assumes that there is a peak likelihood of the serial perpetrator's home base at some optimal distance from the crime incidents. The function rises in likelihood to that distance and then declines at an equal rate (both prior to and after the peak likelhihood) giving the symetrical normal distribution. This model utilizes Euclidean distance for the likelihood calculation. The formula for the distance decay function is:

\begin{equation}
p_{i,j} = a \frac{1}{S \sqrt{2\pi}}e^{\frac{-1}{2}(\frac{(d_{i,j} - \overline{d})}{S})^2}
\end{equation}

```
desalvo <- data.frame(rgeoprofile:::boston_strangler)
test <- norm_profile(desalvo$lat, desalvo$lon)
## See Mapping the Resultant Geographic Profile above for plotting the results
```

#### trun_neg_exp_profile

The *trun_neg_exp_profile* function applies the truncated negative exponential decay model for serial crime analysis within CrimeStat. This is a joint function composed of both the linear and the negative exponential. For distances proximal to the incidents, a positive linear function is defined from zero likelihood at distance zero to a location of peak likelihood. At the peak likelihood the function takes the form of a negative exponential, rapidly declining as distance increases. This model utilizes Euclidean distance for the likelihood calculation. The formula for the distance decay function is:

\begin{equation}
\begin{aligned}
p_{i,j} & = bd_{i,j}\ for\ 0 \leq d_{i,j} \leq d_{p} \\
& = ae^{-cd_{i,j}}\ for\ d_{i,j} > d_{p}
\end{aligned}
\end{equation}

```
desalvo <- data.frame(rgeoprofile:::boston_strangler)
test <- trun_neg_exp_profile(desalvo$lat, desalvo$lon)
## See Mapping the Resultant Geographic Profile above for plotting the results
```

## Centrogrpahic Model Usage

The following examples illustrate the usage of the various centrographic functions within the **rgeoprofile** package.

### cmd_pred

The *cmd_pred* function calculates the center of minimum distance (CMD) for a set of incident coordinates. This function returns a singular coordinate for the CMD which can be used for geographic profiling or area prioritization under the assumption that the home base of the perpetrator is centrally located among the incidents.

```
desalvo <- data.frame(rgeoprofile:::boston_strangler)
cmd_pred(desalvo$lat, desalvo$lon)
```

### geom_mean_pred

The *geom_mean_pred* function calculates the geometric mean for a set of incident coordinates. This function returns a singular coordinate which can be used for geographic profiling or area prioritization under the assumption that the home base of the perpetrator is centrally located among the incidents.

```
desalvo <- data.frame(rgeoprofile:::boston_strangler)
geom_mean_pred(desalvo$lat, desalvo$lon)
```

### harm_mean_pred

The *harm_mean_pred* function calculates the harmonic mean for a set of incident coordinates. This function returns a singular coordinate which can be used for geographic profiling or area prioritization under the assumption that the home base of the perpetrator is centrally located among the incidents.

```
desalvo <- data.frame(rgeoprofile:::boston_strangler)
harm_mean_pred(desalvo$lat, desalvo$lon)
```


## Authors

* **Jamie Spaulding** - *Author, Contributor*

* **Keith Morris** - *Contributor*
