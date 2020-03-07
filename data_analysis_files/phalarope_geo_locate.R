# Load librarieis
library(GeoLocTools)
#setupGeolocation()
library(dplyr)


# Migrate Technology Ltd logger
# Type: 4.12.4, current settings: Approx 0.5'C resolution. Conductivity >63 for 'wets' count. Light range 4, ambient. XT. 
# Logger number: BA055

# MODE: 3 (full range light, temperature, wet/dry and conductivity recorded)
# LIGHT: Sampled every minute with max light recorded every 5mins.
# TEMPERATURE: Sampled every 5mins with max and min recorded every 4hrs.
# WET/DRY and CONDUCTIVITY: Sampled every 30secs with number of samples wet  recorded every hr (capped at 7mins wet per hour), max conductivity recorded every 4hrs.
# Max record length = 15 months. Total battery life upto 23 months. Logger is currently 14 months old.

# Programmed: 07/07/2016 22:13:14. Start of logging (DD/MM/YYYY HH:MM:SS): 07/07/2016 22:13:14
# Age at start of logging (secs): 4488620, approx 2 months
# End of logging (DD/MM/YYYY HH:MM:SS): 01/07/2017 08:02:17
# Age at end of logging (secs): 35454822, approx 14 months
# Timer (DDDDD:HH:MM:SS): 00358:09:43:22
# Drift (secs): 341. Memory not full. 
# Pointers: 64826,123696,3158064,3158064,0,0,61% light mem used,70% wet/dry/temp mem used
# Tcals (Ax^3+Bx^2+Cx+D): -35407.051,64980.930,-40521.832,8621.896
# Approx 0.5'C resolution. Conductivity >63 for 'wets' count. Light range 4, ambient. XT. 
# Load raw data

# Read in raw data
raw <- read.csv("Phalarope.csv")

# Create datasets for geo-location and additional measurements
d.lux <- raw %>% 
  # Select omly data related too light levels
  dplyr::filter(sensor.type == "solar-geolocator-raw")%>%
  # Convert data too POSIX format
  mutate(Date = as.POSIXct(.$Date, tz = "GMT")) %>%
  # Select only only one track (there are two in the file)
  filter(tag.local.identifier == "BA055") %>%
  # Log transfrom light intensity
  mutate(Light = log(.$Light)) %>%
  # Selecy only those columns required
  dplyr::select(Date, Light)

head(d.lux)


ID <- "BA055"
Species <- "Phalarope"

# BA055's departure/capture point
lon.calib <- 177.05
lat.calib <- 62.55


  
# Twilight annotation
threshold <- 2

col = colorRampPalette(c('black',"purple",'orange'))(50)[as.numeric(cut(geo_locator[3000:10000,2],breaks = 50))]

par(mfrow = c(1, 1), 
mar = c(2, 2, 2, 2) )
with(geo_locator[3000:10000,], plot(Date, Light, type = "o", pch=16,  col = col, cex = 0.5)) 
abline(h=threshold, col="orange", lty = 2, lwd = 2)

# Plot light image
offset <- 4 # adjusts the y-axis to put night (dark shades) in the middle
lmax <- 12

# Light equals night time
# Dark equals day time
par(mfrow = c(1, 1), 
    mar = c(2, 2, 2, 2) )
lightImage(tagdata = geo_locator,
           offset = offset,     
           zlim = c(0,lmax),
           ylab = "Hour (GMT)")


# Calculate transitions
BA055_transitions <- twilightCalc(datetime = geo_locator$Date,
                                  light= geo_locator$Light,
                                  LightThreshold= threshold,    
                                  ask=FALSE) 

head(BA055_transitions)


# Calculate Twilight data 
twl <- preprocessLight(geo_locator, 
                       threshold = threshold,
                       offset = offset, 
                       lmax = lmax,       # max. light value
                       gr.Device = "x11") # MacOS version

# Account for storing maximum values for each 5 minute period. This step will depend on tag type. 
twl_adj <- twilightAdjust(twl, 150)

### Show results of defining twilights against raw light data
lightImage(geo_locator, offset = offset, zlim = c(0,lmax))
tsimagePoints(twl_adj$Twilight, offset = offset, pch = 16, cex = 0.5,
              col = ifelse(twl_adj$Rise, "dodgerblue", "firebrick"))









# Clean twilight times
twl <- twilightEdit(twilights = twl,
                    offset = offset,
                    window = 4,           # two days before and two days after
                    outlier.mins = 45,    # difference in mins
                    stationary.mins = 25, # are the other surrounding twilights within 25 mins of one another
                    plot = TRUE)

# Let's have a look at the combined image set
offset <- 4 # adjusts the y-axis to put night (dark shades) in the middle

lightImage( tagdata = geo_locator,
            offset = offset,     
            zlim = c(0, 4))

tsimagePoints(twl$Twilight, 
              offset = offset, 
              pch = 16, 
              cex = 1.2,
              col = ifelse(twl$Deleted, "grey20", ifelse(twl$Rise, "firebrick", "cornflowerblue")))

# Filter out those twilight events that were marked as deleted
twl <- subset(twl, !Deleted)

# Convert the format that it can be used with Geolight
twl.gl  <- export2GeoLight(twl) 
head(twl.gl)

# Calibration
lightImage(tagdata = geo_locator,
           offset = offset,     
           zlim = c(0, 12))

tsimageDeploymentLines(twl$Twilight, lon.calib, lat.calib, offset = offset,
                       lwd = 2, col = adjustcolor("orange", alpha.f = 0.8))

# Calibration times
tm1 <- c(as.POSIXct("2016-07-08"), as.POSIXct("2016-07-23"))
tm2 <- c(as.POSIXct("2017-06-21"), as.POSIXct("2017-07-01"))

abline(v = tm1, lty = c(1,2), col = "firebrick", lwd = 1.5)
abline(v = tm2, lty = c(1,2), col = "firebrick", lwd = 1.5)

# Subset the twilight data with calibration periods
d.calib <- subset(twl, (tFirst>=tm1[1] & tSecond<=tm1[2]) | (tFirst>=tm2[1] & tSecond<=tm2[2]))

# Calibrate the data
gE <- getElevation(twl = twl, 
                   known.coord = c(lon.calib, lat.calib), 
                   method = "gamma")

