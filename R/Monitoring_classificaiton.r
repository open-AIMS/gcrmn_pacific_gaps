### -------------------------------------------------------------------------
### Script: Monitoring_classificaiton.r
### Author: Manuel González-Rivero
### Organisation: Australian Institute of Marine Science
### Date: 08/06/2026
### Purpose: Classify Pacific ecoregions by monitoring coverage and reef impact
### Approach:
###   1) Aggregate nearby monitoring records into site clusters (100 m buffers)
###   2) Derive coverage metrics per ecoregion (site count and temporal depth)
###   3) Combine coverage with relative reef extent to assign priority classes
###   4) Map resulting priorities and exploratory clusters across the Pacific
### Inputs: data/gcrmn_data-pacific.csv, data/pacific_meow.geojson,
###         data/meow_reef_extent.rds
### Outputs: fig/monitoring_risk_map_meow.png
### -------------------------------------------------------------------------

# Libraries
library(tidyverse)
library(sf)
library(h3jsr)
library(mregions2)
library(rnaturalearth)
library(rnaturalearthdata)

# Load datasets
m.df<-read_csv("data/gcrmn_data-pacific.csv") ##Monitoring metadata to estimate sampling efforts
# re<-read_csv("data/reefs_extent.csv")  |> ##Reef extent data to extract reef extent for each subregion
#   filter(region=="Pacific", subregion!="All")

# subregions<-st_read("../../../gcrmn_regions/data/gcrmn-regions/gcrmn_subregions.shp") |> 
# rename(subregion=subregn)  |> 
# filter(region=="Pacific") 

meow<-st_read("data/pacific_meow.geojson")
re.meow<-readRDS("data/meow_reef_extent.rds") |> rename(ecoregion=ECOREGION)


## SECTION 1: Build site clusters from point records
## Aggregate sites into 100 m buffers to account for spatial uncertainty and
## to group nearby records as a single monitoring location.

# Create buffers around sites
m.df<- m.df  |> st_as_sf(coords = c("decimalLongitude", "decimalLatitude"), crs = 4326) |>
 st_transform(crs = 3857) |>
 st_buffer(dist = 100)

## Cluster sites based on overlapping buffers
cluster<-st_union(m.df) |> # |> st_make_grid(cellsize = c(100, 100)) |> st_sf() |> st_transform(crs = 3857)
st_cast("POLYGON") |> st_sf() |> st_transform(crs = 4326) |> 
mutate(site_id=1:n())  |> 
st_make_valid() 

## Join clusters back to sites and summarize temporal coverage per site cluster
ecoregion_sites <- m.df  |> 
st_transform(crs = 4326)  |> 
st_join( cluster, join = st_intersects) |> 
st_drop_geometry() |>
select(ecoregion, year, site_id) |>
group_by(ecoregion, site_id) |>
summarise(n_years=length(unique(year)), length_surveys=max(year)-min(year)+1) |> 
ungroup() 


## SECTION 2: Summarize monitoring coverage and assign monitoring priority
## Metrics:
##   n.sites = number of unique site clusters in each ecoregion
##   lc      = longest yearly coverage observed in any clustered site
##   tc      = number of clusters with >10 years of surveys
## Priorities combine probability of detecting change (coverage) and impact
## (relative reef extent) into classes A-D.
df<-ecoregion_sites |> 
group_by(ecoregion) |> summarise(n.sites=n(), lc=max(n_years), tc=length(length_surveys[length_surveys>10])) |> 
full_join(re.meow, by="ecoregion") |>
#filter(!(subregion %in% c("All", "PEAll")), region!="Global") |>
mutate(n.sites=replace_na(n.sites, 0), lc=replace_na(lc, 0), tc=replace_na(tc, 0)) |>
mutate(Prob.detect=case_when(
 n.sites>200 & tc >30 ~ "High",
 n.sites<100 ~ "Low",
 tc<10 ~ "Low",
 TRUE ~ "Medium"), 
Impact=case_when(
  p.reefarea > quantile(p.reefarea, 0.75) ~ "High",
  p.reefarea <= quantile(p.reefarea, 0.25) ~ "Low",
 TRUE ~ "Medium"),
 Priority=case_when(
  Prob.detect=="High" ~ "D",
  Prob.detect=="Low" & Impact=="High" ~ "A",
  Prob.detect=="Low" & Impact=="Medium" ~ "B",
  Prob.detect=="Low" & Impact=="Low" ~ "C",
Prob.detect=="Medium" & Impact=="High" ~ "B",
  Prob.detect=="Medium" & Impact=="Medium" ~ "C",
  Prob.detect=="Medium" & Impact=="Low" ~ "D",
  TRUE ~ "D")) |> 
  arrange(Priority) |> 
  rename(p.extent=p.reefarea) |>
  select(ecoregion, n.sites, lc, tc, p.extent, Prob.detect, Impact, Priority)


## SECTION 3: Map ecoregion priorities
worldMap <- ne_countries(scale = "medium", returnclass = "sf")
target_crs <- st_crs("+proj=eqc +x_0=0 +y_0=0 +lat_0=0 +lon_0=155")

# define a long & slim polygon that overlaps the meridian line & set its CRS to match
# that of world
# Centered in lon 133

offset <- 180 - 155
polygon <- st_polygon(x = list(rbind(
  c(-0.0001 - offset, 90),
  c(0 - offset, 90),
  c(0 - offset, -90),
  c(-0.0001 - offset, -90),
  c(-0.0001 - offset, 90)
))) %>%
  st_sfc() %>%
  st_set_crs(4326)

world2 <- worldMap %>% st_difference(polygon)
#> Warning: attribute variables are assumed to be spatially constant throughout all
#> geometries

# Transform
world3 <- world2 %>% st_shift_longitude() |> st_transform(crs = target_crs)

# Bounding box of Pacific ecoregions in the display CRS for map cropping
pacific_bbox <- meow |>
  st_shift_longitude() |>
  st_transform(crs = "+proj=robin +lon_0=155") |>
  st_buffer(1000000) |>  # Add a buffer to ensure we capture all relevant areas
  st_bbox()

risk_map <- meow  |> 
rename(ecoregion=ECOREGION) |> 
left_join(df, by="ecoregion")|> 
mutate(Priority=factor(Priority, levels=c("A", "B", "C", "D"))) |>
st_as_sf() |>
 st_shift_longitude() |>
 st_transform(crs=target_crs) |>
    ggplot()+  
    geom_sf(aes(fill=Priority), color="black")+
    geom_sf(data = world3, fill = "lightgray", color = "darkgray") +
    coord_sf(
      crs = "+proj=robin +lon_0=155",
      datum = st_crs(4326),
      xlim = c(pacific_bbox["xmin"], pacific_bbox["xmax"]),
      ylim = c(pacific_bbox["ymin"], pacific_bbox["ymax"]),
      expand = FALSE,
      label_axes = list(bottom = "E", left = "N")
    ) +
    scale_fill_manual(name="Classificaiton", values=c( "#ff1010", "#ffc400" , "#a5eb03", "#1acc02"), na.value = "lightgrey")+
    theme_minimal()+
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
      axis.ticks = element_line(color = "black"),
      axis.ticks.length = grid::unit(-0.15, "cm")
    ) +
    labs(#title="Priorities for additional monitoring efforts across ecoregions",
          #subtitle="Based on monitoring coverage (time and space) vs reef extent",
         fill="Priority")

ggsave("fig/monitoring_risk_map_meow.png", risk_map, width=10, height=6)

### Cluster ecoregions based on monitoring coverage and reef extent to identify priority areas for additional monitoring efforts
## SECTION 4: Exploratory clustering of ecoregions
## This unsupervised step groups ecoregions with similar monitoring and reef
## extent profiles to support interpretation of priority patterns.
hclust<-df |> 
select(n.sites, lc, tc, p.extent) |>  
scale() |>
dist() |>
hclust(method="ward.D2")
#plot(hclust)

cuts<-cutree(hclust, k=4)
df$cluster<-cuts

clusr_map <- meow  |> 
rename(ecoregion=ECOREGION) |> 
left_join(df, by="ecoregion")|> 
mutate(cluster=factor(cluster, levels=c(1, 2, 3, 4))) |>
st_as_sf() |>
 st_shift_longitude() |>
 st_transform(crs=target_crs) |>
    ggplot()+  
    geom_sf(aes(fill=cluster), color="black")+
    geom_sf(data = world3, fill = "lightgray", color = "darkgray") +
    coord_sf(
      crs = "+proj=robin +lon_0=155",
      datum = st_crs(4326),
      xlim = c(pacific_bbox["xmin"], pacific_bbox["xmax"]),
      ylim = c(pacific_bbox["ymin"], pacific_bbox["ymax"]),
      expand = FALSE,
      label_axes = list(bottom = "E", left = "N")
    ) +
    scale_fill_manual(name="Classificaiton", values=c( "#ff1010", "#ffc400" , "#a5eb03", "#1acc02"), na.value = "lightgrey")+
    theme_minimal()+
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
      axis.ticks = element_line(color = "black"),
      axis.ticks.length = grid::unit(-0.15, "cm")
    ) +
    labs(title="Monitoring efforts across ecoregions",
          subtitle="Hierarchical clustering based on number of sites, temporal coverage, survey frequency, and reef extent",
         fill="Cluster")

##Summarise cluster profiles
cluster_summary<-df |>  
group_by(cluster) |>
summarise(n.ecoregions=n(), n.sites=median(n.sites), lc=median(lc), tc=median(tc), p.extent=median(p.extent)) |>
arrange(cluster)  
