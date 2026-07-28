## PACIFIC CORAL REEF EXTENT

#load libraries
library(tidyverse)
library(sf)
library(h3)
library(h3js)

#define ancilliary functions
get_reef_area <- function(id, reef) {
  hex <- id %>% h3_to_geo_boundary_sf()
  r.area <- st_intersection(reef, hex$geometry) |> 
    st_area()
  if (length(r.area) == 0) {
    return(units::set_units(0, "m^2"))
  }
  return(sum(r.area))
}

#Load data
# reefs<-st_read("C:/Users/mgonzale/OneDrive - Australian Institute of Marine Science/GIS_Datasets/Reefs/MCRMP_reef500/WCMC008_CoralReef2018_Py_v4_1_valid.geojson") 
# ree
# subregions<-st_read("../../../gcrmn_regions/data/gcrmn-regions/gcrmn_subregions.shp") |> 
# rename(subregion=subregn)  |> 
# filter(region=="Pacific") 
# meow<-st_read("../../../gcrmn_regions/data/meow/Marine_Ecoregions_Of_the_World__MEOW_.shp")  |> 
# st_transform(st_crs(subregions)) |>
# st_intersection(subregions)
# st_write(meow, "data/pacific_meow.geojson", delete_dsn = TRUE )
reefs<-st_read("data/pacific_reefs.geojson")
meow<-st_read("data/pacific_meow.geojson")

# reefs<-reefs |> st_make_valid() |> 
# st_crop(subregions |> st_shift_longitude())

# st_write(reefs, "data/pacific_reefs.geojson", delete_dsn = TRUE)

# subregions.H3<-subregions |> group_by(subregion, geometry) |> 
# reframe(H3=polyfill(geometry, res=5)) |> ungroup() |> 
# group_by(subregion, H3) |> 
# st_drop_geometry()
# summarise(reef_area=get_reef_area(H3, reefs)) 

test<- st_drop_geometry(subregions.H3 |> ungroup()) |> as_tibble() |> 
group_by(subregion, H3) |>
summarise(H3=unique(H3))

meow.h3<-meow  |>  group_by(ECOREGION, geometry) |> select(ECOREGION, geometry) |>
st_shift_longitude() |> 
reframe(H3=polyfill(geometry, res=5)) |> ungroup() |> 
group_by(ECOREGION, H3) |> 
st_drop_geometry() |> 
mutate(reef_area=get_reef_area(H3, reefs)) 

# test<-meow |>  filter(ECOREGION=="Bismarck Sea")
# plot(test[,1])


meow_reef<-meow |> group_by(ECOREGION, geometry) |>
st_shift_longitude() |>
st_intersection(reefs) |> 
mutate(reef_area=st_area(geometry)) 



meow_reef<-meow_reef |> select(ECO_CODE, ECOREGION, reef_area) |> 
st_drop_geometry() |> 
group_by(ECO_CODE, ECOREGION) |> 
summarise(reef_area=sum(reef_area)) |>
ungroup() |>  
filter(!(ECO_CODE %in% c(20142,20051,20127,20129,20130,20139,20142))) |>  #Remove some ecoregions from East Asia (Eastern Phillipines, Kuroshi,Halmahera, Papua, Arafura Sea, Torres Strait)
 mutate(p.reefarea=as.numeric(reef_area)/sum(as.numeric(reef_area))) |> 
 arrange(-reef_area)

saveRDS(meow_reef, "data/meow_reef_extent.RDS")


## Convert reef polygons to H3 cells and extract area
h3_res <- 6

to_h3_cells <- function(poly_geom, res) {
  poly_sf <- st_sf(geometry = st_sfc(poly_geom, crs = 4326))

  if ("polyfill" %in% getNamespaceExports("h3jsr")) {
    return(h3jsr::polyfill(poly_sf, res = res, simple = TRUE))
  }

  if ("polygon_to_cells" %in% getNamespaceExports("h3jsr")) {
    return(h3jsr::polygon_to_cells(poly_sf, resolution = res, simple = TRUE))
  }

  stop("No polygon-to-H3 function found in h3jsr package.")
}

reef_h3_cells <- subregions |>
  st_transform(4326) |>
  st_make_valid() |>
  select(subregion) |>
  mutate(h3 = purrr::map(geometry, ~to_h3_cells(.x, h3_res))) |>
  st_drop_geometry() |>
  tidyr::unnest(h3) |>
  distinct(subregion, h3) |>
  mutate(h3_area_km2 = h3jsr::cell_area(h3, unit = "km2"))

reef_h3_area <- reef_h3_cells |>
  group_by(subregion) |>
  summarise(n_h3_cells = n(), reef_h3_area_km2 = sum(h3_area_km2), .groups = "drop")

write_csv(reef_h3_cells, "data/reefs_h3_cells.csv")
write_csv(reef_h3_area, "data/reefs_h3_area_by_subregion.csv")
