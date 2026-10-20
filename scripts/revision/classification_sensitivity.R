source(file.path("scripts","revision","common.R"))
library(sf)
directory <- file.path(output_dir,"classification")
dir.create(directory,recursive=TRUE,showWarnings=FALSE)
occ <- read.csv(file.path(root,"recovery_audit","published_2022_occurrences.csv"))
stopifnot(nrow(occ)==37,all(is.finite(occ$longitude)),all(is.finite(occ$latitude)))
occ_points <- st_as_sf(occ,coords=c("longitude","latitude"),crs=4326)
pop_points <- st_as_sf(pop_meta,coords=c("longitude","latitude"),crs=4326,remove=FALSE)
projected_occ <- st_transform(occ_points,32648)
projected_pop <- st_transform(pop_points,32648)
# This is an explicitly occurrence-derived boundary sensitivity, not a
# certified range boundary. It uses coordinates without climate/offset data.
hull <- st_convex_hull(st_union(projected_occ))
inside <- lengths(st_intersects(projected_pop,hull))>0
distance <- as.numeric(st_distance(projected_pop,st_boundary(hull)))/1000
signed_interior_distance <- ifelse(inside,distance,-distance)
membership <- list()
for(region in c("core","peri")) {
  polygon <- st_read(file.path(bundle,"ammo","gf",paste0("ammo_",region,".shp")),quiet=TRUE)
  if(is.na(st_crs(polygon))) st_crs(polygon) <- 4326
  valid_polygon <- st_make_valid(st_transform(polygon,32648))
  membership[[region]] <- lengths(st_intersects(projected_pop,valid_polygon))>0
}
polygon_class <- ifelse(membership$core & membership$peri,"overlap",ifelse(membership$core,"core",ifelse(membership$peri,"peripheral","other")))
definition <- data.frame(pop_meta,habitat_polygon_class=polygon_class,inside_occurrence_hull=inside,
  signed_distance_from_occurrence_hull_boundary_km=signed_interior_distance,
  geographic_marginality_km=-signed_interior_distance,
  geographic_boundary_definition="Convex hull of 37 published 2022 occurrence coordinates, EPSG:32648",
  limitation="Occurrence envelope sensitivity, not complete species range; uneven sampling and shared occurrence provenance remain")
write_csv(definition,file.path(directory,"population_class_and_geographic_marginality.csv"))
st_write(st_sf(geometry=hull),file.path(directory,"occurrence_hull_sensitivity.gpkg"),delete_dsn=TRUE,quiet=TRUE)
rona <- read.csv(file.path(output_dir,"rona","rona_population_long.csv"))
gdm <- read.csv(file.path(output_dir,"gdm","population_finite_domain_offsets.csv"))
gdm <- gdm[gdm$model=="geo_TRUE" & gdm$search_scoring=="environment_only",]
rows <- list();ri <- 0L
calculate <- function(d,y,analysis,variable,ssp,period,gcm) {
 x <- definition$geographic_marginality_km[match(d$population,definition$population)]
 stopifnot(length(x)==19,all(is.finite(x)),all(is.finite(y)))
 test <- suppressWarnings(cor.test(x,y,method="spearman",exact=FALSE))
 if(gcm=="ENSEMBLE") {
  boot <- replicate(2000,{j <- sample(seq_along(y),19,TRUE);suppressWarnings(cor(x[j],y[j],method="spearman"))})
  interval <- quantile(boot,c(.025,.975),na.rm=TRUE)
 } else interval <- c(NA_real_,NA_real_)
 ri <<- ri+1L
 rows[[ri]] <<- data.frame(analysis=analysis,variable=variable,ssp=ssp,period=period,gcm=gcm,n=19,rho=unname(test$estimate),CI_low=interval[1],CI_high=interval[2],p=test$p.value,
  marginality_definition="Negative signed interior distance to occurrence-only convex hull in km",interpretation="Exploratory coordinate-based sensitivity; no claim of certified range coverage or removal of spatial dependence")
}
for(ssp in unique(rona$ssp)) for(period in unique(rona$period)) for(gcm in unique(rona$gcm)) {
 for(bio in bio_names) {
  d <- rona[rona$ssp==ssp & rona$period==period & rona$gcm==gcm & rona$bio==bio,]
  calculate(d,d$rona,"RONA",bio,ssp,period,gcm)
 }
 d <- gdm[gdm$ssp==ssp & gdm$period==period & gdm$gcm==gcm,]
 for(metric in c("local","forward","reverse")) calculate(d,d[[metric]],"GDM",metric,ssp,period,gcm)
}
result <- do.call(rbind,rows);result$q_BH <- NA_real_
for(analysis in unique(result$analysis)) for(ssp in unique(result$ssp)) for(period in unique(result$period)) for(gcm in unique(result$gcm)) {
 idx <- result$analysis==analysis & result$ssp==ssp & result$period==period & result$gcm==gcm
 result$q_BH[idx] <- p.adjust(result$p[idx],"BH")
}
write_csv(result,file.path(directory,"continuous_geographic_marginality_sensitivity.csv"))
save_session(directory)
marker(directory,list(occurrences=37,populations=19,association_rows=nrow(result),population_polygon_disagreements=sum(polygon_class!=pop_meta$group),
  primary_population_definition="Verified 2022 published classes retained",alternative="Continuous occurrence-derived geographic marginality; no outcome-selected binary split",
  habitat_geometry_handling="Duplicate vertices repaired on projected working copies with GEOS; original polygon files retained",
  limitations="Climate-independent construction does not make uneven occurrence sampling, geographic population structure or uncertain true range boundary disappear."))
