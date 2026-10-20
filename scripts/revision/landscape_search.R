source(file.path("scripts","revision","common.R"))
library(terra)
library(gdm)
library(geosphere)
directory <- file.path(output_dir,"landscape")
dir.create(directory,recursive=TRUE,showWarnings=FALSE)
model <- readRDS(file.path(output_dir,"gdm","model_geo_TRUE.rds"))
base <- file.path(bundle,"ammo","gf","chelsa2_crop_selected")
read_stack <- function(folder,native=FALSE) {
  x <- rast(file.path(folder,if(native) paste0(bio_names,".tif") else paste0("b",bio_numbers,".asc")))
  names(x) <- bio_names;x
}
current <- read_stack(file.path(base,"current"))
range <- vect(file.path(bundle,"ammo","gf","ammo_binary.shp"))
if(!nzchar(crs(range))) {
  stopifnot(all(abs(ext(range)[1:2])<=180),all(abs(ext(range)[3:4])<=90))
  crs(range) <- "EPSG:4326"
  range_crs_note <- "EPSG:4326 assigned from documented legacy code and degree-coordinate bounds; original shapefile has no CRS"
} else range_crs_note <- "Original shapefile CRS retained"
masked <- mask(current,range)
values_current <- values(masked,mat=TRUE)
valid_cells <- which(apply(is.finite(values_current),1,all))
xy <- xyFromCell(masked,valid_cells)
current_grid <- data.frame(longitude=xy[,1],latitude=xy[,2],values_current[valid_cells,,drop=FALSE])
names(current_grid)[3:7] <- bio_names
pop_current <- read.csv(file.path(output_dir,"climate","current_population_climate.csv"))
pop_current <- pop_current[match(pop_meta$population,pop_current$population),]
# Known sampled sites augment the habitat candidate domain. Preserve exact
# coordinates and explicit membership rather than silently snapping or dropping.
candidate <- rbind(current_grid,pop_current[,c("longitude","latitude",bio_names)])
n_grid <- nrow(current_grid);n_candidate <- nrow(candidate)
write_csv(data.frame(candidate_index=seq_len(n_candidate),candidate,
    candidate_type=c(rep("habitat_grid_cell",n_grid),rep("verified_sampled_population",19)),
    population=c(rep("",n_grid),pop_meta$population)),file.path(directory,"current_candidate_domain.csv"))
distances <- lapply(seq_len(19),function(i) geosphere::distGeo(matrix(c(pop_meta$longitude[i],pop_meta$latitude[i]),n_candidate,2,byrow=TRUE),as.matrix(candidate[,c("longitude","latitude")]))/1000)
pair_frame <- function(first,second,penalty) {
  n <- nrow(second)
  cbind(data.frame(distance=rep(0,n),weights=rep(1,n),s1.xCoord=rep(first$longitude,n),s1.yCoord=rep(first$latitude,n),
      s2.xCoord=if(penalty) second$longitude else rep(first$longitude,n),s2.yCoord=if(penalty) second$latitude else rep(first$latitude,n)),
      setNames(first[rep(1,n),bio_names],paste0("s1.",bio_names)),setNames(second[,bio_names],paste0("s2.",bio_names)))
}
score <- function(first,second,penalty) {
  raw <- as.numeric(predict(model,pair_frame(first,second,penalty)))
  stopifnot(all(is.finite(raw)),all(raw<1))
  eta <- -log1p(-raw)-model$intercept
  stopifnot(all(eta> -1e-8))
  1-exp(-pmax(0,eta))
}
radii <- c(Inf,25,50,100)
gcms <- c("ENSEMBLE","GFDL-ESM4","IPSL-CM6A-LR","MPI-ESM1-2-HR","MRI-ESM2-0","UKESM1-0-LL")
all_results <- list(); ai <- 0L
for(ssp in c("ssp126","ssp585")) for(period in c("2041-2070","2071-2100")) for(gcm in gcms) {
  name <- paste(ssp,period,gcm,sep="_")
  checkpoint <- file.path(directory,paste0(name,"_complete.json"))
  csv_path <- file.path(directory,paste0(name,"_matches.csv"))
  if(file.exists(checkpoint) && file.exists(csv_path)) {
    ai <- ai+1L;all_results[[ai]] <- read.csv(csv_path)
    message("Reusing validated landscape checkpoint ",name);next
  }
  short <- if(period=="2041-2070") "2040" else "2070"
  future <- if(gcm=="ENSEMBLE") read_stack(file.path(base,paste0(ssp,"_",short))) else read_stack(file.path(root,"public_dependencies","climate_individual","CHELSA_v2.1_regional",period,gcm,ssp),TRUE)
  fvalues <- terra::extract(future,candidate[,c("longitude","latitude")],method="simple",ID=FALSE)
  stopifnot(all(is.finite(as.matrix(fvalues))))
  future_candidate <- data.frame(candidate[,c("longitude","latitude")],fvalues)
  local <- vapply(seq_len(19),function(i) score(pop_current[i,],future_candidate[n_grid+i,,drop=FALSE],FALSE),numeric(1))
  results <- list();ri <- 0L
  for(i in seq_len(19)) for(direction in c("forward","reverse")) for(penalty in c(TRUE,FALSE)) {
    first <- if(direction=="forward") pop_current[i,] else future_candidate[n_grid+i,,drop=FALSE]
    targets <- if(direction=="forward") future_candidate else candidate
    minimum <- rep(Inf,length(radii));best <- rep(NA_integer_,length(radii));ties <- rep(0L,length(radii));nearest <- rep(Inf,length(radii))
    for(start in seq.int(1,n_candidate,by=25000L)) {
      end <- min(n_candidate,start+24999L);idx <- start:end
      scores <- score(first,targets[idx,,drop=FALSE],penalty)
      for(radius in seq_along(radii)) {
        eligible <- distances[[i]][idx]<=radii[radius]+1e-10
        if(!any(eligible)) next
        m <- min(scores[eligible])
        if(m<minimum[radius]-1e-10) {
          minimum[radius] <- m;ties[radius] <- 0L;best[radius] <- NA_integer_;nearest[radius] <- Inf
        }
        if(abs(m-minimum[radius])<=1e-10) {
          tied <- idx[eligible & abs(scores-minimum[radius])<=1e-10]
          ties[radius] <- ties[radius]+length(tied)
          win <- tied[order(distances[[i]][tied],tied)][1]
          if(distances[[i]][win]<nearest[radius]-1e-10 || (abs(distances[[i]][win]-nearest[radius])<=1e-10 && (is.na(best[radius]) || win<best[radius]))) {
            best[radius] <- win;nearest[radius] <- distances[[i]][win]
          }
        }
      }
    }
    stopifnot(all(is.finite(minimum)),all(!is.na(best)),all(minimum<=local[i]+1e-8))
    for(radius in seq_along(radii)) {
      j <- best[radius]
      ri <- ri+1L
      results[[ri]] <- data.frame(population=pop_meta$population[i],group=pop_meta$group[i],ssp=ssp,period=period,gcm=gcm,direction=direction,
        search_scoring=if(penalty) "inherited_geographic_penalty_degrees" else "environment_only",radius_km=if(is.infinite(radii[radius])) "unrestricted" else as.character(radii[radius]),
        offset=minimum[radius],local_offset=local[i],distance_km=nearest[radius],tie_count=ties[radius],
        endpoint_index=j,endpoint_longitude=targets$longitude[j],endpoint_latitude=targets$latitude[j],
        endpoint_population=if(j>n_grid) pop_meta$population[j-n_grid] else "",
        endpoint_population_class=if(j>n_grid) pop_meta$group[j-n_grid] else "",candidate_domain="fixed_current_habitat_grid_plus_19_verified_population_coordinates")
    }
    message(name," ",pop_meta$population[i]," ",direction," ",if(penalty) "geographic" else "environment")
  }
  result <- do.call(rbind,results)
  write_csv(result,csv_path)
  jsonlite::write_json(list(validated=TRUE,rows=nrow(result),candidate_count=n_candidate,domain=unique(result$candidate_domain),chunk_size=25000),checkpoint,pretty=TRUE,auto_unbox=TRUE)
  ai <- ai+1L;all_results[[ai]] <- result
}
all_results <- do.call(rbind,all_results)
library(sf)
points <- st_as_sf(all_results,coords=c("endpoint_longitude","endpoint_latitude"),crs=4326,remove=FALSE)
membership <- function(filename) {
  polygon <- st_read(file.path(bundle,"ammo","gf",filename),quiet=TRUE)
  if(is.na(st_crs(polygon))) st_crs(polygon) <- 4326
  projected_polygon <- st_transform(polygon,32648)
  valid_polygon <- st_make_valid(projected_polygon)
  lengths(st_intersects(st_transform(points,32648),valid_polygon))>0
}
incore <- membership("ammo_core.shp");inperi <- membership("ammo_peri.shp")
all_results$endpoint_habitat_class <- ifelse(incore & inperi,"overlap",ifelse(incore,"core",ifelse(inperi,"peripheral","other")))
write_csv(all_results,file.path(directory,"full_habitat_domain_population_matches.csv"))
summary <- aggregate(list(n=all_results$offset),all_results[,c("group","ssp","period","gcm","direction","search_scoring","radius_km","endpoint_habitat_class")],length)
summary$denominator <- ifelse(summary$group=="core",9,10)
summary$percentage_of_focal_populations <- 100*summary$n/summary$denominator
write_csv(summary,file.path(directory,"source_destination_population_percentages.csv"))
save_session(directory)
marker(directory,list(rows=nrow(all_results),candidate_grid_cells=n_grid,known_population_points_added=19,
  scope="Full fixed-habitat candidate search for 19 biological focal populations; no all-pixel pseudoreplicated inference. Landscape-wide all-focal maps remain separate.",
  native_climate_handling="Point extraction at identical fixed-grid and known-population coordinates; future habitat contraction not silently applied as a domain change",
  geographic_model_note="Inherited geo_TRUE fit uses Euclidean longitude/latitude degrees; degree-distance penalty retained only as explicitly labeled inherited specification; environmental-only and distance/domain sensitivities reported. Metric-projected refit still pending.",
  illustrative_radii_km=c(25,50,100),radius_caveat="Illustrative limits, not empirically estimated species dispersal distances",range_crs_note=range_crs_note,
  habitat_geometry_handling="Archived rings contain duplicate vertices; transform to EPSG:32648 and GEOS st_make_valid before membership. Original polygons remain unchanged.",
  direction="Forward endpoints are future destinations; reverse endpoints are present candidate analogues/sources",percentage_denominator="Nine or ten focal populations per class, not source area or inferred migration probability"))
