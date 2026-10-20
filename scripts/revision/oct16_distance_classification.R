source(file.path("scripts","revision","common.R"))
library(terra)
library(gdm)
library(geosphere)
library(sf)
directory <- file.path(output_dir,"oct16_revision","distance_classification")
dir.create(directory,recursive=TRUE,showWarnings=FALSE)
writeLines(c("Explicit authorization: uploaded author response comments16and18,16October2026. Exact text archived revision_oct16/analysis/REQUEST_BASIS_DISTANCE_CLASSIFICATION.md.","Fixed before results:10/200km missing matches only,19focals,4ensemblecases,2directions,2scorings. Existing50/100/unrestricted reused. Class sensitivity:9nearest occurrence-hull centroid versus10outer;12inside-core-polygon versus7outside. Saved primary functional indices,Ho,RONA,GDM; no model refits."),file.path(directory,"REQUEST_BASIS.txt"))

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
radii <- c(10,200)
gcms <- "ENSEMBLE"
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
new_results <- do.call(rbind,all_results)
old_path <- file.path(output_dir,"landscape","full_habitat_domain_population_matches.csv")
old <- read.csv(old_path)
old <- old[old$gcm=="ENSEMBLE" & old$radius_km %in% c("50","100","unrestricted"),]
stopifnot(nrow(old)==19*4*2*2*3)
new_results$endpoint_habitat_class <- NA_character_
points <- sf::st_as_sf(new_results,coords=c("endpoint_longitude","endpoint_latitude"),crs=4326,remove=FALSE)
membership <- function(filename) {
 polygon <- sf::st_read(file.path(bundle,"ammo","gf",filename),quiet=TRUE)
 if(is.na(sf::st_crs(polygon))) sf::st_crs(polygon) <- 4326
 lengths(sf::st_intersects(sf::st_transform(points,32648),sf::st_make_valid(sf::st_transform(polygon,32648))))>0
}
incore <- membership("ammo_core.shp");inperi <- membership("ammo_peri.shp")
new_results$endpoint_habitat_class <- ifelse(incore&inperi,"overlap",ifelse(incore,"core",ifelse(inperi,"peripheral","other")))
all_results <- rbind(old[,names(new_results)],new_results)
stopifnot(nrow(all_results)==19*4*2*2*5)
all_results$radius_order <- match(all_results$radius_km,c("10","50","100","200","unrestricted"))
all_results <- all_results[order(all_results$population,all_results$ssp,all_results$period,all_results$direction,all_results$search_scoring,all_results$radius_order),]
key <- interaction(all_results[,c("population","ssp","period","direction","search_scoring")],drop=TRUE)
for(k in levels(key)) {
 d <- all_results[key==k,]
 stopifnot(nrow(d)==5,all(diff(d$offset)<=1e-9),all(d$offset<=d$local_offset+1e-8))
 finite <- d$radius_km!="unrestricted"
 stopifnot(all(d$distance_km[finite]<=as.numeric(d$radius_km[finite])+1e-8))
}
# Re-evaluate saved winner scores for 50/100/unrestricted with exactly this model
# and extraction, without rerunning the already completed search domains.
checks <- list();ci <- 0L
for(ssp in c("ssp126","ssp585")) for(period in c("2041-2070","2071-2100")) {
 future <- read_stack(file.path(base,paste0(ssp,"_",if(period=="2041-2070") "2040" else "2070")))
 fv <- terra::extract(future,candidate[,c("longitude","latitude")],method="simple",ID=FALSE)
 future_candidate <- data.frame(candidate[,c("longitude","latitude")],fv)
 d <- old[old$ssp==ssp & old$period==period,]
 for(i in seq_len(nrow(d))) {
  focal <- match(d$population[i],pop_meta$population);j <- d$endpoint_index[i]
  first <- if(d$direction[i]=="forward") pop_current[focal,] else future_candidate[n_grid+focal,,drop=FALSE]
  targets <- if(d$direction[i]=="forward") future_candidate else candidate
  value <- score(first,targets[j,,drop=FALSE],d$search_scoring[i]=="inherited_geographic_penalty_degrees")
  ci <- ci+1L;checks[[ci]] <- data.frame(population=d$population[i],ssp=ssp,period=period,direction=d$direction[i],search_scoring=d$search_scoring[i],radius_km=d$radius_km[i],absolute_error=abs(value-d$offset[i]))
 }
}
checks <- do.call(rbind,checks);stopifnot(max(checks$absolute_error)<1e-10)
write_csv(checks,file.path(directory,"reused_winner_score_validation.csv"))
write_csv(all_results,file.path(directory,"all_requested_radius_population_matches.csv"))
# Complete zero-count categories so absent pie slices are explicitly recoverable.
groups <- c("core","peripheral");cats <- c("core","peripheral","other","overlap")
summary <- expand.grid(group=groups,ssp=c("ssp126","ssp585"),period=c("2041-2070","2071-2100"),gcm="ENSEMBLE",direction=c("forward","reverse"),search_scoring=unique(all_results$search_scoring),radius_km=c("10","50","100","200","unrestricted"),endpoint_habitat_class=cats,stringsAsFactors=FALSE)
summary$n <- vapply(seq_len(nrow(summary)),function(i) sum(Reduce(`&`,lapply(names(summary),function(n)all_results[[n]]==summary[[n]][i]))),integer(1))
summary$denominator <- ifelse(summary$group=="core",9,10)
summary$percentage_of_focal_populations <- 100*summary$n/summary$denominator
stopifnot(all(abs(aggregate(percentage_of_focal_populations~group+ssp+period+direction+search_scoring+radius_km,summary,sum)$percentage_of_focal_populations-100)<1e-10))
write_csv(summary,file.path(directory,"requested_radius_endpoint_percentages.csv"))

# New, prespecified alternative classes; no outcome is consulted in assignment.
occ <- read.csv(file.path(root,"recovery_audit","published_2022_occurrences.csv"))
stopifnot(nrow(occ)==37)
occ_xy <- st_transform(st_as_sf(occ,coords=c("longitude","latitude"),crs=4326),32648)
hull <- st_convex_hull(st_union(occ_xy))
centre <- st_centroid(hull)
pop_xy <- st_transform(st_as_sf(pop_meta,coords=c("longitude","latitude"),crs=4326,remove=FALSE),32648)
centre_distance <- as.numeric(st_distance(pop_xy,centre))/1000
rank_order <- order(centre_distance,pop_meta$population)
geographic_group <- rep("peripheral",19);geographic_group[rank_order[seq_len(9)]] <- "core"
old_classes <- read.csv(file.path(output_dir,"classification","population_class_and_geographic_marginality.csv"))
old_classes <- old_classes[match(pop_meta$population,old_classes$population),]
stopifnot(identical(old_classes$population,pop_meta$population))
polygon_group <- ifelse(old_classes$habitat_polygon_class %in% c("core","overlap"),"core","peripheral")
classes <- data.frame(pop_meta,published_group=pop_meta$group,habitat_polygon_class=old_classes$habitat_polygon_class,polygon_core_vs_outside_group=polygon_group,centroid_distance_km=centre_distance,centroid_distance_rank=match(seq_len(19),rank_order),geographic_inner9_vs_outer10_group=geographic_group)
centre_ll <- st_coordinates(st_transform(centre,4326))[1,]
classes$centre_longitude <- centre_ll[1];classes$centre_latitude <- centre_ll[2]
write_csv(classes,file.path(directory,"alternative_population_classifications.csv"))
st_write(st_sf(definition="Area centroid of occurrence convex hull",geometry=centre),file.path(directory,"geographic_centre.gpkg"),delete_dsn=TRUE,quiet=TRUE)
# Collect measurements already computed; no models or annotation are refit.
measurements <- list();mi <- 0L
add_values <- function(d,val,analysis,metric,ssp="",period="",scoring="",family) {
 mi <<- mi+1L
 measurements[[mi]] <<- data.frame(population=d$population,analysis=analysis,metric=metric,ssp=ssp,period=period,search_scoring=scoring,value=val,test_family=family)
}
rona <- read.csv(file.path(output_dir,"rona","rona_population_long.csv"));rona <- rona[rona$gcm=="ENSEMBLE",]
stopifnot(all(rona$association_filter=="inherited_p_lt_1"),all(rona$K==4))
for(ssp in c("ssp126","ssp585")) for(period in c("2041-2070","2071-2100")) for(bio in bio_names) {
 d <- rona[rona$ssp==ssp & rona$period==period & rona$bio==bio,]
 stopifnot(nrow(d)==19)
 add_values(d,d$rona,"RONA",bio,ssp,period,family="RONA_20_ensemble_tests")
}
gdm <- read.csv(file.path(output_dir,"landscape_statistics","population_offsets_wide.csv"));gdm <- gdm[gdm$gcm=="ENSEMBLE",]
for(sc in unique(gdm$search_scoring)) for(ssp in c("ssp126","ssp585")) for(period in c("2041-2070","2071-2100")) for(metric in c("local","forward","reverse")) {
 d <- gdm[gdm$search_scoring==sc & gdm$ssp==ssp & gdm$period==period,];stopifnot(nrow(d)==19)
 add_values(d,d[[metric]],"GDM_full_habitat_unrestricted",metric,ssp,period,sc,family="GDM_24_ensemble_tests_both_scorings")
}
functional <- read.csv(file.path(output_dir,"functional_reconstruction","postprocess","population_ratios.csv"))
functional <- functional[functional$metric=="summed_ALT_frequency_ratio_common_mask" & functional$numerator_class %in% c("strict_LOF","SIFT_deleterious"),]
for(cl in unique(functional$numerator_class)) {
 d <- functional[functional$numerator_class==cl,];stopifnot(nrow(d)==19)
 add_values(d,d$value,"reference_relative_functional_index",paste0(cl,"_SYN_common_mask_frequency_ratio"),family="Functional_2_primary_ratios")
}
diversity <- read.csv(file.path(output_dir,"statistics","population_diversity.csv"));diversity <- diversity[diversity$subset=="full",];stopifnot(nrow(diversity)==19)
add_values(diversity,diversity$Ho,"diversity","Ho_full38124_panel",family="Ho_1_test")
measurements <- do.call(rbind,measurements)
write_csv(measurements,file.path(directory,"saved_population_measurements_for_classification.csv"))
set.seed(20261016)
B <- 2000L;effects <- list();ei <- 0L
class_defs <- c(published="published_group",polygon_core_vs_outside="polygon_core_vs_outside_group",geographic_inner9_vs_outer10="geographic_inner9_vs_outer10_group")
measure_keys <- interaction(measurements[,c("analysis","metric","ssp","period","search_scoring","test_family")],drop=TRUE)
for(def in names(class_defs)) for(k in levels(measure_keys)) {
 d <- measurements[measure_keys==k,];stopifnot(nrow(d)==19,!anyDuplicated(d$population),all(is.finite(d$value)))
 group <- classes[[class_defs[[def]]]][match(d$population,classes$population)]
 a <- d$value[group=="core"];b <- d$value[group=="peripheral"]
 # Polygon "peripheral" here is deliberately a comparison-group coding;
 # display labels must read outside-core polygon, not peripheral habitat.
 test <- suppressWarnings(wilcox.test(b,a,exact=FALSE,correct=TRUE))
 ci <- unname(quantile(replicate(B,mean(sample(b,length(b),TRUE))-mean(sample(a,length(a),TRUE))),c(.025,.975)))
 ei <- ei+1L
 effects[[ei]] <- data.frame(classification=def,d[1,c("analysis","metric","ssp","period","search_scoring","test_family")],n_core=length(a),n_comparison=length(b),mean_core=mean(a),mean_comparison=mean(b),comparison_minus_core=mean(b)-mean(a),CI_low=ci[1],CI_high=ci[2],rank_biserial=2*as.numeric(test$statistic)/(length(a)*length(b))-1,p_rank_sum=test$p.value)
}
effects <- do.call(rbind,effects);effects$q_BH <- NA_real_
family <- interaction(effects$classification,effects$test_family,drop=TRUE)
for(f in levels(family)){idx <- family==f;effects$q_BH[idx] <- p.adjust(effects$p_rank_sum[idx],"BH")}
write_csv(effects,file.path(directory,"alternative_classification_population_comparisons.csv"))
save_session(directory)
marker(directory,list(scope="Uploaded author comments16/18; missing10/200km search, existing50/100/unrestricted reuse, saved-measurement class sensitivity",missing_radius_match_rows=nrow(new_results),all_requested_radius_rows=nrow(all_results),endpoint_summary_rows=nrow(summary),classification_effect_rows=nrow(effects),candidate_grid_cells=n_grid,added_population_coordinates=19,focal_populations=19,ENSEMBLE_cases=4,new_radii_km=c(10,200),reused_radii_km=c(50,100,"unrestricted"),new_model_fits=0,maximum_reused_winner_score_error=max(checks$absolute_error),geographic_class_definition="9nearest sampled populations to area centroid of convex hull of37published occurrences inEPSG32648;10outer; lexical-ID tie break; sample-relative sensitivity",centre_longitude=centre_ll[1],centre_latitude=centre_ll[2],polygon_comparison="Inside core polygon(12) versus outside core polygon(7); outside includes4outside both polygons, not equivalent to peripheral habitat",bootstrap_repetitions=B,seed=20261016,test="Unpaired population Wilcoxon rank-sum with continuity correction; mean difference population bootstrap interval",BH_families="Within each classification:20RONA ensemble tests,24GDM ensemble tests(both scorings),2primary functional ratios,1Ho test",limitations="Reference-relative functional indices; illustrative distances not measured dispersal; occurrence-centre classification not true full-range delineation; exchangeable population tests do not remove spatial dependence"))
