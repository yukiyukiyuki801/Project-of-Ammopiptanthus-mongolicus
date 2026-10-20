#!/usr/bin/env Rscript
# Targeted recovery check, not a replacement GDM fit. See README_scope.md.
source("scripts/revision/common.R")
directory <- file.path(output_dir, "fst_subset_comparison")
dir.create(directory, recursive=TRUE, showWarnings=FALSE)
stopifnot(nrow(sample_meta)==217L, nrow(pop_meta)==19L,
          identical(sample_meta$vcf_index, seq_len(217L)))
input <- file.path(output_dir,"figures","genetic_pca","full_genotypes_conversion.tsv")
g <- matrix(scan(input, what=integer(), quiet=TRUE), nrow=217L, byrow=TRUE)
stopifnot(identical(dim(g),c(217L,38124L)), all(g %in% c(0L,1L,2L,9L)))
g[g==9L] <- NA_integer_
loci <- read.csv(file.path(output_dir,"inputs","full_loci.csv"))
adaptive <- read.csv(file.path(output_dir,"inputs","adaptive_loci.csv"))
stopifnot(identical(loci$locus_index,seq_len(ncol(g))),nrow(adaptive)==1637L,
          !anyDuplicated(adaptive$full_locus_index))
candidate <- adaptive$full_locus_index
stopifnot(all(loci$chrom[candidate]==adaptive$chrom),
          all(loci$pos[candidate]==adaptive$pos),
          all(loci$ref[candidate]==adaptive$ref),all(loci$alt[candidate]==adaptive$alt))
sets <- list(all_SNPs=seq_len(ncol(g)),candidate=candidate,
             noncandidate=setdiff(seq_len(ncol(g)),candidate))
stopifnot(identical(vapply(sets,length,integer(1)),
                    c(all_SNPs=38124L,candidate=1637L,noncandidate=36487L)))
pops <- pop_meta$population
N <- ALT <- HET <- matrix(NA_real_,nrow=length(pops),ncol=ncol(g),dimnames=list(pops,NULL))
for(i in seq_along(pops)) {
  z <- g[sample_meta$population==pops[i],,drop=FALSE]
  stopifnot(nrow(z)==pop_meta$n[i])
  N[i,] <- colSums(!is.na(z))
  ALT[i,] <- colSums(z,na.rm=TRUE)
  HET[i,] <- colSums(z==1L,na.rm=TRUE)
}
# Biallelic simplification of the unchanged official hierfstat::wc ANOVA:
# use one allele's components; the complementary allele duplicates all three
# components and therefore cancels from the multilocus ratio of sums.
wc_components <- function(n1,n2,x1,x2,h1,h2) {
  nt <- n1+n2
  p1 <- x1/(2*n1); p2 <- x2/(2*n2)
  pb <- (x1+x2)/(2*nt)
  nc <- nt-(n1^2+n2^2)/nt
  MSG <- (h1+h2)/(2*nt)
  MSI <- (2*n1*p1*(1-p1)+2*n2*p2*(1-p2)-(h1+h2)/2)/(nt-2)
  MSP <- 2*(n1*(p1-pb)^2+n2*(p2-pb)^2)
  a <- (MSP-MSI)/(2*nc)
  b <- (MSI-MSG)/2
  c <- MSG
  valid <- n1>0 & n2>0 & nt>2 & is.finite(a+b+c)
  list(a=a,b=b,c=c,valid=valid)
}
ratio <- function(v,idx) {
  keep <- idx[v$valid[idx]]
  denom <- sum(v$a[keep]+v$b[keep]+v$c[keep])
  stopifnot(length(keep)>0,is.finite(denom),denom>0)
  c(fst=sum(v$a[keep])/denom,usable=length(keep),excluded=length(idx)-length(keep))
}
# Oracle uses unchanged official sources pinned to a recorded upstream commit.
# Only is.genind's class guard is supplied locally because the input is a plain
# data.frame and adegenet is not installed; no wc algorithm is modified.
oracle <- new.env(parent=globalenv())
oracle$is.genind <- function(x) inherits(x,"genind")
for(f in c("misc.R","nwc.R","pairwise.fst.R"))
  sys.source(file.path(directory,"primary_sources",f),envir=oracle)
set.seed(20261015)
oracle_rows <- list()
oi <- 0L
for(pair in list(c(1L,2L),c(1L,19L),c(8L,12L))) {
  i <- pair[1]; j <- pair[2]
  v <- wc_components(N[i,],N[j,],ALT[i,],ALT[j,],HET[i,],HET[j,])
  for(tag in names(sets)) {
    eligible <- sets[[tag]][v$valid[sets[[tag]]]]
    idx <- sample(eligible,min(300L,length(eligible)))
    rows <- which(sample_meta$population %in% pops[c(i,j)])
    # hierfstat's diploid allele encoding 11=0/0,12=0/1,22=1/1.
    encoded <- matrix(c(11,12,22)[g[rows,idx,drop=FALSE]+1L],nrow=length(rows))
    dat <- data.frame(pop=match(sample_meta$population[rows],pops[c(i,j)]),encoded)
    expected <- oracle$pairwise.WCfst(dat)[1,2]
    actual <- unname(ratio(v,idx)["fst"])
    oi <- oi+1L
    oracle_rows[[oi]] <- data.frame(population_1=pops[i],population_2=pops[j],subset=tag,
                                   loci=length(idx),official_hierfstat=expected,
                                   simplified=actual,absolute_difference=abs(actual-expected))
  }
}
oracle_table <- do.call(rbind,oracle_rows)
write_csv(oracle_table,file.path(directory,"official_source_oracle_validation.csv"))
stopifnot(all(is.finite(oracle_table$absolute_difference)),
          max(oracle_table$absolute_difference)<1e-12)
arch <- read.csv(file.path(output_dir,"inputs","adaptive_fst.csv"),check.names=FALSE)
stopifnot(setequal(arch$population,pops),setequal(names(arch)[-1],pops),!anyDuplicated(arch$population))
arch <- as.matrix(arch[match(pops,arch$population),pops])
stopifnot(max(abs(arch-t(arch)))<1e-12,all(diag(arch)==0))
pair_rows <- list(); k <- 0L
matrices <- lapply(sets,function(x) matrix(0,length(pops),length(pops),dimnames=list(pops,pops)))
for(i in seq_len(length(pops)-1L)) for(j in (i+1L):length(pops)) {
  v <- wc_components(N[i,],N[j,],ALT[i,],ALT[j,],HET[i,],HET[j,])
  for(tag in names(sets)) {
    z <- ratio(v,sets[[tag]])
    matrices[[tag]][i,j] <- matrices[[tag]][j,i] <- z["fst"]
    k <- k+1L
    pair_rows[[k]] <- data.frame(subset=tag,population_1=pops[i],population_2=pops[j],
      loci_in_panel=length(sets[[tag]]),usable_loci=unname(z["usable"]),
      excluded_loci=unname(z["excluded"]),fst=unname(z["fst"]),
      archived_fst=arch[i,j],difference=unname(z["fst"])-arch[i,j])
  }
}
pairs <- do.call(rbind,pair_rows)
write_csv(pairs,file.path(directory,"pairwise_comparison.csv"))
for(tag in names(matrices)) write_csv(data.frame(population=pops,matrices[[tag]],check.names=FALSE),
                                         file.path(directory,paste0("fst_",tag,".csv")))
summaries <- do.call(rbind,lapply(names(sets),function(tag) {
  z <- pairs[pairs$subset==tag,]
  data.frame(subset=tag,loci=length(sets[[tag]]),pairs=nrow(z),
    max_absolute_difference=max(abs(z$difference)),rmse=sqrt(mean(z$difference^2)),
    mean_absolute_difference=mean(abs(z$difference)),mean_difference=mean(z$difference),
    pearson=cor(z$fst,z$archived_fst),matches_absolute_1e_12=sum(abs(z$difference)<1e-12),
    matches_absolute_1e_8=sum(abs(z$difference)<1e-8),
    matches_rounded_6dp=sum(round(z$fst,6)==round(z$archived_fst,6)),
    matches_rounded_4dp=sum(round(z$fst,4)==round(z$archived_fst,4)),
    matches_rounded_3dp=sum(round(z$fst,3)==round(z$archived_fst,3)),
    minimum_usable_loci=min(z$usable_loci),maximum_excluded_loci=max(z$excluded_loci))
}))
summaries <- summaries[order(summaries$rmse),]
write_csv(summaries,file.path(directory,"subset_match_summary.csv"))
writeLines(c("# Archived FST subset comparison", "",capture.output(print(summaries,row.names=FALSE)),"",
  "All comparisons align the 19 population labels and use all 171 unique pairs.",
  "Noncandidate means the complement of the candidate set; it is not a certified neutral panel.",
  "No negative FST truncation, genotype imputation or additional locus filtering was applied.",
  "A numerical closest match supports a subset hypothesis but does not establish historical execution/version.",
  "This job does not refit GDM or replace its archived input matrix."),file.path(directory,"RESULTS.md"))
save_session(directory)
marker(directory,list(scope="Explicitly authorized FST subset recovery comparison only; no GDM refit",
  samples=nrow(g),populations=length(pops),unique_pairs=171,
  panels=as.list(vapply(sets,length,integer(1))),
  implementation="Biallelic simplification of official hierfstat WC ANOVA, diploid ratio of summed variance components",
  oracle_commit="fbc8b8c34f62389213ca0a1b721533744b4735c1",
  oracle_checks=nrow(oracle_table),oracle_max_absolute_difference=max(oracle_table$absolute_difference),
  closest_rmse_panel=summaries$subset[1],historical_execution_proved=FALSE,
  missing_encoding_input=9,missing_handling="Pair/locus-specific called individuals; require both populations typed and total called individuals >2",
  source_manifest="primary_sources/SOURCE_MANIFEST.json"))
