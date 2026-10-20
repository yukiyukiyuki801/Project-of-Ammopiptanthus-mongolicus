#!/usr/bin/env Rscript
# Exact basis: R3.M1; user requests population segregating-site ratios and dosage alternatives.
source('scripts/revision/common.R')
suppressPackageStartupMessages({library(jsonlite)})
directory <- file.path(output_dir,'functional_reconstruction/postprocess')
stopifnot(fromJSON(file.path(directory,'score_parsing_validation.json'))$validated)
s <- read.csv(file.path(directory,'population_variant_counts.csv'));ind <- read.csv(file.path(directory,'individual_functional_counts.csv'))
stopifnot(nrow(ind)==217,length(unique(ind$population))==19)
set.seed(20261014)
# Population site inclusion is >=80% genotyped individuals. Report excluded sites;
# this criterion cannot recover unavailable invariant sites or mapping quality.
s$qualified <- s$called_individuals>=ceiling(.8*s$population_individuals)
s$common_qualified <- ave(as.integer(s$qualified),s$index,FUN=function(x)as.integer(length(x)==19 && all(x==1)))==1
stopifnot(any(s$qualified));m <- min(s$called_chromosomes[s$qualified]);stopifnot(m>=2)
classes <- c('synonymous','strict_LOF','SIFT_deleterious','SIFT_tolerated','SIFT_deleterious_strict_quality','SIFT_tolerated_strict_quality','missense_unscored','missense_all')
p0 <- function(a,n,m) ifelse(a<m,0,exp(lchoose(a,m)-lchoose(n,m)))
rows <- list();k<-0
for(pop in sort(unique(s$population)))for(cls in classes){
 x<-s[s$population==pop & vapply(strsplit(s$class_labels,';',fixed=TRUE),function(v)cls%in%v,logical(1)),]
 q<-x[x$qualified,];N<-q$called_chromosomes;A<-q$ALT_copies
 poly<-if(nrow(q))sum(A>0&A<N)else 0
 rare<-if(nrow(q))sum(pmax(0,1-p0(N-A,N,m)-p0(A,N,m)))else 0
 k<-k+1;rows[[k]]<-data.frame(population=pop,group=unique(s$group[s$population==pop]),class=cls,total_sites=nrow(x),qualified_sites=nrow(q),excluded_low_call_sites=nrow(x)-nrow(q),segregating_sites=poly,population_fixed_ALT_sites=sum(A==N & N>0),population_fixed_REF_sites=sum(A==0 & N>0),expected_segregating_sites_common_m=rare,common_called_chromosomes=m,common_mask_sites=sum(x$common_qualified),ALT_frequency_sum_common_mask=sum(x$ALT_copies[x$common_qualified]/x$called_chromosomes[x$common_qualified]))
}
counts<-do.call(rbind,rows);write_csv(counts,file.path(directory,'population_site_composition.csv'))
ratio<-function(a,b)ifelse(is.finite(b)&b>0,a/b,NA_real_)
numclasses<-c('strict_LOF','SIFT_deleterious','SIFT_deleterious_strict_quality');metrics<-list();k<-0
for(pop in sort(unique(ind$population))){
 x<-ind[ind$population==pop,];syn<-counts[counts$population==pop&counts$class=='synonymous',]
 for(cls in numclasses){
  y<-counts[counts$population==pop&counts$class==cls,]
  vals<-c(summed_ALT_frequency_ratio_common_mask=ratio(y$ALT_frequency_sum_common_mask,syn$ALT_frequency_sum_common_mask),segregating_variant_ratio=ratio(y$segregating_sites,syn$segregating_sites),rarefied_segregating_variant_ratio=ratio(y$expected_segregating_sites_common_m,syn$expected_segregating_sites_common_m),pooled_ALT_dosage_ratio=ratio(sum(x[[paste0(cls,'_ALT_dosage')]]),sum(x$synonymous_ALT_dosage)),mean_individual_ALT_dosage_ratio=mean(ratio(x[[paste0(cls,'_ALT_dosage')]],x$synonymous_ALT_dosage),na.rm=TRUE),pooled_homozygous_ALT_ratio=ratio(sum(x[[paste0(cls,'_homozygous_ALT')]]),sum(x$synonymous_homozygous_ALT)),ALT_dosage_per_synonymous_allele_opportunity=ratio(sum(x[[paste0(cls,'_ALT_dosage')]]),2*sum(x$synonymous_called_sites)))
  for(metric in names(vals)){k<-k+1;metrics[[k]]<-data.frame(population=pop,group=unique(x$group),numerator_class=cls,metric=metric,value=unname(vals[metric]),sampled_individuals=nrow(x),synonymous_segregating_sites=syn$segregating_sites,synonymous_expected_segregating_sites=syn$expected_segregating_sites_common_m,synonymous_ALT_copies=sum(x$synonymous_ALT_dosage),synonymous_homozygous_ALT=sum(x$synonymous_homozygous_ALT),synonymous_called_allele_opportunities=2*sum(x$synonymous_called_sites))}
 }
}
metrics<-do.call(rbind,metrics);write_csv(metrics,file.path(directory,'population_ratios.csv'))
# Each of 19 populations is ONE inferential unit, irrespective of individual n.
res<-list();k<-0
for(cls in numclasses)for(metric in unique(metrics$metric)){
 x<-metrics[metrics$numerator_class==cls&metrics$metric==metric,];a<-x$value[x$group=='core'];b<-x$value[x$group=='peripheral'];stopifnot(length(a)==9,length(b)==10)
 if(any(!is.finite(c(a,b)))){p<-NA_real_;ci<-c(NA,NA);boot<-c(NA,NA)}else{
  p<-tryCatch(t.test(b,a,var.equal=FALSE)$p.value,error=function(e)NA_real_)
  boot<-replicate(2000,mean(sample(b,replace=TRUE))-mean(sample(a,replace=TRUE)));ci<-quantile(boot,c(.025,.975))
 }
 k<-k+1;res[[k]]<-data.frame(numerator_class=cls,metric=metric,n_core=length(a),n_peripheral=length(b),mean_core=mean(a),mean_peripheral=mean(b),difference_peripheral_minus_core=mean(b)-mean(a),bootstrap_population_CI_lower=ci[1],bootstrap_population_CI_upper=ci[2],Welch_p=p)
}
res<-do.call(rbind,res);res$BH_q<-p.adjust(res$Welch_p,'BH');write_csv(res,file.path(directory,'population_class_comparisons.csv'))
coverage<-do.call(rbind,lapply(split(counts,counts$population),function(x){get<-function(cls,col)x[x$class==cls,col];data.frame(population=x$population[1],group=x$group[1],qualified_missense_sites=get('missense_all','qualified_sites'),scored_missense_sites=sum(get('SIFT_deleterious','qualified_sites'),get('SIFT_tolerated','qualified_sites')),unscored_missense_sites=get('missense_unscored','qualified_sites'),segregating_scored_missense_sites=sum(get('SIFT_deleterious','segregating_sites'),get('SIFT_tolerated','segregating_sites')),DEL_fraction_among_scored_segregating_missense=ratio(get('SIFT_deleterious','segregating_sites'),sum(get('SIFT_deleterious','segregating_sites'),get('SIFT_tolerated','segregating_sites'))),strict_quality_DEL_fraction_among_scored_segregating_missense=ratio(get('SIFT_deleterious_strict_quality','segregating_sites'),sum(get('SIFT_deleterious_strict_quality','segregating_sites'),get('SIFT_tolerated_strict_quality','segregating_sites'))))}))
write_csv(coverage,file.path(directory,'population_SIFT_coverage.csv'))
writeLines(c('All ratios describe annotated variants in the supplied variant-only RAD panel, not fitness load, genome-wide mutation burden or derived alleles.',
'Population-segregating LOF/SYN and DEL/SYN count sites with 0<ALT copies<called chromosomes and at least80% called individuals; fixed REF and fixed ALT excluded only from these composition ratios. Numerator and denominator counts are reported. Unequal sample sizes affect observed segregating counts.',
sprintf('Rarefied ratios use exact hypergeometric expected segregating-site counts at common m=%d called chromosomes, the minimum over qualified population-site observations. No random subsampling, genotype imputation or model refitting. E[S]=sum(1-C(N-A,m)/C(N,m)-C(A,m)/C(N,m)). Low-call exclusions and rarefaction size are recorded.',m),
'Primary ALT-frequency ratio sums per-site ALT frequencies within each consequence class, then divides by the analogous synonymous sum. Uses a common site mask requiring at least80% called individuals in ALL19 populations, correcting per-site call-count variation and keeping the site set consistent. Retains fixed ALT. This remains reference-relative, not ancestral.',
'Pooled ALT-dosage ratio sums numerator ALT copies divided by synonymous ALT copies across called genotypes in each population; retains fixed ALT sites. Mean individual ALT-dosage ratio averages individual ratios within population; it is distinct from pooled ratios. Homozygous-ALT ratios likewise retain fixed ALT. Opportunity-normalized dosage uses twice the observed called synonymous-site count, not a callable-genome denominator.',
'SIFT all-scored emitted labels are primary; missing predictions excluded and reported, not assigned tolerated. Strict-quality sensitivity requires reported median sequence information<=3.25, at least2 sequences at position, and no position/reference warning. This is an explicitly defined quality sensitivity, not a recovered historical setting. Printed scores are rounded to2decimals; compiled engine classifies raw<0.05 deleterious and equality tolerated, so emitted labels determine classification.',
'Class comparisons use9core versus10peripheral population values, Welch tests,2000 independent population bootstrap draws for peripheral-minus-core difference, and BH adjustment across21 explicitly reported class/metric comparisons. Population exchangeability and spatial independence remain limitations; no individual-level pseudoreplication. Ratios and tests are conditional on this panel, annotation and SIFT scoring coverage.',
'Legacy plotting code labels imported ratios Deleterious sites/CDS sites and LoF sites/CDS sites, whereas published2025Methods2.6 specifies synonymous variants. The legacy code imports output.csv and does not calculate the numerator or denominator, so this inconsistency cannot establish which denominator was historically executed. New ratios explicitly use synonymous-site or synonymous-ALT denominators and are not labeled exact historical reproduction.',
'Published2025Methods2.6 describes individual derived LOF/DEL-to-synonymous proxies usingUniRef100; does not define population polymorphic-site exclusion. New calculations are not exact reproduction: ancestral polarization unavailable and current SIFT reference database is SwissProt. New variant-composition ratios and ALT-relative dosage summaries must be named explicitly.'),file.path(directory,'METHODS_AND_LIMITS.md'))
save_session(directory);marker(directory,list(scope='R3.M1 and explicit user requested population ratio alternatives',populations=19,individuals=217,common_rarefaction_chromosomes=m,comparison_tests=nrow(res),bootstrap_replicates=2000,ratio_rows=nrow(metrics),figure_status='Pending R saved-result rendering'))
