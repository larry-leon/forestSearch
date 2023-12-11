
source_fs_functions<-function(file_loc="/R/"){
files_toload<-c("forestsearch_functions_v0.R","subgroup_search_v0.R","subgroup_consistency_v0.R","forest_search_v0.R",
"virtual-twins_functions_v0.R","grf_functions_v0.R",
"sim_aft_gbsg-mod4_v0.R","oc_analyses_m4_FS4_v0.R",
"forestsearch_bootstrap_functions_v0.R","forestsearch_bootstrap-parallel_v0.R","plotting_functions_v0.R","km_resampling_functions.R",
"get_FSdata_v0.R","bootstrap_summary_v0.R","subgroup_results_summary_v0.R")

for(ff in 1:length(files_toload)){
source(c(paste0(file_loc,files_toload[ff])))
}
}
      
#source_fs_functions(file_loc="~/leolarr2/KN-062/R/")