
#file_loc<-c("/Users/larryleon/Library/CloudStorage/OneDrive-Personal/Documents/Myfiles/Projects-Current/subgroup identification and analysis/R/")
#source("//synology_nas/My Projects/Projects-Current/subgroup identification and analysis/R/source_forestsearch_v0.R")

source_fs_functions<-function(file_loc="../../R/"){
files_toload<-c("forestsearch_functions_v0.R","subgroup_search_v0.R","subgroup_consistency_v0.R","forest_search_v0.R",
                "virtual-twins_functions_v0.R","grf_functions_v0.R",
                "sim_aft_gbsg-mod4_v0.R","oc_analyses_m4_v0.R","plotting_functions_v0.R",
                "forestsearch_bootstrap_v0.R","forestsearch_bootstrap-parallel_v0.R")

for(ff in 1:length(files_toload)){
source(c(paste0(file_loc,files_toload[ff])))
}
}
      