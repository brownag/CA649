#workflowr.R

library(workflowr)

# wflow_git_config(user.name = "Andrew Gene Brown", user.email = "andrew.g.brown@usda.gov")
setwd("workflow")

wflow_start(".", git = FALSE, existing = TRUE)

wflow_open("analysis/00-get_project.Rmd")
wflow_open("analysis/01-project_extent.Rmd")
wflow_open("analysis/02-project_linework.Rmd")
wflow_open("analysis/03-project_terrain.Rmd")
wflow_open("analysis/04-project_spectral.Rmd")
wflow_open("analysis/05-project_pedons.Rmd")
wflow_open("analysis/06-project_geomorph.Rmd")

wflow_build()
