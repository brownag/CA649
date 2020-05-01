# Implementation of checks for Technical Note 18, Soil Survey Region 2
# author: andrew brown

f <- fetchNASIS('components')

coiid.hz.problem <- get('component.hz.problems', envir=soilDB.env)
if(length(coiid.hz.problem))
   stop("You have components with horizon problems: ", paste0(coiid.hz.problem, collapse=","))

# parse restrictions by component
res <- slot(f, 'restrictions')
coiid.res <- split(res, res$coiid)

## VALIDATE RESTRICTIONS
# check that all restriction top depths match a horizon top depth
all(as.logical(lapply(coiid.res, function(r) {
  all(r$resdept_r %in% filter(f, coiid == unique(r$coiid))$hzdept_r)
})))

# check that all restriction top depths match a diagnostic top depth
all(as.logical(lapply(coiid.res, function(r) {
  all(r$resdept_r %in% slot(filter(f, coiid == unique(r$coiid)), 'diagnostic')$featdept_r)
})))
