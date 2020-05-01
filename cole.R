library(aqp)
library(soilDB)
library(magrittr)

# justification for use of Vertic subgroup in 7083b, 7085b 

# get data in MLRA 15 and 18 - coast range and sierra nevada foothills
f <- fetchKSSL(mlra=c("15","18"))

# optional: filter by taxonomy
#f.sub <- filter(f, taxpartsize=="fine" | taxpartsize=="clayey")

# or use whole set
f.sub <- f

# sqrt transform coefficient of linear extensibility (determined good transformation by inspection)
f.sub$COLEws_sqrt <- sqrt(f.sub$COLEws)

# create linear model with transformed dependent var
fm <- formula(COLEws_sqrt ~ clay)
m <- lm(fm, data=horizons(f.sub))

# make predictions using model coefficeints and 80% prediction interval
clayz <- data.frame(clay=f.sub$clay, COLEws=f.sub$COLEws)
clayz <- cbind(clayz, as.data.frame(predict(m, newdata=clayz, interval = "prediction", level = 0.8)))
clayz <- clayz[order(clayz$clay),]
clayz[,c('fit','lwr','upr')] <- apply(clayz[,c('fit','lwr','upr')], 2, function(x) x^2)

# inspect model
plot(COLEws ~ clay, data=clayz, pch=".", cex=2)
lines(fit ~ clay, data=clayz, col="RED")
lines(lwr ~ clay, data=clayz, col="BLUE", lty=2)
lines(upr ~ clay, data=clayz, col="BLUE", lty=2)

f.sub$cole_pred <- clayz$fit

# get NASIS pedons
nasis <- fetchNASIS()

# merge predictions from model for NASIS data into horizon slot
horizons(nasis) <- cbind(horizons(nasis)[,c("peiid", "phiid")],
                         predict(m, newdata=horizons(nasis),interval = "prediction", level = 0.8))

# take upper 100cm, or to contact, and sum up COLE
f.new <- nasis %>%
  glomApply(function(x) c(0, min(estimateSoilDepth(x), 100), na.rm=T), truncate=TRUE) %>% 
  mutate(thickness = hzdepb - hzdept, 
         margin_cole_pred = lwr*thickness,
         margin_clay = clay*thickness) %>%
  mutate_profile(thick_sum =  sum(thickness, na.rm=T),
                 cole_pred_sum = sum(margin_cole_pred, na.rm=T),
                 clay_wtd_avg = sum(margin_clay, na.rm=T) / thick_sum)

# take anything that makes 6cm or more by sum and has total thickness<75
# this is to get rid of some very deep soils with heavy clay substrata
# honing in on the shallow-dunstone to mod-deep argonaut concept
vertic <- f.new %>% 
  filter(cole_pred_sum >= 6, thick_sum < 75)

# make a grouped profile plot
vertic$taxonname <- factor(vertic$taxonname)
groupedProfilePlot(vertic, groups="taxonname", label="pedon_id")

# tabulate: how many in each taxonname?
res <- table(vertic$taxonname)
res <- res[rev(order(res))]
res

# are heavier / higher COLE soils on lower slopes?
plot(density(filter(f.new, cole_pred_sum < 6)$slope_field, 
             na.rm=T, from=0, bw = pi), ylim=c(0,0.05))
lines(density(vertic$slope_field, na.rm=T, from=0, bw = pi), col="red")

# distributions are similar, considering number of observations (except for steep end?)
quantile(f.new$slope_field, na.rm=T, probs=c(0,0.01,0.05,0.1,0.25,0.5,0.75,0.9,0.95,0.99,1))
quantile(vertic$slope_field, na.rm=T, probs=c(0,0.01,0.05,0.1,0.25,0.5,0.75,0.9,0.95,0.99,1))
