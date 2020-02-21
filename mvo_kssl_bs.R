#mvo_kssl_bs.R
library(aqp)
library(soilDB)
par(mfrow=c(1,1), mar=c(4,4,2,2))

# series <- c("auburn","bonanza","dunstone","sunnyslope",
#             "jasperpeak", "sobrante","loafercreek",
#             "argonaut","las posas","motherlode", "trabuco",
#             "boomer")
# f <- aqp::union(lapply(series, soilDB::fetchKSSL))
my.mlra <- "18"
my.bs.var <- "bs82"
my.analyte <- "Sum of Cations Base Saturation (pH 8.2)"
my.threshold <- 75

# my.mlra <- "22A"
# my.bs.var <- "bs7"
# my.analyte <- "Ammonium Acetate Base Saturation (pH 7.0)"
# my.threshold <- 60

f <- soilDB::fetchKSSL(mlra = my.mlra)
plot(horizons(f)[[my.bs.var]] ~ f$ph_h2o)
abline(h = my.threshold)
abline(v = 6)

aqp::site(f) <- aqp::getSoilDepthClass(f,  name =  "hzn_desgn", top = "hzn_top", bottom="hzn_bot", p = "Cr|R|Cd|^C\\d?$")

f$depthcl <- factor(aqp::denormalize(f, 'depth.class'))
argi.hz <- aqp::horizons(f)[0,]

res <- aqp::profileApply(f, function(p) {
  argi <- aqp::getArgillicBounds(p, hzdesgn = "hzn_desgn", texcl.attr = 'lab_texture_class')
  dept <- aqp::estimateSoilDepth(p, name =  "hzn_desgn", top = "hzn_top", bottom="hzn_bot")
  
  # if it has an argillic
  if(!any(is.na(argi))) {
    # one or more subhorizons within upper 75cm of argillic, or above contact, whichever is shallower
    crit.depth.haploxeralfs <- min(c(argi[1] + 75, dept))    
    hz <- horizons(aqp::glom(p, argi[1], crit.depth.haploxeralfs))
    argi.hz <<- rbind(argi.hz, hz)
    #NB: "any" used for haploxeralfs, but "all" used for palexeralfs
    return(any(hz[[my.bs.var]] < my.threshold))
  }
  return(-1)
})
res <- res[res != -1]
# 0 = no argillic; 1 = argillic
table(na.omit(res))

# remove values that are completely saturated with bases (no carbonates, in general)
argi.hz <- argi.hz[which(argi.hz[[my.bs.var]] < 100),]
plot(argi.hz[[my.bs.var]] ~ argi.hz$ph_h2o)
abline(v = 6)
abline(h=75)

argi.hz$meets_thresh <- (argi.hz[[my.bs.var]] > my.threshold) 
points(argi.hz[argi.hz$meets_thresh, my.bs.var] ~ argi.hz[argi.hz$meets_thresh, "ph_h2o"], pch=19, col="BLUE")
argi.hz$pedon_key <- factor(argi.hz$pedon_key)

# fit glm model, considering pH 1:1 water, horizon bottom depth, and conditioned on depth class
mylogit <- glm(meets_thresh ~ ph_h2o + hzn_bot + factor(depthcl), data = argi.hz, family = "binomial")
summary(mylogit)

png(filename = "foo.png", width = 500 ,height = 500)
par(mfrow=c(2, 2), oma = c(4, 0, 0, 0), mar=c(4,4,2,2))
depth.classes <- levels(argi.hz$depthcl)[2:5]
bdepths <- c(50, 100, 150, 200)
xposition <-  c(0.8, 0.8, 0.8, 0.8)
for(dc in 1:length(depth.classes)) {
  plot(y=argi.hz$hzn_bot, x=runif(nrow(argi.hz), 0, 1), ylim=c(200,0), xlim=c(0,1), type="n",
       xlab="Probability", ylab="Soil Depth, cm")
  pH_levels <- seq(5,7,0.5)
  colorz <- viridis::viridis(length(pH_levels))
  names(colorz) <- as.character(pH_levels)
  for(pH in pH_levels) {
    newdata <- data.frame(ph_h2o = pH, hzn_bot = 0:bdepths[dc], depthcl=depth.classes[dc])
    newdata$rankP <- predict(mylogit, newdata = newdata, type = "response")
    lines(y=newdata$hzn_bot, x=newdata$rankP, ylim=c(150,0), col=colorz[as.character(pH)], lwd=2)
  }
  legend("bottomleft", title="Horizon pH", legend = pH_levels, col=colorz, lwd=2, cex=0.75)
  abline(v=0.5, col="RED", lty=2, lwd=3)
  text(xposition[dc], 10, depth.classes[dc], cex=1.2, font=4)
}
mtext(paste("Probability of",my.analyte,"greater than",my.threshold,"%\nAll Argillic/Kandic in MLRA",my.mlra,"KSSL DB"), outer = TRUE, side = 1, line = 2)
dev.off()

my.pedons <- fetchNASIS()
site(my.pedons) <- getSoilDepthClass(my.pedons)
my.pedons <- my.pedons[grep(my.pedons$pedon_id, pattern="2020CA043|2019CA043"),]

my.pedons <- subsetProfiles(my.pedons, s='depth.class != "very.shallow"')
my.pedons$depth.class <- factor(my.pedons$depth.class)

odds <- as.numeric(exp(predict(mylogit, newdata=data.frame(ph_h2o=my.pedons$phfield, hzn_bot = my.pedons$hzdepb, depthcl=denormalize(my.pedons, 'depth.class')))))

# compute probability from odds ratio
my.pedons$is_typic_p <- odds / (1 + odds)

# "arbitrary" p threshold of 0.5 ~"more likely than not"
my.pedons$is_typic <- my.pedons$is_typic_p > 0.5

# default margins are too wide, make smaller
par(mar=c(1,1,4.1,1))
# subset just the prediction for recently collected CA043 pedons
plotSPC(my.pedons, color='is_typic_p', label="pedon_id", name = "",
        main="horizon is has BS 8.2 > 75%")        
addDiagnosticBracket(my.pedons, kind='argillic horizon',
                     tick.length = 0, lwd=10, col=rgb(red=0, green=0, blue=1, alpha=0.25))
addVolumeFraction(my.pedons, 'total_frags_pct')

my.pedons$hzdepthclass <- denormalize(my.pedons, 'depth.class')
lapply(split(horizons(my.pedons), f=my.pedons$hzdepthclass), function(d) {
  fm <- d$is_typic_p ~ d$phfield
  plot(fm, ylim=c(0,1), xlim=c(3,8))
  m <- lm(fm)
  abline(m)
  abline(v=6)
  abline(h=0.5)
  print(summary(m))
})

# now, summarise to make subgroup call based on subhorizons of argillic
subgrp.lut <- list("FALSE"='ultic',"TRUE"='typic')
taxsubgrp_est <- subgrp.lut[as.character(profileApply(my.pedons, function(p) {
  b <- getArgillicBounds(p)
  if(all(!is.na(b))) {
    sub.hz <- glom(p, b[1], b[2], as.data.frame = TRUE)
    # return true if all subhorizons have Typic probability >0.5
    return(all(na.omit(p$is_typic_p) > .75))
  } else {
    return("NA")
  }
}))]
names(taxsubgrp_est) <- my.pedons$pedon_id
table(unlist(taxsubgrp_est), useNA = "ifany")


