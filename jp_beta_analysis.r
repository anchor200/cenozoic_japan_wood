library(ClassDiscovery)
library(maps)
library(ggplot2)
library(reshape2)
library(betapart)
library(vegan)
library(geosphere)
library(dplyr)
library(stringr)
library(tidyr)
library(RColorBrewer)

setwd("C:/Users/ikari/Documents/laboratory/projects/202105_kaseki_paper/cenozoic_beta_plant_jp-main")

# reading geographic information
regions <- read.delim("data/region_mask.tsv", sep='\t', header=F, stringsAsFactor=FALSE)
colnames(regions) <- 120:149
rownames(regions) <- 45:20

seas <- read.delim("data/sep_sea.tsv", sep='\t', header=F, stringsAsFactor=FALSE)
colnames(seas) <- 120:149
rownames(seas) <- 45:20

islands <- read.delim("data/sep_island.tsv", sep='\t', header=F, stringsAsFactor=FALSE)
colnames(islands) <- 120:149
rownames(islands) <- 45:20

culdata <- data.frame(matrix(rep(NA, 5), nrow=1))[numeric(0), ]
colnames(culdata) <- c("lat", "lon", "region", "sea", "island")

for (i in 120:149){
    for (j in 20:45){
        temp <- regions[as.character(j), as.character(i)]
        temps <- seas[as.character(j), as.character(i)]
        tempi <- islands[as.character(j), as.character(i)]
        if (temp!=-1){
            tempRow <- data.frame("lat"=j, "lon"=i, "region"=temp, "sea"=temps, "island"=tempi)
            culdata <- rbind(culdata, tempRow)
        }
    }
}

# naming regions
culdata$Okhotsk <- 0
culdata$Okhotsk[culdata$sea==0] <- 1
culdata$Pacific1 <- 0
culdata$Pacific1[culdata$sea==2] <- 1
culdata$Pacific2 <- 0
culdata$Pacific2[culdata$sea==3] <- 1
culdata$Pacific3 <- 0
culdata$Pacific2[culdata$sea==6] <- 1
culdata$Japan <- 0
culdata$Japan[culdata$sea==4] <- 1
culdata$China <- 0
culdata$China[culdata$sea==5] <- 1

culdata$Hokkaido <- 0
culdata$Hokkaido[culdata$island==0] <- 1
culdata$Honshu <- 0
culdata$Honshu[culdata$island==1] <- 1
culdata$Shikoku <- 0
culdata$Shikoku[culdata$island==3] <- 1
culdata$Kushu <- 0
culdata$Kushu[culdata$island==2] <- 1
culdata$Okinawa <- 0
culdata$Okinawa[culdata$island==4] <- 1
culdata$Okinawa2 <- 0
culdata$Okinawa2[culdata$island==5] <- 1

culdata$key <- paste(culdata$lat, culdata$lon, sep="_")
culdata <- culdata[, colnames(culdata)!="sea" & colnames(culdata)!="island"]
culdata_resh <- culdata
culdata_resh$region_and_key <- paste(culdata_resh$key, culdata_resh$region, sep="_")
culdata_resh <- culdata_resh[order(culdata_resh$key),]


## analysis setting
mode <- c(3,14,17)
suppl <- ""
#mode <- c(2,10,12)
#suppl <- "_suppl"



# devision of Hokkaido
culdata_resh_Hokkaido <- subset(culdata_resh, culdata_resh$Hokkaido==1)
culdata_resh_kai <- culdata_resh_Hokkaido[, colnames(culdata_resh_Hokkaido)!="key" & colnames(culdata_resh_Hokkaido)!="region" & colnames(culdata_resh_Hokkaido)!="region_and_key"]
culdata_resh_mat <- data.frame(t(as.matrix(culdata_resh_kai)))
colnames(culdata_resh_mat) <- culdata_resh_Hokkaido$key
site_dist <- distanceMatrix(culdata_resh_mat, metric="euclidian")
site.cluster <- hclust(site_dist)
d <- as.dendrogram(site.cluster)
site.dend <- as.dendrogram(d)

## Figure A3(a)
plot(site.dend)


clusters<-cutree(site.cluster,k=mode[1])
cluster_resh <- data.frame(clusters)
rownames(cluster_resh) <- culdata_resh_Hokkaido$region_and_key
cluster_resh$lat <- 0
cluster_resh$lon <- 0
cluster_resh$clusters <- mode[1] - cluster_resh$clusters
for (i in 1:nrow(cluster_resh)){
    cluster_resh$lat[i] <- as.numeric(strsplit(rownames(cluster_resh), "_")[[i]][1])
    cluster_resh$lon[i] <- as.numeric(strsplit(rownames(cluster_resh), "_")[[i]][2])
}
cluster_resh_n <- cluster_resh

# division of Honshu, Shikoku, Kyushu
culdata_resh_Honshu <- subset(culdata_resh, culdata_resh$Hokkaido==0 & culdata_resh$Okinawa==0& culdata_resh$Okinawa2==0)
culdata_resh_Honshu <- subset(culdata_resh_Honshu, !(culdata_resh_Honshu$lat==41 & lon==139))
culdata_resh_kai <- culdata_resh_Honshu[, colnames(culdata_resh_Honshu)!="key" & colnames(culdata_resh_Honshu)!="region" & colnames(culdata_resh_Honshu)!="region_and_key"]
culdata_resh_mat <- data.frame(t(as.matrix(culdata_resh_kai)))
colnames(culdata_resh_mat) <- culdata_resh_Honshu$key

site_dist <- distanceMatrix(culdata_resh_mat, metric="euclidian")
site.cluster <- hclust(site_dist)
d <- as.dendrogram(site.cluster)
site.dend <- as.dendrogram(d)

## Figure A4(a)
plot(site.dend)


clusters<-cutree(site.cluster,k=mode[2])

cluster_resh <- data.frame(clusters)
rownames(cluster_resh) <- culdata_resh_Honshu$region_and_key
cluster_resh$clusters <- mode[2] - cluster_resh$clusters + mode[1]
#cluster_resh$clusters <- cluster_resh$clusters + 2
cluster_resh$lat <- 0
cluster_resh$lon <- 0
for (i in 1:nrow(cluster_resh)){
    cluster_resh$lat[i] <- as.numeric(strsplit(rownames(cluster_resh), "_")[[i]][1])
    cluster_resh$lon[i] <- as.numeric(strsplit(rownames(cluster_resh), "_")[[i]][2])
}
cluster_resh_c <- cluster_resh

# merge & write division table
site_division <- rbind(cluster_resh_n, cluster_resh_c)
write.table(site_division, paste("data/out/site_division", suppl, ".tsv", sep=""), sep="\t", quote=F, row.names=F)

world.map <- map_data ("world")
    japan <- world.map[world.map$long >= 120 & world.map$long < 150 & world.map$lat >= 20 & world.map$lat < 46 & world.map$region=="Japan",]


## Figure A2
g <- ggplot()+
geom_tile(aes(x=lon+0.5, y=lat+0.5, fill=factor(clusters), color=factor(clusters)), data=site_division) + ylim(30, 46)+
geom_path(data=japan, aes(x = long, y = lat, group = group), color="black")+ theme(text = element_text(size = 50))+
xlab("longitude") + ylab("latitude")+theme_bw()+ theme(text = element_text(size = 24))+xlim(127, 146)+
  guides(colour=guide_legend(override.aes = list(alpha=1,size=8)))+theme(legend.position = 'none') +
  geom_text(aes(x=lon+0.5, y=lat+0.5, label=factor(clusters)), size=7, data=site_division)
g





geo.ages <- c("Oligocene", "Miocene", "Pliocene", "Pleistocene", "LastGlacial", "Holocene", "Present")

# Tertiary genus occurrence
data_resh<-read.delim("data/Japan_Tertiary.csv", sep=',', header=T, stringsAsFactor=FALSE)
data_resh$Genus <- tolower(data_resh$Genus)
data_resh$lat_floor <- floor(data_resh$Latitude)
data_resh$long_floor <- floor(data_resh$Longitude)
data_resh$lat_long <- paste(data_resh$lat_floor, data_resh$long_floor, sep="_")

# Quartenary genus occurrence
data_qua<-read.delim("data/Japan_Quaternary.tsv", sep='\t', header=T, fileEncoding = "UTF-8-BOM")
data_qua$Genus <- tolower(data_qua$genus)
data_qua$Latitude <- data_qua$lat
data_qua$Longitude <- data_qua$lon
data_qua$lat_floor <- floor(data_qua$Latitude)
data_qua$long_floor <- floor(data_qua$Longitude)
data_qua$lat_long <- paste(data_qua$lat_floor, data_qua$long_floor, sep="_")


# Current genus occurrence
data_present_2d<-readRDS("data/present_tree_mesh2.rds")
data_present_2d$Latitude <- data_present_2d$lat
data_present_2d$Longitude <- data_present_2d$lon
data_present_2d$lat_floor <- floor(data_present_2d$Latitude)
data_present_2d$long_floor <- floor(data_present_2d$Longitude)
data_present_2d$Genus <- tolower(data_present_2d$genus)
data_present_2d$lat_long <- paste(data_present_2d$lat_floor, data_present_2d$long_floor, sep="_")


# Data description
print(c("Tertiary data points", nrow(data_resh)))
print(c("Tertiary sites", nrow(table(data_resh$Plot.name))))
print(c("Tertiary n genus", nrow(table(data_resh$Genus))))

print(c("Quartenary data points", nrow(data_qua)))
print(c("Quartenary sites", nrow(table(data_qua$chimei))))
print(c("Quartenary n genus", nrow(table(data_qua$Genus))))

table(data_qua$fossil_type)

data_resh$sitecode <- paste(data_resh$Latitude, data_resh$Longitude, sep="_")
data_resh_plots_Oligocene <- as.data.frame(table(data_resh$sitecode[data_resh$Geological.period=="Oligocene"]))
data_resh_plots_Miocene <- as.data.frame(table(data_resh$sitecode[data_resh$Geological.period=="Miocene"]))
data_resh_plots_Pliocene <- as.data.frame(table(data_resh$sitecode[data_resh$Geological.period=="Pliocene"]))

data_qua$sitecode <- paste(data_qua$lat, data_qua$lon, sep="_")
data_qua_plots_Pleistocene <- as.data.frame(table(data_qua$sitecode[data_qua$Geological.period=="Pleistocene"]))
data_qua_plots_LGP <- as.data.frame(table(data_qua$sitecode[data_qua$Geological.period=="Last glacial period"]))
data_qua_plots_Holocene <- as.data.frame(table(data_qua$sitecode[data_qua$Geological.period=="Holocene"]))

sites_num <- data.frame("period"=geo.ages[1:(length(geo.ages)-1)],
"sites"=c(nrow(data_resh_plots_Oligocene), nrow(data_resh_plots_Miocene), nrow(data_resh_plots_Pliocene),
nrow(data_qua_plots_Pleistocene), nrow(data_qua_plots_LGP), nrow(data_qua_plots_Holocene)),
"data"=c(length(data_resh$sitecode[data_resh$Geological.period=="Oligocene"]),
length(data_resh$sitecode[data_resh$Geological.period=="Miocene"]), 
length(data_resh$sitecode[data_resh$Geological.period=="Pliocene"]), 
length(data_qua$sitecode[data_qua$Geological.period=="Pleistocene"]),
length(data_qua$sitecode[data_qua$Geological.period=="Last glacial period"]),
length(data_qua$sitecode[data_qua$Geological.period=="Holocene"])))

sites_num









## options(repr.plot.width=4, repr.plot.height=4)
plot_co_places <- function(cene){

    data_cene <- subset(data_resh, data_resh$Geological.period==cene)    
    
    places <- rownames(table(data_cene$Plot.name))
    n_places <- nrow(table(data_cene$Plot.name))

    world.map <- map_data ("world")
    japan <- world.map[world.map$long >= 120 & world.map$long < 150 & world.map$lat >= 20 & world.map$lat < 46 & world.map$region=="Japan",]
    g <- ggplot(japan, aes(x = long, y = lat)) + geom_path(aes(group = group), color="blue")+
      layer(
        data=data_cene, 
        mapping=aes(x=Longitude, y=Latitude), 
        geom="point", 
        stat="identity", 
        position="identity"
      ) + ggtitle(paste(cene, "(",as.character(n_places), "sites)"))+theme_bw()+ theme(text = element_text(size = 24))
    return(g)
}

plot_co_places_qua <- function(cene){

    data_cene <- subset(data_qua, data_qua$Geological.period==cene)    
    
    places <- rownames(table(data_cene$Latitude))
    n_places <- nrow(table(data_cene$Latitude))
    world.map <- map_data ("world")
    japan <- world.map[world.map$long >= 120 & world.map$long < 150 & world.map$lat >= 20 & world.map$lat < 46 & world.map$region=="Japan",]
    g <- ggplot(japan, aes(x = long, y = lat)) + geom_path(aes(group = group), color="blue")+
      layer(
        data=data_cene, 
        mapping=aes(x=Longitude, y=Latitude), 
        geom="point", 
        stat="identity", 
        position="identity"
      ) + ggtitle(paste(cene, "(", as.character(n_places), "sites)"))+theme_bw()+ theme(text = element_text(size = 24))
    return(g)

}

#Tertiary data processing func
plot_co_region_matrix <- function(cene){

    data_cene <- subset(data_resh, data_resh$Geological.period==cene) 
    temp_aggr <- matrix(0, nrow=nGen, ncol=n_regions)
    rownames(temp_aggr) <- genus_list
    colnames(temp_aggr) <- regions

    for (reg in regions){
        place = subset(regions_resh, regions_resh$value==reg)$lat_long
        for (p in place){
            temp <- subset(data_cene, data_cene$lat_long==p)
            for (sp in temp$Genus){
                temp_aggr[sp, as.character(reg)] = temp_aggr[sp, as.character(reg)]  + 1
            }
        }
    }
    return(temp_aggr)

}

#Quartenary data processing func
plot_co_region_matrix_qua <- function(cene){

    data_cene <- subset(data_qua, data_qua$Geological.period==cene)    
    
    temp_aggr <- matrix(0, nrow=nGen, ncol=n_regions)
    rownames(temp_aggr) <- genus_list
    colnames(temp_aggr) <- regions

    for (reg in regions){
        place = subset(regions_resh, regions_resh$value==reg)$lat_long
        for (p in place){
            temp <- subset(data_cene, data_cene$lat_long==p)
            for (sp in temp$Genus){
                temp_aggr[sp, as.character(reg)] = temp_aggr[sp, as.character(reg)]  + 1
            }
        }
    }
    return(temp_aggr)

}


# extended default functions
merge2 <- function(dfs, ...)
{
 base <- dfs[1]
 lapply(dfs[-1], function(i) base <<- merge(base, i, ...)) # [1]
  return(base)
}
mean_ = function(x){
    return(mean(x, na.rm=TRUE))
}
sd_ = function(x){
    return(sd(x, na.rm=TRUE))
}

colname_correction_region <- function(cene){

    pa_by_region <- get(paste(cene, "_pa_by_region", sep=""))
    names <- colnames(pa_by_region)
    for (n in 1:length(names)){
        names[n] <- paste(cene, names[n], sep="")
    }
    colnames(pa_by_region) <- names
    temp <- as.data.frame(pa_by_region)
    temp$genus <- rownames(temp)
    return(temp)
}


# function for beta div calc
beta_order <- function(res, info=0){
    size <- ncol(as.matrix(res$beta.sim))
    n_pairs <- size * (size - 1) / 2
    output <- matrix(0, ncol=6, nrow=n_pairs)
    e_regions <- union(rownames(as.matrix(res$beta.sim)), colnames(as.matrix(res$beta.sim)))
    pairs <- combn(e_regions, m=2)        
    colnames(output) <- c("s_lat", "n_lat", "beta.sim", "beta.sne", "beta.sor", "info")

    for (i in 1:(n_pairs)){
        output[i,"s_lat"] <- gsub("\\D", "", pairs[,i][1])
        output[i,"n_lat"] <- gsub("\\D", "", pairs[,i][2])
        output[i,"beta.sim"] <- as.matrix(res$beta.sim)[pairs[,i][1],pairs[,i][2]]
        output[i,"beta.sne"] <- as.matrix(res$beta.sne)[pairs[,i][1],pairs[,i][2]]
        output[i,"beta.sor"] <- as.matrix(res$beta.sor)[pairs[,i][1],pairs[,i][2]]
        output[i,"info"] <- info
    }
    output <- as.data.frame(output)
    output$beta.sim <- as.numeric(output$beta.sim)
    output$beta.sne <- as.numeric(output$beta.sne)
    output$beta.sor <- as.numeric(output$beta.sor)
    return(output)
}

# function for regression
beta_show4 <- function(lat, ana=for_ana){
    res <- beta.pair(t(as.matrix(dplyr::select(ana, contains(lat)))))
    beta_order(res, lat)
}





beta.regr <- function(cene, data_aggr, method, clim=FALSE){

    da <- subset(data_aggr, data_aggr$info == cene)
    da$temp <- da$beta.sim
    
    if(clim==FALSE){

    
    if (method == "pow"){
        da$sp_dist <- log(da$sp_dist + .Machine$double.eps)
        res <- glm(1 - temp+.Machine$double.eps ~ sp_dist , family = gaussian(log), data=da)
    }else{
        res <- glm(1 - temp+.Machine$double.eps ~ sp_dist , family = gaussian(log), data=da)
    }
    res.temp <- as.data.frame(coef(summary(res)))
    res.temp$geo_peri <- cene
    
    return(list("coef"=res.temp, "pseudo_r2"=1 - res$deviance / res$null.deviance, "AIC"=AIC(res), "beta.mean"=mean_(da$beta.sim), "beta.sd"=sd_(da$beta.sim), "result"=res))

    }
    else{

    if (cene == "LastGlacial"){
        da$clim_dist <- da[,"clim_dist_LGP"]
    }
    else if(cene == "Present"){
        da$clim_dist <- da[,"clim_dist_Present"]
    }
    else{
        da$clim_dist <- da[,paste("clim_dist_", gsub("cene", "", cene), sep="")]        
    }


    if (method == "pow"){
        da$clim_dist <- log(da$clim_dist + .Machine$double.eps)
        res <- glm(1 - temp+.Machine$double.eps ~ clim_dist , family = gaussian(log), data=da)
    }else{
            res <- glm(1 - temp+.Machine$double.eps ~ clim_dist , family = gaussian(log), data=da)
    }
    res.temp <- as.data.frame(coef(summary(res)))
    res.temp$geo_peri <- cene
    
    return(list("coef"=res.temp, "pseudo_r2"=1 - res$deviance / res$null.deviance, "AIC"=AIC(res), "beta.mean"=mean_(da$beta.sim), "beta.sd"=sd_(da$beta.sim), "result"=res))
    }
}



# regression result reshape
beta.calc.by.ages <- function(method, xseq, ret, clim=FALSE){
    if (xseq == "clim"){
        return
    }
    else{
        dd <- plot_data_region_dist
    }
    
    if (ret == "AIC"){
        temp.ret <- rbind(
        beta.regr("Oligocene", dd, method, clim)$AIC,
        beta.regr("Miocene", dd, method, clim)$AIC,
        beta.regr("Pliocene", dd, method, clim)$AIC,
        beta.regr("Pleistocene", dd, method, clim)$AIC,
        beta.regr("LastGlacial", dd, method, clim)$AIC,
        beta.regr("Holocene", dd, method, clim)$AIC,
        beta.regr("Present", dd, method, clim)$AIC
     )
    }else if (ret == "pseudo_r2"){
        temp.ret <- rbind(
        beta.regr("Oligocene", dd, method, clim)$pseudo_r2,
        beta.regr("Miocene", dd, method, clim)$pseudo_r2,
        beta.regr("Pliocene", dd, method, clim)$pseudo_r2,
        beta.regr("Pleistocene", dd, method, clim)$pseudo_r2,
        beta.regr("LastGlacial", dd, method, clim)$pseudo_r2,
        beta.regr("Holocene", dd, method, clim)$pseudo_r2,
        beta.regr("Present", dd, method, clim)$pseudo_r2
         )
    }else{
        temp.ret <- rbind(
        beta.regr("Oligocene", dd, method, clim)$coef,
        beta.regr("Miocene", dd, method, clim)$coef,
        beta.regr("Pliocene", dd, method, clim)$coef,
        beta.regr("Pleistocene", dd, method, clim)$coef,
        beta.regr("LastGlacial", dd, method, clim)$coef,
        beta.regr("Holocene", dd, method, clim)$coef,
        beta.regr("Present", dd, method, clim)$coef
     )
    }

return(temp.ret)

}


mean__ <- function(x){
    return(list("mean"=mean_(x), "sd"=sd_(x)))
}
































# Genus list
genus_list <-union(union(union(rownames(table(data_present_2d$Genus)), rownames(table(data_qua$Genus))), rownames(table(data_resh$Genus))), rownames(table(data_present_2d$Genus)))
nGen <- length(genus_list)
genus_list_ter <-rownames(table(data_resh$Genus))
nGen_ter <- length(genus_list_ter)

# plot of fossile sites
options(repr.plot.width=24, repr.plot.height=16)
p1 <- plot_co_places("Oligocene")
p2 <- plot_co_places("Miocene")
p3 <- plot_co_places("Pliocene")
p4 <- plot_co_places_qua("Pleistocene")
p5 <- plot_co_places_qua("Last glacial period")
p6 <- plot_co_places_qua("Holocene")

gridExtra::grid.arrange( p1, p2, p3, p4, p5, p6, nrow = 2)

# read a sites data
sepfile <- paste("data/out/site_division_matrix", suppl, ".tsv", sep="")
regions <- read.delim(sepfile, sep='\t', header=F, stringsAsFactor=FALSE)

colnames(regions) <- 120:149
rownames(regions) <- 45:20
regions$lat <- rownames(regions)
regions_resh <- melt(regions, id.vars=c("lat"), variable.name="long")
regions_resh$lat_long <- paste(regions_resh$lat, regions_resh$long, sep="_")
regions_resh$lat <- as.numeric(regions_resh$lat)
regions_resh$long <- as.numeric(as.character(regions_resh$long))

## lat / lon values represent lower boundary of grid.


n_regions = 100 # just a big number for produce a loop
regions = 0:(n_regions-1)

regions_info <- data.frame(matrix(0, nrow=n_regions))
colnames <- c("region_id")
regions_info$region_id <- regions
regions_info <- regions_info[-1]
##  centering grid axes
for (id in regions){
    regions_info$lat_cent[regions_info$region_id==id] <- mean(subset(regions_resh, regions_resh$value==id)$lat) + 0.5
    regions_info$long_cent[regions_info$region_id==id] <- mean(subset(regions_resh, regions_resh$value==id)$long) + 0.5
}

rownames(regions_info) <- regions_info$region_id
regions_info <- regions_info[-1]

ggplot(data=regions_info, aes(x=long_cent, y=lat_cent)) + geom_text(aes(label=rownames(regions_info)))



Oligocene_pa_by_region <- plot_co_region_matrix("Oligocene")
Miocene_pa_by_region <- plot_co_region_matrix("Miocene")
Pliocene_pa_by_region <- plot_co_region_matrix("Pliocene")

Pleistocene_pa_by_region <- plot_co_region_matrix_qua("Pleistocene")
LastGlacial_pa_by_region <- plot_co_region_matrix_qua("Last glacial period")
Holocene_pa_by_region <- plot_co_region_matrix_qua("Holocene")

#Current data processing
    data_cene <-   data_present_2d
    temp_aggr <- matrix(0, nrow=nGen, ncol=n_regions)
    rownames(temp_aggr) <- genus_list
    colnames(temp_aggr) <- regions
    for (reg in regions){
        place = subset(regions_resh, regions_resh$value==reg)$lat_long
        for (p in place){
            temp <- subset(data_cene, data_cene$lat_long==p)
            for (sp in temp$Genus){
                temp_aggr[sp, as.character(reg)] = temp_aggr[sp, as.character(reg)]  +1
            }
        }
    }
Present_pa_by_region <- temp_aggr




whole_samples_region <- merge2(list(
    colname_correction_region("Oligocene"),
    colname_correction_region("Miocene"),
    colname_correction_region("Pliocene"),
    colname_correction_region("Pleistocene"),
    colname_correction_region("LastGlacial"),
    colname_correction_region("Holocene"),
    colname_correction_region("Present")
),by=c("genus"))

## data reshape for analyses
rownames(whole_samples_region) <- whole_samples_region$genus
for_ana_region <- whole_samples_region[-1]

## sample completeness calculation
## Table A1
pool_1 <- t(for_ana_region)
pool_1 <- pool_1[apply(pool_1, 1, sum)>1,]
pool_ <- gsub("\\d", "", rownames(pool_1))
pool <- specpool(pool_1, pool_)
pool$coverage <- pool$Species / pool$chao
pool

ppp <- pool_1[1:4,]
write.table(pool_1, "data/out/occurrence_freq.tsv", sep="\t", row.names=T)

## filtering with n of genus par subregion
CUT <- 10 # for results of Table A3, change this number.
for_ana_region[for_ana_region > 1] <- 1

for_ana_region <- for_ana_region[apply(for_ana_region, 2, sum)>CUT]

# beta diversity calculation
plot_data_region_dist <- rbind(
beta_show4("Oligocene", for_ana_region),
beta_show4("Miocene", for_ana_region),
beta_show4("Pliocene", for_ana_region),
beta_show4("Pleistocene", for_ana_region),
beta_show4("LastGlacial", for_ana_region),
beta_show4("Holocene", for_ana_region),
beta_show4("Present", for_ana_region)
) 

plot_data_region_dist$info <- factor(plot_data_region_dist$info, levels=geo.ages)
colnames(plot_data_region_dist) <- c("regionA", "regionB", "beta.sim", "beta.sne",  "beta.sor", "info")
plot_data_region_dist$regions <- paste(plot_data_region_dist$regionA, plot_data_region_dist$regionB, sep="_")
plot_data_region_dist$regionA_cent_lat<- regions_info[plot_data_region_dist$regionA,1]
plot_data_region_dist$regionA_cent_long<- regions_info[plot_data_region_dist$regionA,2]
plot_data_region_dist$regionB_cent_lat<- regions_info[plot_data_region_dist$regionB,1]
plot_data_region_dist$regionB_cent_long<- regions_info[plot_data_region_dist$regionB,2]
plot_data_region_dist$dist <- 0

for (k in 1:nrow(plot_data_region_dist)){
    plot_data_region_dist$dist[k] <- distGeo(c(plot_data_region_dist$regionA_cent_long[k], plot_data_region_dist$regionA_cent_lat[k]), 
                                             c(plot_data_region_dist$regionB_cent_long[k], plot_data_region_dist$regionB_cent_lat[k])) / 1000
}

plot_data_region_dist$label <- paste(plot_data_region_dist$regionA, plot_data_region_dist$regionB, sep="_")
plot_data_region_dist$info <- factor(plot_data_region_dist$info, levels=geo.ages)
plot_data_region_dist$sp_dist <- plot_data_region_dist$dist



table_resher <- function(coefs, fits, clim=F){

	if(!clim){
		rowss <- c("Intercept", "std.error", "t-value", "p-value", "Spatial distance (x1000)", "std.error", "t-value", "p-value", "Pseudo-R2", "AIC")
	}else{
		rowss <- c("Intercept", "std.error", "t-value", "p-value", "Climatic distance", "std.error", "t-value", "p-value", "Pseudo-R2", "AIC")
	}

	temp_table <- cbind(coefs[coefs$predictor=="(Intercept)",c(1,2,3,4)],
	coefs[coefs$predictor!="(Intercept)",c(1,2,3,4)])
	temp_table <- data.frame(t(as.matrix(temp_table)))
	colnames(temp_table) <- geo.ages
	neg_exp_fit2 <- data.frame(t(as.matrix(fits)))
	neg_exp_fit2 <- neg_exp_fit2[c(2,3),]
	colnames(neg_exp_fit2) <- geo.ages
	temp <- rbind(temp_table, neg_exp_fit2)
	temp$rn <- rowss
	temp <- temp[, c(ncol(temp), 1:(ncol(temp)-1))]
	return(temp)
}



# spatial distance-decay
## negative exp
## Table A4(a) # A2(a)
neg_exp_coefs_distonly_region <- beta.calc.by.ages("neg", "space", "ret")
neg_exp_coefs_distonly_region$predictor <- gsub("\\d", "", rownames(neg_exp_coefs_distonly_region))
colnames(neg_exp_coefs_distonly_region) <- gsub("\\s", "",colnames(neg_exp_coefs_distonly_region))
#neg_exp_coefs_distonly_region

# fit measures
## Table A4(a) # A2(a)
neg_exp_aic <- beta.calc.by.ages("neg", "space", "AIC", clim=FALSE)
neg_exp_r2 <- beta.calc.by.ages("neg", "space", "pseudo_r2", clim=FALSE)
neg_exp_fit <- data.frame("info"=geo.ages, "pseudo_r2"=neg_exp_r2, "AIC"=neg_exp_aic)
#neg_exp_fit

a4a <- table_resher(neg_exp_coefs_distonly_region, neg_exp_fit)
write.table(a4a, "data/out/table_a4a.tsv", sep="\t", row.names=F, quote=F)


## Figure 2(b)
neg_exp_coef_dist <- subset(neg_exp_coefs_distonly_region, neg_exp_coefs_distonly_region$predictor!="(Intercept)")
neg_exp_coef_dist$geological_age <- factor(neg_exp_coef_dist$geo_peri, levels=geo.ages)

pdf("data/out/fig_2b.pdf", width = 9, height = 6)
gr <- barplot(-neg_exp_coef_dist$Estimate*1000 ~ neg_exp_coef_dist$geological_age, ylim=c(-0.4, 0.7), main="spacial distance-tunover relationship", xlab="geological age", ylab="turnover rate (/1000km)")
arrows(gr, -neg_exp_coef_dist$Estimate*1000+neg_exp_coef_dist$Std.Error*1.96*1000, gr, -neg_exp_coef_dist$Estimate*1000-neg_exp_coef_dist$Std.Error*1.96*1000, code = 3, angle = 90, lwd = 1, length = 0.1, data=neg_exp_coef_dist)
text(x=7.9,y=0.6,"***", cex=2)
text(x=5.5,y=0.6,"*", cex=2)
text(x=6.7,y=0.6,"*", cex=2)
abline(h=0)
dev.off()

## power law
## Table A4(b) # A2(b)
pow_coefs_distonly_region <- beta.calc.by.ages("pow", "space", "ret")
pow_coefs_distonly_region$predictor <- gsub("\\d", "", rownames(pow_coefs_distonly_region))
colnames(pow_coefs_distonly_region) <- gsub("\\s", "",colnames(pow_coefs_distonly_region))
pow_coefs_distonly_region

pow_exp_aic <- beta.calc.by.ages("pow", "space", "AIC", clim=FALSE)
pow_exp_r2 <- beta.calc.by.ages("pow", "space", "pseudo_r2", clim=FALSE)
pow_exp_fit <- data.frame("info"=geo.ages, "pseudo_r2"=pow_exp_r2, "AIC"=pow_exp_aic)

a4b <- table_resher(pow_coefs_distonly_region, pow_exp_fit)
write.table(a4b, "data/out/table_a4b.tsv", sep="\t", row.names=F, quote=F)




# climatic distance-decay
paleoclim<-read.delim("data/paleoclim_degree.csv", sep=',', header=T, stringsAsFactor=FALSE)
paleoclim$int_gridkey <- paste(floor(paleoclim$y), floor(paleoclim$x), sep="_")
paleoclim_resh <- paleoclim[c("int_gridkey", "Oligo", "Mio", "Plio", "Pleist", "LGP", "Holo", "current")]
temp_merged <- merge(regions_resh,paleoclim_resh,by.x="lat_long", by.y="int_gridkey", all=T)
clim_region <- subset(temp_merged, temp_merged$value>=0)

clim_resh <- aggregate(list(
               "Oligo"=clim_region$Oligo,
               "Mio"=clim_region$Mio,
               "Plio"=clim_region$Plio,
               "Pleist"=clim_region$Pleist,
               "LGP"=clim_region$LGP,
               "Holo"=clim_region$Holo,
               "Present"=clim_region$current),by=list("region"=clim_region$value), mean__)

# compile climate data
for (k in 1:nrow(plot_data_region_dist)){
    plot_data_region_dist$clim_dist_Oligo[k] <- abs(clim_resh$Oligo[[as.numeric(plot_data_region_dist$regionA[k])+1,1]]-
    clim_resh$Oligo[[as.numeric(plot_data_region_dist$regionB[k])+1,1]])
    plot_data_region_dist$clim_dist_Mio[k] <- abs(clim_resh$Mio[[as.numeric(plot_data_region_dist$regionA[k])+1,1]] -
    clim_resh$Mio[[as.numeric(plot_data_region_dist$regionB[k])+1,1]])
    plot_data_region_dist$clim_dist_Plio[k] <- abs(clim_resh$Plio[[as.numeric(plot_data_region_dist$regionA[k])+1,1]]-
    clim_resh$Plio[[as.numeric(plot_data_region_dist$regionB[k])+1,1]])
    plot_data_region_dist$clim_dist_Pleisto[k] <- abs(clim_resh$Pleist[[as.numeric(plot_data_region_dist$regionA[k])+1,1]] -
    clim_resh$Pleist[[as.numeric(plot_data_region_dist$regionB[k])+1,1]])
    plot_data_region_dist$clim_dist_LGP[k] <- abs(clim_resh$LGP[[as.numeric(plot_data_region_dist$regionA[k])+1,1]] -
    clim_resh$LGP[[as.numeric(plot_data_region_dist$regionB[k])+1,1]])
    plot_data_region_dist$clim_dist_Holo[k] <- abs(clim_resh$Holo[[as.numeric(plot_data_region_dist$regionA[k])+1,1]] -
    clim_resh$Holo[[as.numeric(plot_data_region_dist$regionB[k])+1,1]])
    plot_data_region_dist$clim_dist_Present[k] <- abs(clim_resh$Present[[as.numeric(plot_data_region_dist$regionA[k])+1,1]] -
    clim_resh$Present[[as.numeric(plot_data_region_dist$regionB[k])+1,1]])
}

## negative exp
## Table A4(c) #A2(c)
neg_exp_coefs_clim_region <- beta.calc.by.ages("neg", "space", "ret", clim=TRUE)
neg_exp_coefs_clim_region$predictor <- gsub("\\d", "", rownames(neg_exp_coefs_clim_region))
colnames(neg_exp_coefs_clim_region) <- gsub("\\s", "",colnames(neg_exp_coefs_clim_region))
neg_exp_coefs_clim_region

neg_exp_coef_dist_c <- subset(neg_exp_coefs_clim_region, neg_exp_coefs_clim_region$predictor!="(Intercept)")
neg_exp_coef_dist_c$geological_age <- factor(neg_exp_coef_dist_c$geo_peri, levels=geo.ages)

# fit measures
## Table A4(c, d) #A2(c, d)
neg_exp_aic_clim <- beta.calc.by.ages("neg", "space", "AIC", clim=TRUE)
neg_exp_r2_clim <- beta.calc.by.ages("neg", "space", "pseudo_r2", clim=TRUE)
neg_exp_fit_clim <- data.frame("info"=c("Oligocene", "Miocene", "Pliocene", "Pleistocene", "LastGlacial", "Holocene", "Present"), "pseudo_r2"=neg_exp_r2_clim, "AIC"=neg_exp_aic_clim)
#neg_exp_fit_clim

a4c <- table_resher(neg_exp_coefs_clim_region, neg_exp_fit_clim, TRUE)
write.table(a4c, "data/out/table_a4c.tsv", sep="\t", row.names=F, quote=F)



## Figure 3(b)
pdf("data/out/fig_3b.pdf", width = 9, height = 6)
gr <- barplot(-Estimate ~ geological_age, ylim=c(-0.06, 0.07), data=neg_exp_coef_dist_c, main="climatic distance-tunover relationship", xlab="geological age", ylab="turnover rate (/degree)", las=1)
arrows(gr, -neg_exp_coef_dist_c$Estimate+neg_exp_coef_dist_c$Std.Error*1.96, gr, -neg_exp_coef_dist_c$Estimate-neg_exp_coef_dist_c$Std.Error*1.96, code = 3, angle = 90, lwd = 1, length = 0.1, data=neg_exp_coef_dist_c)
text(x=7.9,y=0.06,"***", cex=2)
text(x=5.5,y=0.06,"*", cex=2)
text(x=6.7,y=0.06,"*", cex=2)
abline(h=0)
dev.off()

## power law
## Table A4(d) #A2(d)
pow_coefs_clim_region <- beta.calc.by.ages("pow", "space", "ret", clim=TRUE)
pow_coefs_clim_region$predictor <- gsub("\\d", "", rownames(pow_coefs_clim_region))
colnames(pow_coefs_clim_region) <- gsub("\\s", "",colnames(pow_coefs_clim_region))

pow_aic_clim <- beta.calc.by.ages("pow", "space", "AIC", clim=TRUE)
pow_r2_clim <- beta.calc.by.ages("pow", "space", "pseudo_r2", clim=TRUE)
pow_fit_clim <- data.frame("info"=c("Oligocene", "Miocene", "Pliocene", "Pleistocene", "LastGlacial", "Holocene", "Present"), "pseudo_r2"=pow_r2_clim, "AIC"= pow_aic_clim)

a4d <- table_resher(pow_coefs_clim_region, pow_fit_clim, TRUE)
write.table(a4d, "data/out/table_a4d.tsv", sep="\t", row.names=F, quote=F)




cene_decay_plot <- function(cene, col_, clim=FALSE){

    if (clim == TRUE){
        dist_seq <- seq(from=0, to=30, by=0.2)
        plot_decay <- data.frame("dist"=dist_seq)
        newx <- data.frame("clim_dist" = seq(0, 30, 0.2)) # 下で描画に使う
        icpt <- neg_exp_coefs_clim_region$Estimate[neg_exp_coefs_distonly_region$geo_peri==cene & neg_exp_coefs_clim_region$predictor=="(Intercept)"]
        slop <- neg_exp_coefs_clim_region$Estimate[neg_exp_coefs_distonly_region$geo_peri==cene & neg_exp_coefs_clim_region$predictor=="clim_dist"]
        plot_decay$pred <- 1-exp((slop * plot_decay$dist + icpt))

        plotdata <- subset(clim_dist_plot_clim_region, clim_dist_plot_clim_region$info==cene) 
        plotdata$dist <- plotdata$clim_dist
        
        plot(beta.sim ~ dist, plotdata, main=paste(cene, ""), xlab="delta temperature [degree]", ylab="turnover", col=col_, pch=20,
        cex.lab  = 2, 
        cex.axis = 1.8, 
        cex.main = 1.8,
        ylim=c(0.0, 0.65), las=1, xlim=c(0, 15))
        
    }else{
        dist_seq <- seq(from=0, to=3000, by=10)
        plot_decay <- data.frame("dist"=dist_seq)
        newx <- data.frame("sp_dist" = seq(0, 3000, 20))
        icpt <- neg_exp_coefs_distonly_region$Estimate[neg_exp_coefs_distonly_region$geo_peri==cene & neg_exp_coefs_distonly_region$predictor=="(Intercept)"]
        slop <- neg_exp_coefs_distonly_region$Estimate[neg_exp_coefs_distonly_region$geo_peri==cene & neg_exp_coefs_distonly_region$predictor=="sp_dist"]
        plot_decay$pred <- 1-exp((slop * plot_decay$dist + icpt))
        
        plotdata <- subset(plot_data_region_dist, plot_data_region_dist$info==cene) 
        

        plot(beta.sim ~ dist, plotdata, main=paste(cene, ""), xlab="spatial distance [km]", ylab="turnover", col=col_, pch=20,
        cex.lab  = 2,
        cex.axis = 1.8,
        cex.main = 1.8,
        ylim=c(0.0, 0.65), las=1, xlim=c(0, 2000))
        
    }

    if (clim == TRUE){
            m <- beta.regr(cene, plot_data_region_dist, "neg", clim=TRUE)$result   
    }else{
            m <- beta.regr(cene, plot_data_region_dist, "neg")$result
    }

    preds <- predict(m, newdata = newx, se.fit = TRUE, type = "link")
    critval <- 1.96 ## approx 95% CI
    upr <- 1 - exp(preds$fit + (critval * preds$se.fit))
    lwr <- 1 - exp(preds$fit - (critval * preds$se.fit))
    fit <- 1 - exp(preds$fit)
    if (clim == TRUE){
        lines(newx$clim_dist, upr, col = 'darkgreen')
        lines(newx$clim_dist, fit, col = 'orange')
        lines(newx$clim_dist, lwr, col = 'darkgreen')
    }else{
        lines(newx$sp_dist, upr, col = 'darkgreen')
        lines(newx$sp_dist, fit, col = 'orange')
        lines(newx$sp_dist, lwr, col = 'darkgreen')
    }
            
}

pdf("data/out/fig_2a.pdf", width = 12, height = 9)
options(repr.plot.width=16, repr.plot.height=12)
## Figure 2(a)
par(mfrow=c(3,3)) 
for (i in 1:length(geo.ages)){
cene_decay_plot(geo.ages[i], "grey", FALSE)
}
dev.off()


clim_dist_add <- function(cene, data_aggr){
    da <- subset(data_aggr, data_aggr$info == cene)
    if (cene == "LastGlacial"){
        da$clim_dist <- da[,"clim_dist_LGP"]
    }
    else if(cene == "Present"){
        da$clim_dist <- da[,"clim_dist_Present"]
    }
    else{
        da$clim_dist <- da[,paste("clim_dist_", gsub("cene", "", cene), sep="")]        
    }
        return(da)
}

clim_dist_plot_clim_region <- rbind(
clim_dist_add("Oligocene", plot_data_region_dist),
clim_dist_add("Miocene", plot_data_region_dist),
clim_dist_add("Pliocene", plot_data_region_dist),
clim_dist_add("Pleistocene", plot_data_region_dist),
clim_dist_add("LastGlacial", plot_data_region_dist),
clim_dist_add("Holocene", plot_data_region_dist),
clim_dist_add("Present", plot_data_region_dist)
)



pdf("data/out/fig_3a.pdf", width = 12, height = 9)
## Figure 3(a)
par(mfrow=c(3,3)) 
for (i in 1:length(geo.ages)){
cene_decay_plot(geo.ages[i], "grey", TRUE)
}
dev.off()


# Clim etc plot
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
clim_resh$lat_cent <- regions_info$lat_cent[1:nrow(clim_resh)]
clim_resh$long_cent <- regions_info$long_cent[1:nrow(clim_resh)]
#head(clim_resh)

clim_long <- clim_resh %>%
  pivot_longer(-c(region, lat_cent, long_cent), names_to = "age", values_to = "temperature")
head(clim_long)

clim_long$age <- factor(clim_long$age, c("Oligo", "Mio", "Plio", "Pleist", "LGP", "Holo", "Present"))

clim_long$temperature.mean <- as.numeric(clim_long$temperature[,1])
clim_long$temperature.sd <- as.numeric(clim_long$temperature[,2])
head(clim_long)

world.map <- map_data ("world")
    japan <- world.map[world.map$long >= 120 & world.map$long < 150 & world.map$lat >= 20 & world.map$lat < 46 & world.map$region=="Japan",]

regions_resh_plot <- subset(regions_resh, regions_resh$value >= 0)
for (i in regions_resh_plot$value){
regions_resh_plot$Oligocene[regions_resh_plot$value==i] <- clim_long$temperature.mean[clim_long$region==i & clim_long$age=="Oligo"]    
}
for (i in regions_resh_plot$value){
regions_resh_plot$Miocene[regions_resh_plot$value==i] <- clim_long$temperature.mean[clim_long$region==i & clim_long$age=="Mio"]    
}
for (i in regions_resh_plot$value){
regions_resh_plot$Pliocene[regions_resh_plot$value==i] <- clim_long$temperature.mean[clim_long$region==i & clim_long$age=="Plio"]    
}
for (i in regions_resh_plot$value){
regions_resh_plot$Pleistocene[regions_resh_plot$value==i] <- clim_long$temperature.mean[clim_long$region==i & clim_long$age=="Pleist"]    
}
for (i in regions_resh_plot$value){
regions_resh_plot$LastGlacial[regions_resh_plot$value==i] <- clim_long$temperature.mean[clim_long$region==i & clim_long$age=="LGP"]    
}
for (i in regions_resh_plot$value){
regions_resh_plot$Holocene[regions_resh_plot$value==i] <- clim_long$temperature.mean[clim_long$region==i & clim_long$age=="Holo"]    
}
for (i in regions_resh_plot$value){
regions_resh_plot$Present[regions_resh_plot$value==i] <- clim_long$temperature.mean[clim_long$region==i & clim_long$age=="Present"]    
}

head(regions_resh_plot)

regions_resh_plot_long <- regions_resh_plot %>%
  pivot_longer(-c(lat, long, value, lat_long), names_to = "age", values_to = "temperature")

regions_resh_plot_long$age <- factor(regions_resh_plot_long$age, levels=c('Oligocene', 'Miocene', 'Pliocene', 'Pleistocene', 'LastGlacial', 'Holocene', 'Present'))

regions_resh_plot_long$istropics <- as.numeric(regions_resh_plot_long$temperature>=18)

kaseki_sites <- rbind(data.frame("age"=data_resh$Geological.period, "lat"=data_resh$Latitude, "long"=data_resh$Longitude),
data.frame("age"=data_qua$Geological.period, "lat"=data_qua$Latitude, "long"=data_qua$Longitude))
kaseki_sites <- subset(kaseki_sites, age!="Eocene")
kaseki_sites$age[kaseki_sites$age=="Last glacial period"] <- "LastGlacial"
kaseki_sites$age <- factor(kaseki_sites$age, levels=c('Oligocene', 'Miocene', 'Pliocene', 'Pleistocene', 'LastGlacial', 'Holocene', 'Present'))


g <- ggplot()+
geom_tile(aes(x=long+0.5, y=lat+0.5, fill=temperature, color=temperature), data=regions_resh_plot_long) +
scale_fill_gradientn(colours = myPalette(100), limits=c(0, 35)) + 
scale_colour_gradientn(colours = myPalette(100), limits=c(0, 35)) + facet_wrap(~age) +
geom_path(data=japan, aes(x = long, y = lat, group = group), color="black")+theme(legend.title=element_text(size = 18), strip.text.x = element_text(size = 20), axis.text.x=element_text(size=12), axis.text.y=element_text(size=12), axis.title.x=element_text(size=15), axis.title.y=element_text(size=15))+
xlab("longitude") + ylab("latitude")+theme_bw() + xlim(127, 146) + ylim(29, 46)+ theme(text = element_text(size = 24))+
geom_tile(aes(x=long+0.5, y=lat+0.5),alpha=0.2, data=subset(regions_resh_plot_long, istropics==1))+
geom_point(aes(x=long, y=lat), size=1, data=kaseki_sites, color="blue")
ggsave(plot=g, file = "data/out/fig_a1a.png", width=10, height=10, dpi=600)


# Quaternary LDG

lats <- 30:45

data_ <- data_qua[data_qua$Geological.period=="Pleistocene", ]
data_qua_ldg <- NA

for(i in lats){
	temp <- data_[data_$lat_floor==i, ]
	temp_n <- length(unique(temp$genus))
	td <- data.frame("lat"=i, "n_sp"=temp_n, "Period"="d. Pleistocene")
	if(is.na(data_qua_ldg)){
		data_qua_ldg <- td
	}else{
		data_qua_ldg <- rbind(data_qua_ldg, td)
	}
}


data_ <- data_qua[data_qua$Geological.period=="Holocene", ]
for(i in lats){
	temp <- data_[data_$lat_floor==i, ]
	temp_n <- length(unique(temp$genus))
	td <- data.frame("lat"=i, "n_sp"=temp_n, "Period"="f. Holocene")
	if(is.na(data_qua_ldg)){
		data_qua_ldg <- td
	}else{
		data_qua_ldg <- rbind(data_qua_ldg, td)
	}
}

data_ <- data_qua[data_qua$Geological.period=="Last glacial period", ]
for(i in lats){
	temp <- data_[data_$lat_floor==i, ]
	temp_n <- length(unique(temp$genus))
	td <- data.frame("lat"=i, "n_sp"=temp_n, "Period"="e. Last glacial period")
	if(is.na(data_qua_ldg)){
		data_qua_ldg <- td
	}else{
		data_qua_ldg <- rbind(data_qua_ldg, td)
	}
}


data_ <- data_present_2d
for(i in lats){
	temp <- data_[data_$lat_floor==i, ]
	temp_n <- length(unique(temp$genus))
	td <- data.frame("lat"=i, "n_sp"=temp_n, "Period"="g. Present")
	if(is.na(data_qua_ldg)){
		data_qua_ldg <- td
	}else{
		data_qua_ldg <- rbind(data_qua_ldg, td)
	}
}


data_ter_ldg <- NA


data_ <- data_resh[data_resh$Geological.period=="Oligocene", ]

for(i in lats){
	temp <- data_[data_$lat_floor==i, ]
	temp_n <- length(unique(temp$Genus))
	td <- data.frame("lat"=i, "n_sp"=temp_n, "Period"="a. Oligocene")
	if(is.na(data_ter_ldg)){
		data_ter_ldg <- td
	}else{
		data_ter_ldg <- rbind(data_ter_ldg, td)
	}
}

data_ <- data_resh[data_resh$Geological.period=="Miocene", ]

for(i in lats){
	temp <- data_[data_$lat_floor==i, ]
	temp_n <- length(unique(temp$Genus))
	td <- data.frame("lat"=i, "n_sp"=temp_n, "Period"="b. Miocene")
	if(is.na(data_ter_ldg)){
		data_ter_ldg <- td
	}else{
		data_ter_ldg <- rbind(data_ter_ldg, td)
	}
}

data_ <- data_resh[data_resh$Geological.period=="Pliocene", ]

for(i in lats){
	temp <- data_[data_$lat_floor==i, ]
	temp_n <- length(unique(temp$Genus))
	td <- data.frame("lat"=i, "n_sp"=temp_n, "Period"="c. Pliocene")
	if(is.na(data_ter_ldg)){
		data_ter_ldg <- td
	}else{
		data_ter_ldg <- rbind(data_ter_ldg, td)
	}
}



evidences <- data.frame("evidence"=c(
summary(lm(n_sp ~ lat, data=data_ter_ldg[data_ter_ldg$n_sp != 0 & data_ter_ldg$"Period"=="a. Oligocene",]))$adj.r.squared,#$coefficients[2,4],
summary(lm(n_sp ~ lat, data=data_ter_ldg[data_ter_ldg$n_sp != 0 & data_ter_ldg$"Period"=="b. Miocene",]))$adj.r.squared,#$coefficients[2,4],
summary(lm(n_sp ~ lat, data=data_ter_ldg[data_ter_ldg$n_sp != 0 & data_ter_ldg$"Period"=="c. Pliocene",]))$adj.r.squared,#$coefficients[2,4],
summary(lm(n_sp ~ lat, data=data_qua_ldg[data_qua_ldg$n_sp != 0 & data_qua_ldg$"Period"=="d. Pleistocene",]))$adj.r.squared,#$coefficients[2,4],
summary(lm(n_sp ~ lat, data=data_qua_ldg[data_qua_ldg$n_sp != 0 & data_qua_ldg$"Period"=="e. Last glacial period",]))$adj.r.squared,#$coefficients[2,4],
summary(lm(n_sp ~ lat, data=data_qua_ldg[data_qua_ldg$n_sp != 0 & data_qua_ldg$"Period"=="f. Holocene",]))$adj.r.squared,#$coefficients[2,4],
summary(lm(n_sp ~ lat, data=data_qua_ldg[data_qua_ldg$n_sp != 0 & data_qua_ldg$"Period"=="g. Present",]))$adj.r.squared#$coefficients[2,4]
),
"Period"=factor(geo.ages, levels=geo.ages))

gsl1 <- evidences %>%
ggplot(aes(x=factor(geo.ages, levels=geo.ages), y=evidence)) + geom_bar(stat="identity") + theme_bw()+
xlab("geological period") + ylab("evidence (linear regression adjusted R2) of LDG")+ theme(text = element_text(size = 24)) 
ggsave(plot=gsl1, file = "data/out/fig_gsl1.pdf", width=13, height=10)

evidences <- data.frame("evidence"=c(
summary(lm(n_sp ~ lat, data=data_ter_ldg[data_ter_ldg$n_sp != 0 & data_ter_ldg$"Period"=="a. Oligocene",]))$coefficients[2,4],
summary(lm(n_sp ~ lat, data=data_ter_ldg[data_ter_ldg$n_sp != 0 & data_ter_ldg$"Period"=="b. Miocene",]))$coefficients[2,4],
summary(lm(n_sp ~ lat, data=data_ter_ldg[data_ter_ldg$n_sp != 0 & data_ter_ldg$"Period"=="c. Pliocene",]))$coefficients[2,4],
summary(lm(n_sp ~ lat, data=data_qua_ldg[data_qua_ldg$n_sp != 0 & data_qua_ldg$"Period"=="d. Pleistocene",]))$coefficients[2,4],
summary(lm(n_sp ~ lat, data=data_qua_ldg[data_qua_ldg$n_sp != 0 & data_qua_ldg$"Period"=="e. Last glacial period",]))$coefficients[2,4],
summary(lm(n_sp ~ lat, data=data_qua_ldg[data_qua_ldg$n_sp != 0 & data_qua_ldg$"Period"=="f. Holocene",]))$coefficients[2,4],
summary(lm(n_sp ~ lat, data=data_qua_ldg[data_qua_ldg$n_sp != 0 & data_qua_ldg$"Period"=="g. Present",]))$coefficients[2,4]
),
"Period"=factor(geo.ages, levels=geo.ages))

gsl2 <- evidences %>%
ggplot(aes(x=factor(geo.ages, levels=geo.ages), y=evidence)) + geom_bar(stat="identity") + theme_bw()+
xlab("geological period") + ylab("evidence (linear regression p-value) of LDG")+ theme(text = element_text(size = 24)) 
ggsave(plot=gsl2, file = "data/out/fig_gsl2.pdf", width=13, height=10)



temp <- rbind(data_ter_ldg, data_qua_ldg)

gsl3 <- temp[temp$n_sp != 0,] %>%
filter(Period != "g. Present") %>%
ggplot(aes(x=lat+0.5, y=n_sp, color=Period)) + geom_point() + stat_smooth(se=F,method = "loess",method.args= list(degree = 1), span = 2)+ theme_bw()+ theme(text = element_text(size = 24))+
xlab("latitude") + ylab("number of genera")
ggsave(plot=gsl3, file = "data/out/fig_gsl3.pdf", width=10, height=10)

gsl4 <- temp[temp$n_sp != 0,] %>%
filter(Period == "g. Present") %>%
ggplot(aes(x=lat+0.5, y=n_sp, color=Period)) + geom_point() + stat_smooth(se=F,method = "loess",method.args= list(degree = 1), span = 2)+ theme_bw()+ theme(text = element_text(size = 24))+
xlab("latitude") + ylab("number of genera")
ggsave(plot=gsl4, file = "data/out/fig_gsl4.pdf", width=10, height=10)

