# Beta diversity calculation (for robust analyses)
pairwise_beta <- function(lat, ana=for_ana){
    ## for error skip during loop
    if (nrow(as.matrix(dplyr::select(ana, contains(lat)))) >= 3  ){
        res <- try(beta.pair(t(as.matrix(dplyr::select(ana, contains(lat))))), silent = FALSE) 
        if (class(res) == "try-error") {
                return(data.frame(
                    "s_lat"="none",
                    "n_lat"="none",
                    "beta.sim" = 0,
                    "beta.sne" = 0,
                    "beta.sor" = 0,
                    "info" = "none"
                ))
        }  
        return(beta_order(res, lat))
    }else{
        return(data.frame(
                    "s_lat"="none",
                    "n_lat"="none",
                    "beta.sim" = 0,
                    "beta.sne" = 0,
                    "beta.sor" = 0,
                    "info" = "none"
                ))
    }
    ## for error skip during loop
}

glm_negative_exp <- function(cene, data_aggr, clim=FALSE){
    da <- subset(data_aggr, data_aggr$info == cene)
    
    if (nrow(da) >= 2){
        da$beta.sim.nor <- (da$beta.sim - mean_(da$beta.sim)) / sd_(da$beta.sim)
        if (clim){
            res <- glm(1- beta.sim+.Machine$double.eps ~ clim_dist , family = gaussian(log), data=da)            
        }
        else{
            res <- glm(1- beta.sim+.Machine$double.eps ~ sp_dist , family = gaussian(log), data=da)            
        }

        temp <- as.data.frame(coef(summary(res)))
        temp$geo_peri <- cene

        return(list("coef"=temp, "pseudo_r2"=1 - res$deviance / res$null.deviance, "AIC"=AIC(res), "beta.mean"=mean_(da$beta.sim), "beta.sd"=sd_(da$beta.sim), "result"=res))        
    }else{
        return(list("coef"=0.0, "pseudo_r2"= 0.0, "AIC"=0.0, "beta.mean"=0.0, "beta.sd"=0.0, "result"="none"))        
    }
}


sepfile <- paste("data/out/site_division_matrix", suppl, ".tsv", sep="")
filename <- sepfile
regions <- read.delim(filename, sep='\t', header=F, stringsAsFactor=FALSE)
    colnames(regions) <- 120:149
    rownames(regions) <- 45:20
    regions$lat <- rownames(regions)
    regions_resh <- melt(regions, id.vars=c("lat"), variable.name="long")
    regions_resh$lat_long <- paste(regions_resh$lat, regions_resh$long, sep="_")

    regions_resh$lat <- as.numeric(regions_resh$lat)
    regions_resh$long <- as.numeric(as.character(regions_resh$long))

# Regression calculation
beta_div_grad <- function(filename=sepfile, CUT=10, removal=0, clim=FALSE) {


    n_regions = 100 # just a random large number
    regions = 0:(n_regions-1) 

    regions_info <- data.frame(matrix(0, nrow=n_regions))
    colnames <- c("region_id")
    regions_info$region_id <- regions
    regions_info <- regions_info[-1]
    for (id in regions){
        regions_info$lat_cent[regions_info$region_id==id] <- mean(subset(regions_resh, regions_resh$value==id)$lat) + 0.5
        regions_info$long_cent[regions_info$region_id==id] <- mean(subset(regions_resh, regions_resh$value==id)$long) + 0.5
    }

    rownames(regions_info) <- regions_info$region_id
    regions_info <- regions_info[-1]

    Oligocene_pa_by_region <- plot_co_region_matrix("Oligocene")
    Miocene_pa_by_region <- plot_co_region_matrix("Miocene")
    Pliocene_pa_by_region <- plot_co_region_matrix("Pliocene")

    Pleistocene_pa_by_region <- plot_co_region_matrix_qua("Pleistocene")
    LastGlacial_pa_by_region <- plot_co_region_matrix_qua("Last glacial period")
    Holocene_pa_by_region <- plot_co_region_matrix_qua("Holocene")

    # Present data processing
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

    rownames(whole_samples_region) <- whole_samples_region$genus
    for_ana_region <- whole_samples_region[-1]

    for_ana_region[for_ana_region > 0] <- 1 # presence/absenceにする

    eliminate_index <- sample(1:nGen,floor(nGen * (1-removal)))
    for_ana_region <- for_ana_region[eliminate_index,]
    for_ana_region <- for_ana_region[apply(for_ana_region, 2, sum)>CUT]
    apply(for_ana_region, 2, sum)
    plot_data_region_dist <- rbind(
    pairwise_beta("Oligocene", for_ana_region),
    pairwise_beta("Miocene", for_ana_region),
    pairwise_beta("Pliocene", for_ana_region),
    pairwise_beta("Pleistocene", for_ana_region),
    pairwise_beta("LastGlacial", for_ana_region),
    pairwise_beta("Holocene", for_ana_region),
    pairwise_beta("Present", for_ana_region)
    ) 
    plot_data_region_dist$info <- factor(plot_data_region_dist$info, levels=c('Oligocene', 'Miocene', 'Pliocene', 'Pleistocene', 'LastGlacial', 'Holocene', 'Present'))
    colnames(plot_data_region_dist) <- c("regionA", "regionB", "beta.sim", "beta.sne",  "beta.sor", "info")
    plot_data_region_dist$regions <- paste(plot_data_region_dist$regionA,   plot_data_region_dist$regionB, sep="_")
    plot_data_region_dist$regionA_cent_lat<- regions_info[plot_data_region_dist$regionA,1]
    plot_data_region_dist$regionA_cent_long<- regions_info[plot_data_region_dist$regionA,2]
    plot_data_region_dist$regionB_cent_lat<- regions_info[plot_data_region_dist$regionB,1]
    plot_data_region_dist$regionB_cent_long<- regions_info[plot_data_region_dist$regionB,2]
    plot_data_region_dist$dist <- 0
    for (k in 1:nrow(plot_data_region_dist)){
        plot_data_region_dist$dist[k] <- distGeo(c(plot_data_region_dist$regionA_cent_long[k],  plot_data_region_dist$regionA_cent_lat[k]), 
                                             c(plot_data_region_dist$regionB_cent_long[k], plot_data_region_dist$regionB_cent_lat[k])) / 1000
    }
    plot_data_region_dist$label <- paste(plot_data_region_dist$regionA,     plot_data_region_dist$regionB, sep="_")
    plot_data_region_dist$info <- factor(plot_data_region_dist$info, levels=c('Oligocene', 'Miocene', 'Pliocene', 'Pleistocene', 'LastGlacial', 'Holocene', 'Present'))
    plot_data_region_dist$sp_dist <- plot_data_region_dist$dist
    
# compile climate data
    if (nrow(plot_data_region_dist) >= 2){
        for (k in 1:nrow(plot_data_region_dist)){
            plot_data_region_dist$clim_dist_Oligo[k] <- abs(clim_resh$Oligo[[as.numeric(plot_data_region_dist$regionA[k])+1,1]]-
            clim_resh$Oligo[[as.numeric(plot_data_region_dist$regionB[k])+1,1]])
            plot_data_region_dist$clim_dist_Mio[k] <- abs(clim_resh$Mio[[as.numeric(plot_data_region_dist$regionA[k])+1,1]] -
            clim_resh$Mio[[as.numeric(plot_data_region_dist$regionB[k])+1,1]])
            plot_data_region_dist$clim_dist_Plio[k] <- abs(clim_resh$Plio[[as.numeric(plot_data_region_dist$regionA[k])+1,1]]-
            clim_resh$Plio[[as.numeric(plot_data_region_dist$regionB[k])+1,1]])
            plot_data_region_dist$clim_dist_Pleisto[k] <- abs(clim_resh$Pleist[[as.numeric(plot_data_region_dist$regionA[k])+1,1]] -
            clim_resh$Pleist[[as.numeric(plot_data_region_dist$regionB[k])+1,1]])
            plot_data_region_dist$clim_dist_LGP[k] <-   abs(clim_resh$LGP[[as.numeric(plot_data_region_dist$regionA[k])+1,1]] -
            clim_resh$LGP[[as.numeric(plot_data_region_dist$regionB[k])+1,1]])
            plot_data_region_dist$clim_dist_Holo[k] <-  abs(clim_resh$Holo[[as.numeric(plot_data_region_dist$regionA[k])+1,1]] -
            clim_resh$Holo[[as.numeric(plot_data_region_dist$regionB[k])+1,1]])
            plot_data_region_dist$clim_dist_Present[k] <-   abs(clim_resh$Present[[as.numeric(plot_data_region_dist$regionA[k])+1,1]] -
            clim_resh$Present[[as.numeric(plot_data_region_dist$regionB[k])+1,1]])
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
    }else{
        return(c("error"))
    }
    
    plot_data_region_dist <- clim_dist_plot_clim_region
    ages <- c('Oligocene', 'Miocene', 'Pliocene', 'Pleistocene', 'LastGlacial', 'Holocene', 'Present')
    df <- data.frame(matrix(rep(NA, 5), nrow=1))[numeric(0), ]
    colnames(df) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)", "geo_peri")
    
    df2 <- data.frame(matrix(rep(NA, 3), nrow=1))[numeric(0), ]
    colnames(df2) <- c("geo_peri", "pseudoR2", "AIC")

    neg_exp_coefs_distonly_region <- df
    neg_exp_fit_measures <- df2
    
    for (age_ in ages){
        aaa <- try(glm_negative_exp(age_, plot_data_region_dist, clim=clim), silent = FALSE)
          if (class(aaa) == "try-error") {
             next
          }

        if (length(aaa$coef) >= 1 ){
            ttemp <- try(rbind(neg_exp_coefs_distonly_region, aaa$coef), silent = FALSE)
             if (class(ttemp) == "try-error") {
                 next
              }
          neg_exp_coefs_distonly_region <- ttemp
        }
        
        tttemp <- data.frame("geo_peri"=age_, "pseudoR2"=aaa$pseudo_r2, "AIC"=aaa$AIC)
          neg_exp_fit_measures <- rbind(neg_exp_fit_measures, tttemp)

    }
    
    neg_exp_coefs_distonly_region$predictor <- gsub("\\d", "", rownames(neg_exp_coefs_distonly_region))
    colnames(neg_exp_coefs_distonly_region) <- gsub("\\s", "",colnames(neg_exp_coefs_distonly_region))
    
    return(list("coefs"=neg_exp_coefs_distonly_region, "fit"=neg_exp_fit_measures))
    
}

# Robust analyses iteration function
removal_iter <- function(removal=0, iter=10, CUT=10, clim=FALSE){
    
    est_pool <- data.frame()
    
    for(i in 1:iter){
        temp <- try(beta_div_grad(removal=removal, CUT=CUT, clim=clim), silent = FALSE)
        #temp <- beta_div_grad(removal=removal, CUT=CUT, clim=clim)
        if (class(temp) == "try-error") {
            next
        }
        if (ncol(temp$coefs) < 4){
            next
        }
        temp_estimates <- temp$coefs$Estimate#[temp$predictor!="(Intercept)"]
        temp2 <- data.frame("Estimate"=temp_estimates)
        temp2$sd <- temp$coefs[,"Std.Error"]
        temp2$p.value <- temp$coefs[,"Pr(>|t|)"]
        temp2$geo_peri <- factor(temp$coefs$geo_peri, levels=c('Oligocene', 'Miocene', 'Pliocene', 'Pleistocene', 'LastGlacial', 'Holocene', 'Present'))
        temp2$predictor <- temp$coefs$predictor
        temp2$r2 <- as.vector(t(cbind(temp$fit$pseudoR2, temp$fit$pseudoR2)))

        if (nrow(est_pool)==0){
            est_pool <- temp2
        }
        else{
            est_pool <- rbind(est_pool, temp2)
        }
    }
    return(est_pool)
}


# Coefficient reshape func
table_resh <- function(result, CUT, removal, clim=FALSE){
    icp <- show_coefs(subset(result$coefs, predictor=="(Intercept)"))
    slp <- show_coefs(subset(result$coefs, predictor!="(Intercept)"))
    temp <- data.frame(t(test_dicay_coef_$fit))
    colnames(temp) <- temp[1,]
    temp <- temp[c(2,3),]
    temp <- data.frame(t(temp))
    temp$pseudoR2 <- as.numeric(temp$pseudoR2)
    temp$AIC <- as.numeric(temp$AIC)
    temp <- t(temp)
    temp
    regTable <- rbind(icp, slp, temp)
    if(clim){
        clim = "clim"
    }else{
        clim=""
    }
    write.csv(regTable, paste("regrTable_fromR/CUT", CUT, "removal", removal, clim, nre,".csv", sep=""),row.names=TRUE,col.names=TRUE, append=FALSE)
    return(regTable)
}
regions = 0:(n_regions-1)

# Genera remove robust analyses
test_dicay_coef_25 <- removal_iter(removal=0.25, iter=999)

saveRDS(test_dicay_coef_25, "regr_robust.rds")

## Figure A6(a)
ga <- subset(test_dicay_coef_25, test_dicay_coef_25$predictor=="sp_dist" )%>%
ggplot(aes(x=geo_peri, y=-Estimate*1000, fill=predictor)) + geom_boxplot(color="black", fill="grey") + ylim(-0.5, 0.3) + geom_hline(yintercept = 0, color = 'blue',
linetype = 'dotted')+theme_bw()+theme(text = element_text(size = 20))+xlab("")+ylab("spatial turnover (×1000)")
ggsave(plot=ga, file = "data/out/fig_a6a.pdf", width=10, height=10)
ggsave(plot=ga, file = "data/out/fig_a6a.png", width=10, height=10, dpi=600)


## Figure A6(b)
ga <- subset(test_dicay_coef_25, test_dicay_coef_25$predictor=="sp_dist" )%>%
ggplot(aes(x=geo_peri, y=r2, fill=predictor)) + geom_boxplot(color="black", fill="grey") + ylim(-0.1, 0.85) + geom_hline(yintercept = 0, color = 'blue',
linetype = 'dotted')+theme_bw()+theme(text = element_text(size = 20))+xlab("")+ylab("pseudo R2")
ggsave(plot=ga, file = "data/out/fig_a6b.pdf", width=10, height=10)
ggsave(plot=ga, file = "data/out/fig_a6b.png", width=10, height=10, dpi=600)


test_dicay_coef_25_clim <- removal_iter(removal=0.25, iter=999, clim=TRUE)
saveRDS(test_dicay_coef_25_clim, "regr_robust_clim.rds")


## Figure A6(c)
ga <- subset(test_dicay_coef_25_clim, test_dicay_coef_25_clim$predictor=="clim_dist" )%>%
ggplot(aes(x=geo_peri, y=-Estimate, fill=predictor)) + geom_boxplot(color="black", fill="grey") + ylim(-0.05, 0.03) + geom_hline(yintercept = 0, color = 'blue',
 linetype = 'dotted')+theme_bw()+theme(text = element_text(size = 20))+xlab("")+ylab("climatic turnover")
ggsave(plot=ga, file = "data/out/fig_a6c.pdf", width=10, height=10)
ggsave(plot=ga, file = "data/out/fig_a6c.png", width=10, height=10, dpi=600)

 ## Figure A6(d)
ga <- subset(test_dicay_coef_25_clim, test_dicay_coef_25_clim$predictor=="clim_dist" )%>%
ggplot(aes(x=geo_peri, y=r2, fill=predictor))  + geom_boxplot(color="black", fill="grey") + ylim(-0.1, 0.5) + geom_hline(yintercept = 0, color = 'blue',
linetype = 'dotted')+theme_bw()+theme(text = element_text(size = 20))+xlab("")+ylab("pseudo R2")
ggsave(plot=ga, file = "data/out/fig_a6d.pdf", width=10, height=10)
ggsave(plot=ga, file = "data/out/fig_a6d.png", width=10, height=10, dpi=600)




