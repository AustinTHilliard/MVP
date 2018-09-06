setwd("~/Documents/lipids/")
options(stringsAsFactors=FALSE)
rm(list=ls())

# ------------------------------------------------------------------------------
# -- load and prepare PAGE data
# ------------------------------------------------------------------------------
# load data and add column to match format of MVP snps
snps_page = read.csv("PAGE_Lipids_Replication_Table.csv")
snps_page$MarkerName = paste0(snps_page$CHR, ":", snps_page$POS_hg19)

# load proxy snps and add columns to match format of MVP snps
snps_page_proxies = read.csv("PAGE_Lipids_Replication_Proxies.csv")
snps_page_proxies$Index.MarkerName = paste0(snps_page_proxies$Chr_indexSNP, ":", 
                                            snps_page_proxies$Pos_indexSNP)
snps_page_proxies$Proxy.MarkerName = paste0(snps_page_proxies$Chr_proxy, ":", 
                                            snps_page_proxies$Pos_proxy)

# get names of all PAGE snps, including proxies
c1 = all(snps_page$MarkerName %in% 
         snps_page_proxies$Index.MarkerName)
c2 = all(snps_page_proxies$Index.MarkerName %in% 
         snps_page_proxies$Proxy.MarkerName)
if (c1 && c2) {
    all_snps_page = unique(c(snps_page_proxies$Proxy.MarkerName))
}
rm(c1, c2)

# ------------------------------------------------------------------------------                                       
# -- MVP - load and filter ethnicity-specific data
# ------------------------------------------------------------------------------
# get names of files downloaded from Box:
#   MVP GWAS/Lipids GWAS/1000G_GWAS_OUT_GENESIS_KLARIN/
mvp_files = grep("^allchr.*gz$", list.files(), value=T)

for (file in mvp_files) {
    xlist = list()
    cat(paste0("---- ", Sys.time(), "\n", file), "...\n")
    line_num = 0
    for (xline in readLines(file)) {
        if (line_num %% 1e6 == 0) 
            cat(line_num, "... ")
        xline_split = strsplit(xline, " ")[[1]]
        snp = xline_split[1]
        if (snp %in% all_snps_page) 
            xlist[[snp]] = xline_split
        line_num = line_num + 1
    }
    xdf = as.data.frame(do.call(rbind, xlist))
    rownames(xdf) = NULL
    names(xdf) = strsplit(readLines(file, n=1), " ")[[1]]
    saveRDS(xdf, file=paste0(file, "_PAGE.rds"))
    cat("\n")
    rm(line_num, xline, xline_split, snp, xdf)
    gc()
}
rm(file)
# 
# ------------------------------------------------------------------------------                                       
# ---- 2018-09-05 17:32:01
# allchr..EUR.LDL.results_10_21_17.gz ...
# 0 ... 1e+06 ... 2e+06 ... 3e+06 ... 4e+06 ... 5e+06 ... 6e+06 ... 7e+06 ... 8e+06 ... 9e+06 ... 1e+07 ... 1.1e+07 ... 1.2e+07 ... 1.3e+07 ... 1.4e+07 ... 1.5e+07 ... 1.6e+07 ... 1.7e+07 ... 1.8e+07 ... 1.9e+07 ... 2e+07 ... 
# ---- 2018-09-05 17:44:17
# allchr.AFR.HDL.results_10_31_17.gz ...
# 0 ... 1e+06 ... 2e+06 ... 3e+06 ... 4e+06 ... 5e+06 ... 6e+06 ... 7e+06 ... 8e+06 ... 9e+06 ... 1e+07 ... 1.1e+07 ... 1.2e+07 ... 1.3e+07 ... 1.4e+07 ... 1.5e+07 ... 1.6e+07 ... 1.7e+07 ... 1.8e+07 ... 1.9e+07 ... 2e+07 ... 2.1e+07 ... 2.2e+07 ... 2.3e+07 ... 2.4e+07 ... 2.5e+07 ... 2.6e+07 ... 2.7e+07 ... 2.8e+07 ... 2.9e+07 ... 3e+07 ... 3.1e+07 ... 3.2e+07 ... 
# ---- 2018-09-05 18:04:20
# allchr.AFR.LDL.results_10_31_17.gz ...
# 0 ... 1e+06 ... 2e+06 ... 3e+06 ... 4e+06 ... 5e+06 ... 6e+06 ... 7e+06 ... 8e+06 ... 9e+06 ... 1e+07 ... 1.1e+07 ... 1.2e+07 ... 1.3e+07 ... 1.4e+07 ... 1.5e+07 ... 1.6e+07 ... 1.7e+07 ... 1.8e+07 ... 1.9e+07 ... 2e+07 ... 2.1e+07 ... 2.2e+07 ... 2.3e+07 ... 2.4e+07 ... 2.5e+07 ... 2.6e+07 ... 2.7e+07 ... 2.8e+07 ... 2.9e+07 ... 3e+07 ... 3.1e+07 ... 3.2e+07 ... 
# ---- 2018-09-05 18:23:25
# allchr.AFR.TC.results_10_31_17.gz ...
# 0 ... 1e+06 ... 2e+06 ... 3e+06 ... 4e+06 ... 5e+06 ... 6e+06 ... 7e+06 ... 8e+06 ... 9e+06 ... 1e+07 ... 1.1e+07 ... 1.2e+07 ... 1.3e+07 ... 1.4e+07 ... 1.5e+07 ... 1.6e+07 ... 1.7e+07 ... 1.8e+07 ... 1.9e+07 ... 2e+07 ... 2.1e+07 ... 2.2e+07 ... 2.3e+07 ... 2.4e+07 ... 2.5e+07 ... 2.6e+07 ... 2.7e+07 ... 2.8e+07 ... 2.9e+07 ... 3e+07 ... 3.1e+07 ... 3.2e+07 ... 
# ---- 2018-09-05 18:43:07
# allchr.AFR.TG.results_10_31_17.gz ...
# 0 ... 1e+06 ... 2e+06 ... 3e+06 ... 4e+06 ... 5e+06 ... 6e+06 ... 7e+06 ... 8e+06 ... 9e+06 ... 1e+07 ... 1.1e+07 ... 1.2e+07 ... 1.3e+07 ... 1.4e+07 ... 1.5e+07 ... 1.6e+07 ... 1.7e+07 ... 1.8e+07 ... 1.9e+07 ... 2e+07 ... 2.1e+07 ... 2.2e+07 ... 2.3e+07 ... 2.4e+07 ... 2.5e+07 ... 2.6e+07 ... 2.7e+07 ... 2.8e+07 ... 2.9e+07 ... 3e+07 ... 3.1e+07 ... 3.2e+07 ... 
# ---- 2018-09-05 19:03:37
# allchr.EUR.HDL.results_10_31_17.gz ...
# 0 ... 1e+06 ... 2e+06 ... 3e+06 ... 4e+06 ... 5e+06 ... 6e+06 ... 7e+06 ... 8e+06 ... 9e+06 ... 1e+07 ... 1.1e+07 ... 1.2e+07 ... 1.3e+07 ... 1.4e+07 ... 1.5e+07 ... 1.6e+07 ... 1.7e+07 ... 1.8e+07 ... 1.9e+07 ... 2e+07 ... 
# ---- 2018-09-05 19:15:02
# allchr.EUR.TC.results_10_31_17.gz ...
# 0 ... 1e+06 ... 2e+06 ... 3e+06 ... 4e+06 ... 5e+06 ... 6e+06 ... 7e+06 ... 8e+06 ... 9e+06 ... 1e+07 ... 1.1e+07 ... 1.2e+07 ... 1.3e+07 ... 1.4e+07 ... 1.5e+07 ... 1.6e+07 ... 1.7e+07 ... 1.8e+07 ... 1.9e+07 ... 
# ---- 2018-09-05 19:26:36
# allchr.EUR.TG.results_10_31_17.gz ...
# 0 ... 1e+06 ... 2e+06 ... 3e+06 ... 4e+06 ... 5e+06 ... 6e+06 ... 7e+06 ... 8e+06 ... 9e+06 ... 1e+07 ... 1.1e+07 ... 1.2e+07 ... 1.3e+07 ... 1.4e+07 ... 1.5e+07 ... 1.6e+07 ... 1.7e+07 ... 1.8e+07 ... 1.9e+07 ... 2e+07 ... 
# ---- 2018-09-05 19:38:33
# allchr.HIS.HDL.results_10_31_17.gz ...
# 0 ... 1e+06 ... 2e+06 ... 3e+06 ... 4e+06 ... 5e+06 ... 6e+06 ... 7e+06 ... 8e+06 ... 9e+06 ... 1e+07 ... 1.1e+07 ... 1.2e+07 ... 1.3e+07 ... 1.4e+07 ... 1.5e+07 ... 1.6e+07 ... 1.7e+07 ... 1.8e+07 ... 1.9e+07 ... 2e+07 ... 2.1e+07 ... 2.2e+07 ... 2.3e+07 ... 2.4e+07 ... 2.5e+07 ... 2.6e+07 ... 2.7e+07 ... 2.8e+07 ... 2.9e+07 ... 3e+07 ... 3.1e+07 ... 3.2e+07 ... 
# ---- 2018-09-05 19:58:14
# allchr.HIS.LDL.results_10_31_17.gz ...
# 0 ... 1e+06 ... 2e+06 ... 3e+06 ... 4e+06 ... 5e+06 ... 6e+06 ... 7e+06 ... 8e+06 ... 9e+06 ... 1e+07 ... 1.1e+07 ... 1.2e+07 ... 1.3e+07 ... 1.4e+07 ... 1.5e+07 ... 1.6e+07 ... 1.7e+07 ... 1.8e+07 ... 1.9e+07 ... 2e+07 ... 2.1e+07 ... 2.2e+07 ... 2.3e+07 ... 2.4e+07 ... 2.5e+07 ... 2.6e+07 ... 2.7e+07 ... 2.8e+07 ... 2.9e+07 ... 3e+07 ... 3.1e+07 ... 3.2e+07 ... 
# ---- 2018-09-05 20:16:27
# allchr.HIS.TC.results_10_31_17.gz ...
# 0 ... 1e+06 ... 2e+06 ... 3e+06 ... 4e+06 ... 5e+06 ... 6e+06 ... 7e+06 ... 8e+06 ... 9e+06 ... 1e+07 ... 1.1e+07 ... 1.2e+07 ... 1.3e+07 ... 1.4e+07 ... 1.5e+07 ... 1.6e+07 ... 1.7e+07 ... 1.8e+07 ... 1.9e+07 ... 2e+07 ... 2.1e+07 ... 2.2e+07 ... 2.3e+07 ... 2.4e+07 ... 2.5e+07 ... 2.6e+07 ... 2.7e+07 ... 2.8e+07 ... 2.9e+07 ... 3e+07 ... 3.1e+07 ... 3.2e+07 ... 
# ---- 2018-09-05 20:35:00
# allchr.HIS.TG.results_10_31_17.gz ...
# 0 ... 1e+06 ... 2e+06 ... 3e+06 ... 4e+06 ... 5e+06 ... 6e+06 ... 7e+06 ... 8e+06 ... 9e+06 ... 1e+07 ... 1.1e+07 ... 1.2e+07 ... 1.3e+07 ... 1.4e+07 ... 1.5e+07 ... 1.6e+07 ... 1.7e+07 ... 1.8e+07 ... 1.9e+07 ... 2e+07 ... 2.1e+07 ... 2.2e+07 ... 2.3e+07 ... 2.4e+07 ... 2.5e+07 ... 2.6e+07 ... 2.7e+07 ... 2.8e+07 ... 2.9e+07 ... 3e+07 ... 3.1e+07 ... 3.2e+07 ... 

# ------------------------------------------------------------------------------                                       
# -- MVP - define functions for matching snps and grouping the results
# ------------------------------------------------------------------------------
# assumes that snps_page exists in workspace
.matchToPAGE = function (xdf) {
    xdf_matched = xdf[match(snps_page$MarkerName, xdf$rsid), ]
    xdf_matched = as.data.frame(cbind(snps_page, xdf_matched[, -c(1:3)]))
    rownames(xdf_matched) = NULL
    return(xdf_matched)
}

# assumes that xlist_matched and xtraits exist in workspace
.groupByTrait = function (trait) {
    xlist_matched_trait = xlist_matched[xtraits == trait]
    xlist_matched_trait = lapply(xlist_matched_trait, 
                                 function(f) subset(f, Trait == trait))
    xdf_matched_trait = as.data.frame(do.call(rbind, xlist_matched_trait))
    split_rownames = strsplit(rownames(xdf_matched_trait), ".", fixed=T)
    xdf_matched_trait$race = sapply(split_rownames, function(f) f[2])
    names(xdf_matched_trait)[6:19] = paste0("MVP.", 
                                            names(xdf_matched_trait)[6:19])
    xdf_matched_trait = xdf_matched_trait[order(xdf_matched_trait$rsID), 
                                          c(1:5, 19, 6:18)]
    rownames(xdf_matched_trait) = NULL
    return(xdf_matched_trait)
}

# assumes miss_snps_proxies, xtraits, and xlist exist in workspace
.getProxyForSnp = function (missing_snp) {
    snp_rows = miss_snps_proxies$Index.MarkerName == missing_snp
    prox = subset(miss_snps_proxies, snp_rows)$Proxy.MarkerName
    prox = prox[prox != missing_snp]
    trait = unique(subset(miss_snps_proxies, snp_rows)$Trait.s.)
    if (grepl("_", trait)) {
        trait = unlist(strsplit(trait, "_"))
    }
    mvp_trait = xlist[xtraits %in% trait]
    mvp_proxies = list()
    for (race_eth in 1:length(mvp_trait)) {
        xmvp = mvp_trait[[race_eth]]
        xmvp_name = names(mvp_trait)[race_eth]
        xmvp = xmvp[match(prox, xmvp$rsid), ]
        xmvp$rsid = prox
        if (!all(is.na(xmvp))) {
            NArows = apply(xmvp, 1, function(f) all(is.na(f)))
            xmvp = xmvp[!NArows, ]
        }
        mvp_proxies[[xmvp_name]] = xmvp
    }
    mvp_proxies = as.data.frame(do.call(rbind, mvp_proxies))
    spl_rownames = strsplit(rownames(mvp_proxies), "\\.")
    
    mvp_proxies$Trait = sapply(spl_rownames, function(f) f[3])
    mvp_proxies$Index.MarkerName = rep(missing_snp, nrow(mvp_proxies))
    mvp_proxies$MVP.race = sapply(spl_rownames, function(f) f[2])
    
    mvp_proxies = mvp_proxies[, c(17:19, 1:16)]
    names(mvp_proxies)[-c(1:3)] = paste0("Proxy.MVP.", 
                                         names(mvp_proxies)[-c(1:3)])
    rownames(mvp_proxies) = NULL
    return(list(Proxy.MarkerName=prox, Trait=trait, MVP=mvp_proxies))
}

# ------------------------------------------------------------------------------                                       
# -- MVP - match and merge with primary PAGE snps and group results by trait
# ------------------------------------------------------------------------------
# get filenames
mvp_page_files = paste0(mvp_files, "_PAGE.rds")

# reload MVP data matched to PAGE snps (including proxies)
if (all(mvp_page_files %in% list.files())) {
    xlist = lapply(as.list(mvp_page_files), readRDS)
    names(xlist) = mvp_page_files
}

# fix names with two dots
names(xlist) = gsub("..", ".", names(xlist), fixed=T)

# match MVP results to just the primary PAGE snps
xlist_matched = lapply(xlist, .matchToPAGE)

# define the trait represented by each dataset
traits = c("HDL", "LDL", "TC", "TG")
xtraits = sapply(as.list(names(xlist)), 
                 function(ff) traits[sapply(as.list(traits), 
                                            function(f) grepl(f, ff))])

# group the matched PAGE snps by trait
xlist_matched_traits = lapply(as.list(traits), .groupByTrait)
names(xlist_matched_traits) = traits

# collapse resulting list to a single data frame
xdf_matched_traits = as.data.frame(do.call(rbind, xlist_matched_traits))
rownames(xdf_matched_traits) = NULL

# ------------------------------------------------------------------------------                                       
# -- get proxies for PAGE snps missing from MVP
# ------------------------------------------------------------------------------
# get names of missing snps and their proxies
miss_snps = unique(subset(xdf_matched_traits, is.na(MVP.alleleA))$MarkerName)
miss_snps_proxies = subset(snps_page_proxies, Index.MarkerName %in% miss_snps)

#
mvp_proxies = lapply(as.list(miss_snps), .getProxyForSnp)
names(mvp_proxies) = miss_snps

for (m in 1:length(mvp_proxies)) {
    if (nrow(mvp_proxies[[m]]$MVP) == 0) {
        mvp_proxies[[m]]$MVP = mvp_proxies[[m]]$MVP[1, ]
        mvp_proxies[[m]]$MVP$Trait = mvp_proxies[[m]]$Trait
        mvp_proxies[[m]]$MVP$Index.MarkerName = names(mvp_proxies)[m]
        
        x = subset(xdf_matched_traits, MarkerName == names(mvp_proxies)[m])
        xNA = is.na(x$MVP.alleleA)
        if (sum(xNA) > 1) {
            mvp_proxies[[m]]$MVP = mvp_proxies[[m]]$MVP[1:sum(xNA), ]
            mvp_proxies[[m]]$MVP[, 1:3] = x[, c(1,5,6)]
        } else {
            mvp_proxies[[m]]$MVP$MVP.race = x$MVP.race[xNA]
        }
    }
}
rm(m, x, xNA)

#
mvp_proxies_df = as.data.frame(do.call(rbind, lapply(mvp_proxies, 
                                                     function(f) f$MVP)))
rownames(mvp_proxies_df) = NULL

#
flt = subset(xdf_matched_traits, MarkerName %in% miss_snps & is.na(MVP.alleleA))
flt = apply(flt[, c("Trait", "MarkerName", "MVP.race")], 1, 
            paste0, collapse="_")

tocheck = apply(mvp_proxies_df[, c("Trait", "Index.MarkerName", "MVP.race")], 1,
                paste0, collapse="_")

mvp_proxies_df = mvp_proxies_df[tocheck %in% flt, ]
    
# ------------------------------------------------------------------------------                                       
# -- write output, then save workspace with session info
# ------------------------------------------------------------------------------
# MVP results for primary PAGE snps
write.csv(xdf_matched_traits, row.names=F,
          file="PAGE_Lipids_Replication_MVP.ethnic-specific_noProxies.csv")

# MVP results for proxies of missing PAGE snps
write.csv(mvp_proxies_df, row.names=F,
          file="PAGE_Lipids_Replication_MVP.ethnic-specific_Proxies_for_missing.csv")

# get session info
session_info = sessionInfo()

# workspace 
save(list=ls(),
     file="lipids_replication_ethnic-specific_WORKSPACE.RData")

