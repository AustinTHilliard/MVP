setwd("~/Documents/lipids/")
options(stringsAsFactors=FALSE)
rm(list=ls())

# ------------------------------------------------------------------------------
# -- load and prepare PAGE data
# ------------------------------------------------------------------------------
# load data and add column to match format of MVP snps
snps_page = read.csv("PAGE_Lipids_Replication_Table.csv")
snps_page$MarkerName = paste0(snps_page$CHR, ":", snps_page$POS_hg19)
dim(snps_page)
# [1] 36  5

# split snps by trait
xpage = split(snps_page, snps_page$Trait)
names(xpage)
# [1] "HDL" "LDL" "TC"  "TG"
sapply(xpage, nrow)
# HDL LDL  TC  TG 
#  14   7  10   5 

# load proxy snps and add columns to match format of MVP snps
snps_page_proxies = read.csv("PAGE_Lipids_Replication_Proxies.csv")
snps_page_proxies$Index.MarkerName = paste0(snps_page_proxies$Chr_indexSNP, ":", 
                                            snps_page_proxies$Pos_indexSNP)
snps_page_proxies$Proxy.MarkerName = paste0(snps_page_proxies$Chr_proxy, ":", 
                                            snps_page_proxies$Pos_proxy)
dim(snps_page_proxies)
# [1] 585  10

# ------------------------------------------------------------------------------                                       
# -- MVP - load transethnic data
# ------------------------------------------------------------------------------
# get names of files downloaded from Box:
#   MVP-CAP other local/Manuscripts/MVP lipids GWAS/Summary_Statistics/
mvp_te_files = grep("^MVP.*gz$", list.files(), value=T)

# load the data
for (file in mvp_te_files) {
    cat(file, "...\n")
    assign(file, read.table(file, header=T, as.is=T))
}
rm(file)

gc()
ls()
# [1] "mvp_te_files" "MVP.te.HDL.gwas.tsv.gz" "MVP.te.LDL.gwas.tsv.gz" 
#     "MVP.te.TC.gwas.tsv.gz" "MVP.te.TG.gwas.tsv.gz" "snps_page" "xpage"

# ------------------------------------------------------------------------------                                       
# -- MVP - work on transethnic data
# ------------------------------------------------------------------------------
# make a list to hold the results, then loop through the traits
mvp_te_page = list()
for (trait in names(xpage)) {
    cat(trait, "...\n")
    
    # get the appropriate MVP dataset and subset based on PAGE snps
    x = get(grep(paste0("^MVP.*", trait), ls(), value=T))
    x = subset(x, MarkerName %in% xpage[[trait]]$MarkerName)
    
    # use match to insert NA rows where MVP data is missing a PAGE snp
    x = x[match(xpage[[trait]]$MarkerName, x$MarkerName), ]
    
    # combine the MVP and PAGE information
    names(x) = paste0("MVP.te_", names(x))
    mvp_te_page[[trait]] = cbind(xpage[[trait]], x[, -1])
}
rm(trait, x)
gc()

# collapse to one data frame
mvp_te_page = as.data.frame(do.call(rbind, mvp_te_page))

# ------------------------------------------------------------------------------                                       
# -- get proxies for missing snps
# ------------------------------------------------------------------------------
# get names of missing snps and their proxies
miss_snps = unique(subset(mvp_te_page, is.na(MVP.te_Allele1))$MarkerName)
miss_snps_proxies = subset(snps_page_proxies, Index.MarkerName %in% miss_snps)

mvp_proxies = list()
for (s in miss_snps) {
    xprox = subset(miss_snps_proxies, Index.MarkerName == s)$Proxy.MarkerName
    xtrait = unique(subset(miss_snps_proxies, Index.MarkerName==s)$Trait.s.)
    xmvp = get(grep(paste0("^MVP.*", xtrait), ls(), value=T))
    mvp_proxies[[s]] = xmvp[match(xprox, xmvp$MarkerName), ]
    if (!all(is.na(mvp_proxies[[s]]))) {
        NArows = apply(mvp_proxies[[s]], 1, function(f) all(is.na(f)))
        mvp_proxies[[s]] = mvp_proxies[[s]][!NArows, ]
    }
}
rm(s, xprox, xtrait, xmvp, NArows)

mvp_proxies = as.data.frame(do.call(rbind, mvp_proxies))
mvp_proxies$Index.MarkerName = sapply(strsplit(rownames(mvp_proxies), "\\."), 
                                      function(f) f[1])

mvp_proxies = mvp_proxies[, c(ncol(mvp_proxies), 1:(ncol(mvp_proxies)-1))]
names(mvp_proxies)[-1] = paste0("Proxy.", names(mvp_proxies)[-1])

# ------------------------------------------------------------------------------                                       
# -- write output, then save workspace with session info
# ------------------------------------------------------------------------------
# MVP results for PAGE snps
write.csv(mvp_te_page, row.names=F,
          file="PAGE_Lipids_Replication_MVP.te_noProxies.csv")

# MVP results for proxies of PAGE snps missing from MVP
write.csv(mvp_proxies, row.names=F,
          file="PAGE_Lipids_Replication_MVP.te_Proxies_for_missing.csv")

# get session info
session_info = sessionInfo()

# workspace including full MVP results
save(list=ls(),
     file="lipids_replication_transethnic_WORKSPACE-with_full_MVP.RData")

# workspace without full MVP results
save(list=ls()[-match(mvp_te_files, ls())],
     file="lipids_replication_transethnic_WORKSPACE-no_full_MVP.RData")
