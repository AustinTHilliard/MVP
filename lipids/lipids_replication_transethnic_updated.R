setwd("~/Documents/lipids/")
options(stringsAsFactors=FALSE)
rm(list=ls())

# ------------------------------------------------------------------------------
# -- load and prepare PAGE data
# ------------------------------------------------------------------------------
# load data and add column to match format of MVP snps
snps_page = read.csv("PAGE_materials/PAGE_Lipids_Replication_Table_updated.csv")
dim(snps_page)
# [1] 54  4

snps_page$MarkerName = paste0(snps_page$CHR, ":", snps_page$POS_hg19)

# split snps by trait
xpage = split(snps_page, snps_page$Trait)
names(xpage)
# [1] "HDL" "LDL" "TC"  "TG"
sapply(xpage, nrow)
# HDL LDL  TC  TG 
#  18   8  14  14

# load proxy snps and add columns to match format of MVP snps
snps_page_proxies = read.csv("PAGE_materials/PAGE_Lipids_Replication_Proxies_updated.csv")
dim(snps_page_proxies)
# [1] 716   8

# order by trait and index SNP position 
xord = order(snps_page_proxies$Trait.s., 
             snps_page_proxies$Chr_indexSNP,
             snps_page_proxies$Pos_indexSNP)
snps_page_proxies = snps_page_proxies[xord, ]
rm(xord)

snps_page_proxies$Index.MarkerName = paste0(snps_page_proxies$Chr_indexSNP, ":", 
                                            snps_page_proxies$Pos_indexSNP)
snps_page_proxies$Proxy.MarkerName = paste0(snps_page_proxies$Chr_proxy, ":", 
                                            snps_page_proxies$Pos_proxy)

# ------------------------------------------------------------------------------                                       
# -- MVP - load transethnic data
# ------------------------------------------------------------------------------
# get names of files downloaded from Box:
#   MVP-CAP other local/Manuscripts/MVP lipids GWAS/Summary_Statistics/
mvp_te_files = grep("^MVP.*gz$", list.files("MVP_materials/"), value=T)
mvp_te_files
# [1] "MVP.te.HDL.gwas.tsv.gz" "MVP.te.LDL.gwas.tsv.gz" "MVP.te.TC.gwas.tsv.gz"  "MVP.te.TG.gwas.tsv.gz" 

# load the data 
#  will take a minute or two for each file since they're each ~570mb
#  if on a machine with not much memory should use readLines approach instead
#   e.g. as in lipids_replication_wthnic-specific.R
for (file in paste0("MVP_materials/", mvp_te_files)) {
    cat(file, "...\n")
    assign(gsub("MVP_materials/", "", file), 
           read.table(file, header=T, as.is=T))
}
rm(file)

gc()
ls()
# [1] "mvp_te_files" 
#     "MVP.te.HDL.gwas.tsv.gz" "MVP.te.LDL.gwas.tsv.gz" "MVP.te.TC.gwas.tsv.gz" "MVP.te.TG.gwas.tsv.gz"
#     "snps_page" "snps_page_proxies" "xpage" 

sapply(as.list(mvp_te_files), function(f) dim(get(f)))
#          [,1]     [,2]     [,3]     [,4]
# [1,] 28946834 29221975 29085659 29184545
# [2,]        9        9        9        9

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

sapply(mvp_te_page, nrow)
# HDL LDL  TC  TG 
#  18   8  14  14 

# collapse to one data frame
mvp_te_page = as.data.frame(do.call(rbind, mvp_te_page))
rownames(mvp_te_page) = NULL

# ------------------------------------------------------------------------------                                       
# -- get proxies for missing snps
# ------------------------------------------------------------------------------
# get names of missing snps and their proxies
miss_snps = unique(subset(mvp_te_page, is.na(MVP.te_Allele1))$MarkerName)

miss_snps_proxies = subset(snps_page_proxies, Index.MarkerName %in% miss_snps)
dim(miss_snps_proxies)
# [1] 23 10

# some of the missing SNPS only have themselves listed as a proxy
sapply(split(miss_snps_proxies, miss_snps_proxies$Index.MarkerName), 
       function(f) 
           all(f$Index.MarkerName == f$Proxy.MarkerName))

# ------------------------------------------------------------------------------                                       
# -- find proxies in MVP data
# ------------------------------------------------------------------------------
# assumes that MVP data frames exist in the workspace:
#  MVP.te.HDL.gwas.tsv.gz, MVP.te.LDL.gwas.tsv.gz, MVP.te.TC.gwas.tsv.gz, MVP.te.TG.gwas.tsv.gz
.getProxiesInMVPte = function (miss_snps=NULL, miss_snps_proxies=NULL) {
    # create list to store results and loop through missing SNPs
    mvp_proxies = list()
    for (s in miss_snps) {
        
        # get proxies and trait for this missing SNP
        xprox0 = subset(miss_snps_proxies, Index.MarkerName == s)
        xprox = xprox0$Proxy.MarkerName
        xtrait = unique(xprox0$Trait.s.)

        # get appropriate set of MVP results from the workspace
        xname = grep(paste0("^MVP.*", xtrait), ls(name=".GlobalEnv"), value=T)
        xmvp = get(xname)
        
        # print updates and sanity check
        cat(paste0("--", s, " missing from ", xtrait, " results\n",
                   "  checking proxy SNPs ", paste0(xprox, collapse=", "), "\n", 
                   "  in ", xname, "\n\n"))
        
        # subset MVP results to the proxy SNPs
        tmp_mvp = xmvp[match(xprox, xmvp$MarkerName), ]
        
        # make sure missing proxies have IDs in the output
        tmp_mvp$MarkerName = xprox
        
        # add columns for missing index SNP and trait
        tmp_mvp$Index.MarkerName = rep(s, nrow(tmp_mvp))
        tmp_mvp$Trait = rep(xtrait, nrow(tmp_mvp))
        
        # remove rows where missing SNP is listed as its own proxy
        #  unless only proxy is itself, i.e. nrow(tmp_mvp) == 1
        if (nrow(tmp_mvp) > 1) {
            sameRow = tmp_mvp$Index.MarkerName == tmp_mvp$MarkerName
            tmp_mvp = tmp_mvp[!sameRow, ]
        }
        
        mvp_proxies[[s]] = tmp_mvp
    }
    
    # collapse to single data frame, reorder columns, set dimnames
    mvp_proxies = as.data.frame(do.call(rbind, mvp_proxies))
    colNum = ncol(mvp_proxies)
    colOrd = c(colNum,   # should be Trait
               colNum-1, # should be Index.MarkerName
               1:(colNum-2))
    mvp_proxies = mvp_proxies[, colOrd]
    names(mvp_proxies)[-c(1:2)] = paste0("Proxy.", 
                                         names(mvp_proxies)[-c(1:2)])
    rownames(mvp_proxies) = NULL
    return(mvp_proxies)
}

mvp_proxies = .getProxiesInMVPte(miss_snps, miss_snps_proxies)

# ------------------------------------------------------------------------------                                       
# -- plot distributions of MVP p-values for each trait
# ------------------------------------------------------------------------------
par(mfrow=c(2,2))
xcuts = c(0, 1e-8, 1e-6, 1e-4, 1e-2, 1e-1, 1)
for (tr in c("HDL","TC","LDL","TG")) {
    barplot(table(cut(as.numeric(mvp_te_page$MVP.te_P.value)[mvp_te_page$Trait == tr], xcuts)),
            main=tr, col='lightgrey')
    abline(h=5, col='darkgrey')
}

# ------------------------------------------------------------------------------                                       
# -- write output, then save workspace with session info
# ------------------------------------------------------------------------------
# MVP results for PAGE snps
write.csv(mvp_te_page, row.names=F,
          file="PAGE_Lipids_Replication_MVP.te_noProxies_updated.csv")

# MVP results for proxies of PAGE snps missing from MVP
write.csv(mvp_proxies, row.names=F,
          file="PAGE_Lipids_Replication_MVP.te_Proxies_for_missing_updated.csv")

# csv files are combined in PAGE_Lipids_Replication_MVP.transethnic_updated.xlsx

# get session and workspace info
session_info = sessionInfo()
full_workspace = ls()

# workspace without full MVP results
save(list=ls()[-match(mvp_te_files, ls())],
     file="lipids_replication_transethnic_WORKSPACE-no_full_MVP_updated.RData")

# workspace including full MVP results
save(list=ls(),
     file="lipids_replication_transethnic_WORKSPACE-with_full_MVP_updated.RData")
