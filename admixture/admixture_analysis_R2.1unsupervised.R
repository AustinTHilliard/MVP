# set working directory to where the data is
dir_dat = "/Users/ath/Dropbox/_PAVIR/admixture/"
setwd(dir_dat)

# set file names
file_race = "R2.1_geno.01_maf.05_hwe.000001.midp_snps.acgt_excl.high.LD.reg_indep-pair.1000.100..03_rsid.fam_race-eth_NO.mvp001_id.txt"
files_subjects = grep("Q$", list.files(), value=TRUE)
files_variants = grep("P$", list.files(), value=TRUE)

# load data
# race/ethnicity file
rac = read.table(file_race, header=FALSE, as.is=TRUE, sep=" ")
evec = rac[,1]
evec[evec == 1] = "HsNo"
evec[evec %in% 2:5] = "Hs"

# 
qlist = lapply(as.list(files_subjects), function(f) read.table(f, header=F))
tmpk = sapply(strsplit(files_subjects, '.', fixed=T), function(f) f[15])
names(qlist) = paste0("k", tmpk)
#x=read.table("R2.1_geno.01_maf.05_hwe.000001.midp_snps.acgt_excl.high.LD.reg_indep-pair.1000.100..03_rsid.4.Q", header=F)



#
#k = 4
for (k in as.numeric(tmpk)) {
  cat(k,"...")
  x = qlist[[which(tmpk==as.character(k))]]
  for (ordby in 1:k) {
    png(file=paste0("k",k ,"_ordered_by_col.",ordby,"_withHs.png") ,width=2200,height=4000)
    par(mfrow=c(10,3))
    xcex = 2
    rtoplot = c("W","Ot","AI","AA","PI","AsCh","AsJp","AsFp","AsIn","AsOt")
    #rtoplot = c("W","AsCh","Ot","AsJp","AI","AsFp","AA","AsIn","PI","AsOt")
    for (r in rtoplot) {
      for(eth in c("","HsNo","Hs")) {
        if (eth=="") {
          toplot = x[rac[,2]==r,]
          xmain = paste0(r,", n=",nrow(toplot))
        } else {
          toplot = x[evec==eth & rac[,2]==r,]
          xmain = paste0(r,":",eth,", n=",nrow(toplot))
        }
        barplot(t(as.matrix(toplot[order(toplot[,ordby]),])),
                border=NA, col=rainbow(k),
                xlab='',main=xmain,
                names.arg=rep('',nrow(toplot)),
                cex.axis=xcex,cex.main=xcex,
                space=0
        )
        xat = quantile(1:nrow(toplot), seq(0,1,.1))
        axis(side=1, at=xat, labels=names(xat), cex.axis=xcex)
        abline(v=xat[c(3,6,9)],lty='dashed')
        
        
      }
    }
    dev.off()
  }
}

qmeans = list()
for (qq in names(qlist)) {
  x = qlist[[qq]]
  tmp = apply(x,2,function(f) tapply(f, rac[,2], mean))
  qmeans[[qq]] = tmp[grep(":",rownames(tmp),invert=T),]
}

qmax = list()
for (qq in names(qlist)) {
  x = qlist[[qq]]
  tmp = apply(x,2,function(f) tapply(f, rac[,2], max))
  qmax[[qq]] = tmp[grep(":",rownames(tmp),invert=T),]
}

qmedians = list()
for (qq in names(qlist)) {
  x = qlist[[qq]]
  tmp = apply(x,2,function(f) tapply(f, rac[,2], median))
  qmedians[[qq]] = tmp[grep(":",rownames(tmp),invert=T),]
}

colwmaxlist = list()
for (k in names(qlist)) {
  colwmax = list()
  for (r in rtoplot) {
    colwmax[[r]] = table(apply(qlist[[k]][rac[,2]==r,], 1, 
                               function(f) which(f==max(f))))
  }
  colwmaxlist[[k]] = colwmax
}
colwmaxlist2 = list()
colwmaxlistnums = as.numeric(gsub("k", "", names(colwmaxlist)))
for (j in 1:length(colwmaxlist)) {
  tmp = colwmaxlist[[j]]
  tmpnum = colwmaxlistnums[j]
  for (i in names(tmp)) {
    if (length(tmp[[i]]) != tmpnum) {
      tmp[[i]] = tmp[[i]][match(1:tmpnum, names(tmp[[i]]))]
      tmp[[i]][is.na(names(tmp[[i]]))] = 0
      names(tmp[[i]]) = 1:tmpnum
    }
  }
  colwmaxlist2[[j]] = do.call(rbind,tmp)
}
names(colwmaxlist2) = names(colwmaxlist)

for (xk in names(colwmaxlist2)) {
  towrite = cbind(signif(qmeans[[xk]], 4), 
                  rep("", nrow(qmeans[[xk]])), 
                  colwmaxlist2[[xk]][(match(rownames(qmeans[[xk]]), 
                                            rownames(colwmaxlist2[[xk]]))),])
  towrite = towrite[rownames(towrite)!="Missing", ]
  write.table(towrite,file=paste0(xk,"_means_nummax",".txt"),col.names=F,row.names=T,quote=F,sep="\t")
}




qmaxmeans = sapply(qmeans, function(f) apply(f, 1, max))
qmaxmeans = qmaxmeans[,c(3:8,1:2)]

qmaxmeancols = sapply(qmeans, function(ff) apply(ff, 1, function(f) which(f==max(f))))
qmaxmeancols = qmaxmeancols[, c(3:8,1:2)]

xo = "order("
for (i in 1:ncol(qmaxmeancols)) {
  xo = paste0(xo, "qmaxmeancols[,",i,"], ")
}
xo = gsub(", $", ")", xo)
qmaxmeancols = qmaxmeancols[eval(parse(text=xo)), ]
qmaxmeancols0 = qmaxmeancols

xknums = as.numeric(gsub("k", "", colnames(qmaxmeancols)))
for (i in 1:ncol(qmaxmeancols)) {
  rcols = rainbow(xknums[i]);#print(rcols)
  for (j in 1:length(rcols)) {
    qmaxmeancols[qmaxmeancols[,i]==j, i] = rcols[j]
  }
}
qmaxmeans2 = signif(qmaxmeans, 2)
qmaxmeans2 = qmaxmeans2[match(rownames(qmaxmeancols0), rownames(qmaxmeans2)), ]
tmp=paste0(qmaxmeancols0,"\n",qmaxmeans2)
dim(tmp) = dim(qmaxmeancols0)

library(plotrix)
color2D.matplot(qmaxmeancols0, cellcolors=qmaxmeancols, axes=FALSE, xlab="", ylab="")
axis(1, at=(1:8)-.5, labels=colnames(qmaxmeancols0))
axis(2, at=(11:1)-.5, labels=rownames(qmaxmeancols0), las=2)

for (xcol0 in 1:ncol(tmp)) {
  for (xrow0 in 1:nrow(tmp)) {
    xrow = nrow(tmp) - xrow0 + .5
    xcol = xcol0 - .5
    text(xcol, xrow, as.character(tmp[xrow0, xcol0]))
  }
}




dev.off()
ylim = c(as.numeric(substring(min(qmaxmeans),1,4)), 1)
lcol = rainbow(nrow(qmaxmeans))
ltype = 'o'
pch = 19
xk = as.numeric(gsub("k", "", colnames(qmaxmeans)))
plot(xk, qmaxmeans[1, ], col=lcol[1], ylim=ylim, type=ltype, pch=pch)
for (i in 2:nrow(qmaxmeans)) {
  par(new=TRUE)
  plot(xk, qmaxmeans[i, ], col=lcol[i], ylim=ylim, type=ltype, pch=pch)
}



