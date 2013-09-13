# Thesis R-code
#This version last checked 12-9-2013

#---- Preliminary
#setwd("~/Documents/thesis/draft")
rm(list = ls())  # clear workspace

library(tseries)
library(lattice)
library(lmtest)
library(zoo)
library(vars)
library(CDVine)
library(aod)
library(dynlm)
library(mvtsplot)
library(gridExtra)
library(plyr)


#------ Functions
subdat = function(x, y){  # subset data by fx
  out = subset(x, fx==y, c(date, price, net, dealer_net, 
                       asset_net, lev_net, other_net))
  out = out[343:1,]  # rearrange time-series to start at earliest date
  return(out)
}


# Read data
dat=read.csv("dat.csv")
dat$date = as.Date(dat$date)


# create data subsets
aud=subdat(dat, "AUD")
head(aud)
cad=subdat(dat, "CAD")
chf=subdat(dat, "CHF")
eur=subdat(dat, "EUR")
gbp=subdat(dat, "GBP")
jpy=subdat(dat, "JPY")

# create combinatorial matrix (list of pairwise combinations of interest)

combi = function(x){  # create combination matrix
  cb = as.vector(colnames(x[,2:7])) # get vector as column names
  t(combn(cb, 2)) # get transpose of matrix of unique pair combinations
}
tr=combi(aud)
tr = tr[-6:-9,]  # create combination matrix (remove combinations of net and its composites)



# 1 ----- Data and Summary Statistics

# functions
summ = function(dat) {  # summary stats function
  nam = names(dat);
  cat(sprintf("%11s  %5s %8s %8s %8s %8s %8s %8s", "name", "obs", "mean", "sd", "min", "max", "skew", "kurt"), "\n");
  for (j in seq.int(ncol(dat))) {
    x = dat[is.finite(dat[,j]),j];
    n = length(x);
    if (n) {
      minmax = range(x);
      if (class(x)=="Date") {
        cat(sprintf("%11s: %5i %8s %8s %8s %8s %8s %8s", nam[j], n, "", "", minmax[1], minmax[2], "", ""), "\n");
      } else {
        mu = mean(x);
        sig = sd(x);
        z = (x - mu)/sig;
        skew = mean(z^3)*n/(n-1);
        kurt = mean(z^4)*n/(n-1);
        cat(sprintf("%11s: %5i %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f", nam[j], n, mu, sig, minmax[1], minmax[2], skew, kurt), "\n");
      }
    } else {
      cat(sprintf("%11s: %5i (%s)", nam[j], sum(!is.na(dat[,j])), class(dat[,j])), "\n");
    }
  }
}

# multivariate time series plot *INVESTOR HETEROGENEITY* - for Australian dollar illustration
#pdf('~/Documents/thesis/draft/heterogeneity.pdf', h=6, w=6)
mvtsplot(aud[,c(4:7)], norm = "global", smooth.df = 200, # spline 100 df 
         levels = 21, main = "", margin=T)
#dev.off()


# Summary stats
summ(aud)
summ(cad)
summ(chf)
summ(eur)
summ(gbp)
summ(jpy)

# 2 ----- Unit Root Testing

#Visualise Data - all fx
xyplot(price ~ date | fx, data = dat, typ = "l", scales = "free", as.table = TRUE)  # price
xyplot(scale(net) + scale(lev_net) + scale(dealer_net) + scale(other_net) + scale(asset_net) ~ date | fx, data = dat, typ = "l", # net positions
       as.table = TRUE, ylab= "", xlab = "", scales = list(y= list(relation = "free")),
       key=simpleKey(text=c("net", "leveraged funds", "dealers", "other large traders", "asset managers"), 
                     lines=TRUE, points = FALSE))

# Dickey- Fuller regressions
d.f = function(x,y){
  options(digits=4, scipen=5)
  l.x = c(NA, y[[x]][1:342])
  l2.x = c(NA, l.x)
  l3.x = c(NA, l2.x)
  l4.x = c(NA, l3.x)
  l5.x = c(NA, l4.x)
  l6.x = c(NA, l5.x) 
  x=diff(y[[x]])
  l.x = l.x[2:343]
  t = 1:342 # time dummy
  
  mod_1 = lm(x ~ t + l.x, na.action = na.exclude)  # dickey-fuller regression  
  mod_2 = lm(x[2:342] ~ t[2:342] + l.x[2:342] + diff(l.x)) # augmented with 1 lag
  mod_3 = lm(x[2:342] ~ t[2:342] + l.x[2:342] + diff(l.x) + diff(l2.x[2:343]))   # 2 lags
  mod_4 = lm(x[2:342] ~ t[2:342] + l.x[2:342] + diff(l.x) + diff(l2.x[2:343]) + diff(l3.x[2:343]))   #3 lags
  mod_5 = lm(x[2:342] ~ t[2:342] + l.x[2:342] + diff(l.x) + diff(l2.x[2:343]) + diff(l3.x[2:343]) + diff(l4.x[2:343]))     #4 lags
  mod_6 = lm(x[2:342] ~ t[2:342] + l.x[2:342] + diff(l.x) + diff(l2.x[2:343]) + diff(l3.x[2:343]) + diff(l4.x[2:343]) + diff(l5.x[2:343]))    #5 lags
  mod_7 = lm(x[2:342] ~ t[2:342] + l.x[2:342] + diff(l.x) + diff(l2.x[2:343]) + diff(l3.x[2:343]) + diff(l4.x[2:343]) + diff(l5.x[2:343]) + diff(l6.x[2:343]))  #6 lags
  
  dfval = function(x){ # extract d-f stat values from model
    summary(x)$coef[,"t value"][3]
  }  
  coval = function(x,y){  # extract reg coefficient values from model
    est =  summary(x)$coef[,"Estimate"][y]
    pval = summary(x)$coef[,"Pr(>|t|)"][y]
    return(c(est, pval))
  }
  
  #Print results
  #cat("lag","DF","DF5%","DF1%","cons","trend", "lag", "d.lag", "d.lag2", "d.lag3", "d.lag4", "dlag5", "dlag6", "\n")
  cat(coval(mod_1, 1), coval(mod_1, 2), coval(mod_1, 3),  dfval(mod_1), "\n" )
  #cat("1  ", dfval(mod_2), -3.42,-3.98,coval(mod_2, 1), coval(mod_2, 2), coval(mod_2, 3), coval(mod_2, 4), "\n" )
  # cat("2  ", dfval(mod_3), -3.42,-3.98,coval(mod_3, 1), coval(mod_3, 2), coval(mod_3, 3), coval(mod_3, 4),coval(mod_3, 5),"\n" )
  #cat("3  ", dfval(mod_4), -3.42, -3.98,coval(mod_4, 1), coval(mod_4, 2), coval(mod_4, 3), coval(mod_4, 4),coval(mod_4, 5),coval(mod_4, 6),"\n" )
  #cat("4  ", dfval(mod_5), -3.42,-3.98,coval(mod_5, 1), coval(mod_5, 2), coval(mod_5, 3), coval(mod_5, 4),coval(mod_5, 5),coval(mod_5, 6),coval(mod_5, 7),"\n" )
  # cat("5  ", dfval(mod_6), -3.42,-3.98,coval(mod_6, 1), coval(mod_6, 2), coval(mod_6, 3), coval(mod_6, 4),coval(mod_6, 5),coval(mod_6, 6),coval(mod_6, 7),coval(mod_6, 8),"\n" )
  cat(coval(mod_7, 1), coval(mod_7, 2), coval(mod_7, 3), coval(mod_7, 4),coval(mod_7, 5),coval(mod_7, 6),coval(mod_7, 7),coval(mod_7, 8),coval(mod_7, 9), dfval(mod_7), "\n" )
} 

#Phillips-Perron Test
PP_test= function(x){
  require(tseries)   
  lst=list(x$price, x$net, x$dealer_net,x$asset_net, x$lev_net,x$other_net)
  sapply(lst, pp.test, type = "Z(t_alpha)")[c(1,4),]
}

#Augmented Dickey-Fuller test statistics
ADF= function(x, k){
  lst=list(x$price, x$net, x$dealer_net,x$asset_net, x$lev_net,x$other_net)
  return(sapply(lst, adf.test, k = k, alternative = "stationary")[c(1,4),])  # return DF stats and pvalues
}
curr = list(aud, cad, chf, eur, gbp, jpy)
sapply(curr, ADF, k = 1)
sapply(curr, ADF, k = 5)

# price
cat("lag","DF","DF5%","DF1%","cons","trend", "lag", "d.lag", "d.lag2", "d.lag3", "d.lag4", "dlag5", "dlag6", "\n")
sapply(curr, d.f, x = "price")
# net
sapply(curr, d.f, x = "net")
#dealer_net
sapply(curr, d.f, x = "dealer_net")
# asset
sapply(curr, d.f, x = "asset_net")
#lev
sapply(curr, d.f, x = "lev_net")
#other
sapply(curr, d.f, x = "other_net")


sapply(curr, PP_test)  # Phillips Perron


# 3 --- Correlation

# Functions

pk = function(x) {  # Pearson's k 
  pears = function(y){ 
    t=tr[y,1] #take value from first column
    u=tr[y,2] # take value from second column
    p = cor.test(x[[t]],x[[u]], exact=FALSE) 
    #return(p)
    return(c(tr[y,1], tr[y,2], p[c(4,3)]))
    
  }
  ref=1:nrow(tr) # repeat over no of columns of combination matrix
  sapply(ref, FUN = pears) # apply this function over all fifteen rows of combination matrix
}
rho = function(x){   # Spearman's Rho 
  spear = function(y) {  # rank order correlation of variables (data not nid?)
    t=tr[y,1] #take value from first column
    u=tr[y,2] # take value from second column
    p = cor.test(x[[t]],x[[u]], method = "spearman", exact=FALSE) 
    #return(p)
    return(c(tr[y,1], tr[y,2], p[c(4,3)]))    
  }
  ref=1:nrow(tr) # repeat over all rows of tr matrix
  sapply(ref, FUN = spear) # apply this function over all fifteen rows of combination matrix
}

ols_bro = function(x, scale){  # OLS on standardised coefficients 
  
  # firstly standarise coefficients  
  standardise = function(a){
    std = a / sd(a)
    return(std)
    
  }
  if(scale == TRUE){x = as.data.frame(apply(x[,2:7], 2, standardise))}
  else{
    if(scale == FALSE){x=x}
  }
    
  justrunols = function(y){
    t=tr[y,1] #take value from first column
    u=tr[y,2] # take value from second column
    reg = lm(x[[t]] ~ x[[u]])
    
    if(scale == TRUE){
      out = c(tr[y,1], tr[y,2],#summary(reg)$coefficients[c(2,8)],  # estimate and p-val (s.e. = 4)
              coeftest(reg, vcov = vcovHAC)[c(2, 8)],  # robust p-value
              summary(reg)$r.squared, dwtest(reg)[c(1,4)])
    }
    else{
      if(scale == FALSE){
        out = c(tr[y,1], tr[y,2],#summary(reg)$coefficients[c(2,8)],  # estimate and p-val (s.e. = 4)
                coeftest(reg, vcov = vcovHAC)[2]*10^7)  # estimate scaled up by 10^7 (all other information is the same)
      }
    }
    return(out)
           
  } 
  ref=1:nrow(tr) #  apply across all rows of combination matrix tr
  if(scale == TRUE){
  sapply(ref, FUN = justrunols)
  } else {   if(scale == FALSE){sapply(ref[1:5], FUN = justrunols)}     }
}



ccf_viz = function(x,y){
  bigt= deparse(substitute(x))
  x1=as.data.frame(x[2:343 ,c(1,3:7)])
  x1 = data.frame(x1[,1], diff(x[,2]), x1[,c(2:6)])
  x=x1
  header=c("date", "diff_price", "net", "dealer_net","asset_net", "lev_net","other_net")
  colnames(x)=header    
  #Open graphics devices
  pdf(y, width = 12, height = 6)
  par(mfrow = c(2,6)) # 6X6 
 
  combi = function(x){  # create combination matrix
    cb = as.vector(colnames(x[,2:7])) # get vector as column names
    t(combn(cb, 2)) # get transpose of matrix of unique pair combinations
  }
  tr=combi(x)
  tr = tr[-6:-9,]  # create combination matrix
  #cat(deparse(substitute(x)))
  
  viz = function(y, z){ 
    t=tr[y,1] #take value from first column
    u=tr[y,2] # take value from second column
    #bind=cbind(x[[t]], x[[u]]) # bind two vectors from data based on above colnames 
    ccf(x[[t]], x[[u]], main = paste(bigt,':', u,"&",t))
    #title(bigt, outer=TRUE)
  }
  ref=1:11 # list of 1-11 (there are 11 combinations excluding net and dealer/asset/lev/other)
  mapply(FUN = viz, ref, title) # apply this function over all fifteen rows of combination matrix
}


# Matrix Scatterplot: Price
p1 = xyplot(net  ~ price | fx, data = dat, scales = list(x= list(relation = "free")), as.table=T, layout = c(6,1),
      panel = function(x,y,...){
         panel.xyplot(x, y, cex = 0.5)
         panel.lmline(x,y,col = "red")
       })
p2 =  xyplot(dealer_net  ~ price | fx, data = dat, scales = list(x= list(relation = "free")),as.table=T, layout = c(6,1),
       panel = function(x,y,...){
         panel.xyplot(x, y, cex = 0.5)
         panel.lmline(x,y,col = "red")
       })
p3 = xyplot(lev_net  ~ price | fx, data = dat, scales = list(x= list(relation = "free")),as.table=T, layout = c(6,1),
       panel = function(x,y,...){
         panel.xyplot(x, y, cex = 0.5)
         panel.lmline(x,y,col = "red")
       })
p4 = xyplot(asset_net  ~ price | fx, data = dat, scales = list(x= list(relation = "free")),as.table=T, layout = c(6,1),
       panel = function(x,y,...){
         panel.xyplot(x, y, cex = 0.5)
         panel.lmline(x,y,col = "red")
       })
p5 = xyplot(other_net  ~ price | fx, data = dat, scales = list(x= list(relation = "free")),as.table=T, layout = c(6,1),
       panel = function(x,y,...){
         panel.xyplot(x, y, cex = 0.5)
         panel.lmline(x,y,col = "red")
       })
#grid.arrange(p1, p2, p3, p4, p5, ncol = 1)

# Matrix Scatterplot: Positions

p6 = xyplot(asset_net  ~ dealer_net | fx, data = dat,  scales = list(x= list(relation = "free")),as.table=T, layout = c(6,1), 
            panel = function(x,y,...){
               panel.xyplot(x, y, cex = 0.5)
              panel.lmline(x,y,col = "red")
            })
p7 = xyplot(lev_net  ~ dealer_net | fx, data = dat, scales = list(x= list(relation = "free")),as.table=T, layout = c(6,1),
            panel = function(x,y,...){
              panel.xyplot(x, y, cex = 0.5)
              panel.lmline(x,y,col = "red")
            })

p8 = xyplot(other_net  ~  dealer_net | fx, data = dat, scales = list(x= list(relation = "free")),as.table=T, layout = c(6,1),
            panel = function(x,y,...){
              panel.xyplot(x, y, cex = 0.5)
              panel.lmline(x,y,col = "red")
            })

p9 = xyplot(lev_net ~ asset_net | fx, data = dat, scales = list(x= list(relation = "free")),as.table=T, layout = c(6,1),
            panel = function(x,y,...){
              panel.xyplot(x, y, cex = 0.5)
              panel.lmline(x,y,col = "red")
            })
p10 = xyplot(other_net  ~ asset_net | fx, data = dat, scales = list(x= list(relation = "free")),as.table=T, layout = c(6,1),
            panel = function(x,y,...){
              panel.xyplot(x, y, cex = 0.5)
              panel.lmline(x,y,col = "red")
            })
p11 = xyplot(other_net  ~ lev_net | fx, data = dat, scales = list(x= list(relation = "free")),as.table=T, layout = c(6,1),
            panel = function(x,y,...){
              panel.xyplot(x, y, cex = 0.5)
              panel.lmline(x,y,col = "red")
            })
#grid.arrange(p6, p7, p8, p9, p10, p11, ncol = 1)

# Pearson's k 

pk(aud)
pk(cad)
pk(chf)
pk(eur)
pk(gbp)
pk(jpy)

# Spearmnan's \rho
rho(aud)
rho(cad)
rho(chf)
rho(eur)
rho(gbp)
rho(jpy)



# OLS estimates

curr = list(aud, cad, chf, eur, gbp, jpy)
sapply(curr, ols_bro, scale = TRUE)
sapply(curr, ols_bro, scale = FALSE)  # note jpy scale is different


#ols_bro(aud)
#ols_bro(cad)
#ols_bro(chf)
#ols_bro(eur)
#ols_bro(gbp)
#ols_bro(jpy)



# CCF for AUD
#ccf_viz(aud, "~/Documents/thesis/draft/ccf.pdf"); dev.off()
ccf_viz(cad, "~/Documents/thesis/draft/ccf_cad.pdf"); dev.off()
ccf_viz(chf, "~/Documents/thesis/draft/ccf_chf.pdf"); dev.off()
ccf_viz(eur, "~/Documents/thesis/draft/ccf_eur.pdf"); dev.off()
ccf_viz(gbp, "~/Documents/thesis/draft/ccf_gbp.pdf"); dev.off()
ccf_viz(jpy, "~/Documents/thesis/draft/ccf_jpy.pdf"); dev.off()



# 4 ----- Multivariate Methods 
# VAR 
vars = function(x, s, z) { 
    vec_ar = function(y){ 
    t=tr[y,1] #take value from first column
    u=tr[y,2] # take value from second column
    
    vec1 = x[[t]]
    vec2 = x[[u]]
    bind=cbind(vec1, vec2) # bind two vectors from data based on above colnames 
    
    criteria=VARselect(bind, lag.max=10)  # select lag for var 
    sel = criteria$selection[c(1,3)] 
    if(s == "max"){pval = as.integer(max(sel))
    } else {
      if(s == "min"){pval = as.integer(min(sel))}
    } 
    
    est = VAR(bind, p = pval )
    granger_test1 = causality(est, cause = "vec1")
    granger_test2 = causality(est, cause = "vec2")
    
    if(z == "granger"){
    return(c(granger_test1$Granger[c("statistic", "p.value")], "lags" = pval, 
               granger_test2$Granger[c("statistic", "p.value")], "lags" = pval))
    } else {
      if(z == 'instant'){
        return(c(granger_test1$Instant[c(1,3,2)]))
                   
      }
    }
    }
  ref=1:nrow(tr) #  apply across all rows of combination matrix tr
  sapply(ref, FUN =vec_ar) # apply this function over all fifteen rows of combination matrix
}
curr = list(aud, cad, chf, eur, gbp, jpy)
sapply(curr, vars, s = "max", z = 'granger')
sapply(curr, vars, s = "min", z = 'granger')
sapply(curr, vars, s = "min", z = 'instant')


# Todo - Yamamoto Correction 
TY.method = function(x,s){
  toda.yam = function(y){ 
    t=tr[y,1] #take value from first column
    u=tr[y,2] # take value from second column
    
    vec1 = zoo(x[[t]])
    vec2 = zoo(x[[u]])
    
    m=1 # order of integration of price    
    
    bind=cbind(vec1, vec2) # bind two vectors from data based on above colnames 
    
    criteria=VARselect(bind, lag.max=10)  # select lag for var
    sel = criteria$selection[c(1,3)] 
    if(s == "max"){pval = as.integer(max(sel))
    } else {
      if(s == "min"){pval = as.integer(min(sel))}
    } 
    
    vec1.l = lag(vec1,-(0:(pval+1)),na.pad=T)  # create lags for pval plus addional TY parameter
    vec2.l = lag(vec2,-(0:(pval+1)),na.pad=T)
    
    lm1=lm(vec1~vec2.l[,2:pval+1]+vec1.l[,2:pval+1])
    lm2=lm(vec2~vec1.l[,2:pval+1]+vec2.l[,2:pval+1])
    
    l=pval+1
    #print(pval)
    #print(l)
    
    modx=dynlm(bind[,1]~L(bind[,2],1:l)+L(bind[,1],1:l))
    mody=dynlm(bind[,2]~L(bind[,1],1:l)+L(bind[,2],1:l))
    
    w1=wald.test(b=coef(modx), Sigma=vcov(modx), Terms= c(2:l))
    w2=wald.test(b=coef(mody), Sigma=vcov(mody), Terms= c(2:l))
    
    res_w1 = w1$result$chi2
    res_w2 = w2$result$chi2       
    return(c(l,res_w2, "",
             res_w1))
  }
  ref=1:5 
  sapply(ref, FUN =toda.yam) 
}
sapply(curr, TY.method, s = "max")
sapply(curr, TY.method, s = "min")


# IRFs 
select_p = function(x){
    irf_var = function(y,z){
    t=tr[y,1]
    u=tr[y,2]
    bind=cbind(x[[t]], x[[u]])
    criteria=VARselect(bind, lag.max=10)
    lag_select = min(criteria$selection[c(1,3)])
    return(lag_select)
   }
  ref=1:nrow(tr)
  sapply(ref, FUN=irf_var)
}
min_p = select_p(aud)

criteria = cbind(tr, min_p)
print(criteria)


imp = function(x){
  IRF = function(fx){   
        vec = cbind(fx[[criteria[x,1]]], fx[[criteria[x,2]]])
        colnames(vec) = c(criteria[x,1], criteria[x,2])  
        imp1 = irf(VAR(vec,  p=2 ), response = criteria[x,1], impulse = criteria[x,2])
        imp2 = irf(VAR(vec,  p=2 ), response = criteria[x,2], impulse = criteria[x,1])
        out = data.frame(imp1$irf[[criteria[x,2]]] , imp1$Upper[[criteria[x,2]]] , imp1$Lower[[criteria[x,2]]], 
                         imp2$irf[[criteria[x,1]]] , imp2$Upper[[criteria[x,1]]] , imp2$Lower[[criteria[x,1]]],                        
                   fx = (deparse(substitute(fx))), ref = 0:10)
        return(out)
        }
  m1 = IRF(aud)
  m2 = IRF(cad)
  m3 = IRF(chf)
  m4 = IRF(eur)
  m5 = IRF(gbp)
  m6 = IRF(jpy)
  bind = rbind(m1, m2, m3, m4, m5, m6)
  
  
  plot_out1 = xyplot(bind[,1] + bind[,2] + bind[,3] ~ ref | fx, data = bind, typ="l", ylab = criteria[x,1],  scales = list(y= list(relation = "free")),
       xlab = paste("Impulse Response from", criteria[x,2]), col = c("blue","red", "red"), lty = c(1,2,2), layout = c(6,1)
                     ,
         panel = function(...){
         panel.abline(a=0, b=0, col = gray(0.5))
         panel.xyplot(...)}
                     )
  plot_out2 = xyplot(bind[,4] + bind[,5] + bind[,6] ~ ref | fx, data = bind, typ="l", ylab = criteria[x,2],  scales = list(y= list(relation = "free")),
                   xlab = paste("Impulse Response from", criteria[x,1]), col = c("blue","red", "red"), lty = c(1,2,2), layout = c(6,1)
                   ,
                   panel = function(...){
                     panel.abline(a=0, b=0, col = gray(0.5))
                     panel.xyplot(...)}
                     )
  return(c(plot(plot_out1), plot(plot_out2)))
}
z = 1:11
pdf('~/Documents/thesis/draft/irf_test.pdf', h=4, w=12)
sapply(z, FUN = imp)
dev.off()

# Subsample robustness


sub_vars = function(x,y,z, type ) { 
  n=3
  dfchunk = split(x, factor(sort(rank(row.names(x))%%n)))
  out = dfchunk[[y]]
  cat(deparse(substitute(x)), 'subsample:', deparse(substitute(y)), '\n')
  
  range = function(z){
    mn = min(z$date)
    mx = max(z$date)
    cat("Time series from", paste(mn), "to", paste(mx), "\n")
  }
  range(out)
  
  vec_ar = function(y){ 
    t=tr[y,1] #take value from first column
    u=tr[y,2] # take value from second column
    
    vec1 = out[[t]]
    vec2 = out[[u]]
    bind=cbind(vec1, vec2) # bind two vectors from data based on above colnames 
    
    criteria=VARselect(bind, lag.max=10)  # select lag for var 
    pval = criteria$selection[1] # use AIC value
       
    est = VAR(bind, p = pval )
    granger_test1 = causality(est, cause = "vec1")
    granger_test2 = causality(est, cause = "vec2")
    
    if(type == 'granger'){
      return(c(granger_test1$Granger[c("statistic", "p.value")], "lags" = pval, 
               granger_test2$Granger[c("statistic", "p.value")], "lags" = pval))
    } else {
      if(type  == 'instant'){
        return(c(granger_test1$Instant[c(1,3,2)]))
        
      }
    }
  }
  toda.yam = function(y){ 
    t=tr[y,1] #take value from first column
    u=tr[y,2] # take value from second column
    
    vec1 = zoo(out[[t]])
    vec2 = zoo(out[[u]])
    
    m=1 # order of integration of price    
    
    bind=cbind(vec1, vec2) # bind two vectors from data based on above colnames 
    
    criteria=VARselect(bind, lag.max=10)  # select lag for var
    pval = criteria$selection[1]     
    
    vec1.l = lag(vec1,-(0:(pval+1)),na.pad=T)  # create lags for pval plus addional TY parameter
    vec2.l = lag(vec2,-(0:(pval+1)),na.pad=T)
    
    lm1=lm(vec1~vec2.l[,2:pval+1]+vec1.l[,2:pval+1])
    lm2=lm(vec2~vec1.l[,2:pval+1]+vec2.l[,2:pval+1])
    
    l=pval+1
    #print(pval)
    #print(l)
    
    modx=dynlm(bind[,1]~L(bind[,2],1:l)+L(bind[,1],1:l))
    mody=dynlm(bind[,2]~L(bind[,1],1:l)+L(bind[,2],1:l))
    
    w1=wald.test(b=coef(modx), Sigma=vcov(modx), Terms= c(2:l))
    w2=wald.test(b=coef(mody), Sigma=vcov(mody), Terms= c(2:l))
    
    res_w1 = w1$result$chi2[c(1,3,2)]
    res_w2 = w2$result$chi2[c(1,3,2)]    
    return(c(res_w2,
             res_w1))
  }
  ref=1:nrow(tr) #  apply across all rows of combination matrix tr
  if(z == FALSE){sapply(ref, FUN =vec_ar) 
  } else {
    if(z==TRUE){sapply(ref[1:5], FUN =toda.yam)
    }
  }
}

curr = list(aud, cad, chf, eur, gbp, jpy)

sapply(curr, FUN = sub_vars, 1, z = FALSE, type = "granger")  # Do for standard Granger-causality
sapply(curr, FUN = sub_vars, 2, z = FALSE, type = "granger") 
sapply(curr, FUN = sub_vars, 3, z = FALSE, type = "granger") 

sapply(curr, FUN = sub_vars, 1, z = FALSE, type = "instant")  # Do for instantaneous-causality
sapply(curr, FUN = sub_vars, 2, z = FALSE, type = "instant")
sapply(curr, FUN = sub_vars, 3, z = FALSE, type = "instant")

sapply(curr, FUN = sub_vars, 1, z = TRUE)  # Do for Toda-Yamamoto method
sapply(curr, FUN = sub_vars, 2, z = TRUE)
sapply(curr, FUN = sub_vars, 3, z = TRUE)



 # 5 --- Copula Stuff

aud=subdat(dat, "AUD")[,c(2,4,6,5,7)] # note we have gotten rid of the 'net' variable here, as we are not interested in a vine copula
# involving the net of all large participants
cad=subdat(dat, "CAD")[,c(2,4,6,5,7)]
chf=subdat(dat, "CHF")[,c(2,4,6,5,7)]
eur=subdat(dat, "EUR")[,c(2,4,6,5,7)]
gbp=subdat(dat, "GBP")[,c(2,4,6,5,7)]
jpy=subdat(dat, "JPY")[,c(2,4,6,5,7)]

par(mfrow = c(1,3))
options("scipen" = 20)
plot(aud$dealer_net, typ = "l", xlab="dealer_net", ylab=NA, main = "Time Series")
plot(ecdf(aud$dealer_net), xlab = "dealer_net", main = "Empirical CDF", ylab = "Fn(dealer_net)")

as.uniform = function(x){
  trans = ecdf(x) # estimate cum dist function 
  out = round(trans(x), 8)
  out[out == 1] = 0.9999  # get rid of any values equal to 1 (data must be on interval [10^(-10), 1 - 10^(-10)])
  # due to stability issues
  return(out) # apply to data
}

aud = apply(aud, 2, as.uniform)  # apply uniform function over cols price and particpants
aud = as.data.frame(aud)  # reformat as data.frame
head(aud)  # we now have price + 4 participants

cad = apply(cad, 2, as.uniform) 
cad = as.data.frame(cad)  
head(cad) 

chf = apply(chf, 2, as.uniform) 
chf = as.data.frame(chf)  
head(chf) 

eur = apply(eur, 2, as.uniform) 
eur = as.data.frame(eur)  
head(eur) 

gbp = apply(gbp, 2, as.uniform) 
gbp = as.data.frame(gbp)  
head(gbp) 

jpy = apply(jpy, 2, as.uniform) 
jpy = as.data.frame(jpy)  
head(jpy) 

apply(aud, 2, max) 
apply(cad, 2, max)
apply(chf, 2, max)
apply(eur, 2, max)
apply(gbp, 2, max)
apply(jpy, 2, max)

apply(aud, 2, min) 
apply(cad, 2, min)
apply(chf, 2, min)
apply(eur, 2, min)
apply(gbp, 2, min)
apply(jpy, 2, min)


### Show transformation in plot
plot(aud$dealer_net, typ = "l", xlab="dealer_net", ylab=NA, main = "Transformed Series")


panel.contour <- function(x, y, bw=2, size=100){
  usr <- par("usr")
  on.exit(par(usr))
  par(usr = c(-3,3,-3,3), new=TRUE)
  BiCopMetaContour(x, y, bw, size, axes=FALSE)
}
pdf('~/Documents/thesis/draft/contour.pdf', h=6, w=6)
pairs(aud, lower.panel=panel.contour, gap=0, main = "Australian Dollar") 
pairs(chf, lower.panel=panel.contour, gap=0, main = "Canadian Dollar") 
pairs(cad, lower.panel=panel.contour, gap=0, main = "Swiss Franc") 
pairs(eur, lower.panel=panel.contour, gap=0, main = "Euro") 
pairs(gbp, lower.panel=panel.contour, gap=0, main = "British Pound") 
pairs(jpy, lower.panel=panel.contour, gap=0, main = "Japanese Yen") 
dev.off()

#Write Copula Algorithm 

copula_algo = function(x){
  
  #x = x[,c(1,3,2,4,5)]  
  family = CDVineCopSelect(x, type = 1, indeptest = F, familyset=c(1:10,13,14,16:20, 23, 24, 
                                                                     26:30, 33, 34, 36:40))$family
  family[family==2] = 1  # get rid of independent cops - replace with Gaussian
  family
  seq_par = CDVineSeqEst(x, type = 1, family = family, method = "mle")
  mle_par = CDVineMLE(x, type = 1, family = family, 
                      start = seq_par$par, start2 = seq_par$par2)
  mle_par$convergence == 0  # check if ML converges
  mle_par
  
  par(mfrow = c(2,2))
  CDVineTreePlot(x, type = 1, par=mle_par$par,  par2=mle_par$par2, tree = 1,
                 family = family, edge.labels=c("family","emptau"))
  CDVineTreePlot(x, type = 1, par=mle_par$par,  par2=mle_par$par2, tree = 2,
               family = family, edge.labels=c("family","emptau"))
  CDVineTreePlot(x, type = 1, par=mle_par$par,  par2=mle_par$par2, tree = 3,
                 family = family, edge.labels=c("family","emptau"))
  CDVineTreePlot(x, type = 1, par=mle_par$par,  par2=mle_par$par2, tree = 4,
                family = family, edge.labels=c("family","emptau"))
  
}
#pdf('~/Documents/thesis/draft/tree.pdf', h=6, w=6)
copula_algo(aud)
#dev.off()
#copula_algo(cad)
#copula_algo(chf)
#copula_algo(eur)
#copula_algo(gbp)
#copula_algo(jpy)

ind_test = function(x){
  cb = as.vector(colnames(aud)) # get vector as column names
  tr= t(combn(cb, 2)) # get transpose of 
  fun = function(y){   
    t=tr[y,1] #take value from first column of comb matrix
    u=tr[y,2] # take value from second column 
    vec1 = x[[t]]
    vec2 = x[[u]]   
    head = paste(t,",",u) # create header for output
    return(c(round(BiCopIndTest(vec1, vec2)$p.value, 3)))
    # print header plus p-values and round to 3 places
  }
  #ref = 2:nrow(tr)
  ref = 1:nrow(tr)
  sapply(ref, fun)
}  
curr= list(aud, cad, chf, eur, gbp, jpy)
sapply(curr, ind_test)



tau_est = function(x){
  family = CDVineCopSelect(x, type = 1, indeptest = F, familyset=c(1:10,13,14,16:20, 23, 24, 
                                                                   26:30, 33, 34, 36:40))$family
  family[family==2] = 1  # get rid of independent cops - replace with Gaussian
  family[family==29] = 30  # get rid of bb7 cops - replace with bb8
  family[family==39] = 40  # get rid of bb7 cops - replace with bb8
  seq_par = CDVineSeqEst(x, type = 1, family = family, method = "mle")
  mle_par = CDVineMLE(x, type = 1, family = family, 
                      start = seq_par$par, start2 = seq_par$par2)
  mle_par$convergence == 0  # check if ML converges
  CDVinePar2Tau(family, par = mle_par$par, par2 = mle_par$par2)
}
curr= list(aud, cad, chf, eur, gbp, jpy)
sapply(curr, tau_est)

#FIN#

  
