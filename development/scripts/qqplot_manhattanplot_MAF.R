args = commandArgs(trailingOnly=TRUE)
summary.file <- args[1]
qqplot <- "./figures/qqplot"
manhattan.plot <- "./figures/manhatten_plot"

####MAF
MAF = 0.01

####input data
dat <- read.csv(summary.file,header = T,stringsAsFactors = F,na.strings=c(NA,""))
dat <- dat[which(!is.na(dat$Score.pval)),]

####qqplot
lambda_logp=function (data, plot = TRUE, proportion = 1, ...) 
{
    data <- data[which(!is.na(data))]
    if (proportion > 1 || proportion <= 0) 
        stop("proportion argument should be greater then zero and less than or equal to one")
    ntp <- round(proportion * length(data))
    if (ntp < 1) 
        stop("no valid measurments")
    if (ntp == 1) {
        warning(paste("One measurment, Lambda = 1 returned"))
        return(list(estimate = 1, se = 999.99))
    }
    if (ntp < 10) 
        warning(paste("number of points is too small:", ntp))
    if (min(data) < 0) 
        stop("data argument has values <0")
    if (max(data) <= 1) {
        data0=data; data <- qchisq(data, 1, low = FALSE)
    }
    data <- sort(data); data0=sort(-log10(data0))
    ppoi <- ppoints(data); ppoi0=sort(-log10(ppoi))
    ppoi <- sort(qchisq(1 - ppoi, 1))
    data <- data[1:ntp]
    ppoi <- ppoi[1:ntp]
    s <- summary(lm(data ~ 0 + ppoi))$coeff
    if (plot) {
        lim <- c(0, max(data0, ppoi0, na.rm = T))

plot(ppoi0, data0, xlab = "Expected (-log10(p))", ylab = "Observed", 
main=paste("Lambda=",signif(s[1,1],5)," SE=",signif(s[1,2],5),sep=""), cex.lab=1.5, cex.main=1.5)
            
        abline(a = 0, b = 1,lwd=2,col="red")
        #abline(a = 0, b = (s[1, 1]), col = "red",lwd=2)
    }
    out <- list()
    out$estimate <- s[1, 1]
    out$se <- s[1, 2]
    out
}

png(paste(qqplot,".png",sep=""))
lambda_logp(dat$Score.pval)
dev.off()

####manhattan plot
data.plot <- c()
for(i in 1:22){
  print(i)
  temp <- dat[which(dat$chr == i),]
  data.plot <- rbind(data.plot,temp[sort(temp$pos,index.return=T)[[2]],])
}

log10.pval <- -log10(data.plot$Score.pval+0.) ## what goes on y-axis
color.length <- c(1,unlist(table(data.plot$chr)))
color.choice <- c(rep(c("darkred", "red"), 11),"darkred")

L.1 <- 0
L.2 <- 0
x.val <- c()
for(i in 1: (length(color.length)-1)){
  L.1 <- L.1 + color.length[i]
  L.2 <- L.2 + color.length[i+1]
  x.val <- c(x.val,(L.1+L.2)/2)
}

lab.chrom <- as.character(c(1:22))
x<- c(1: length(log10.pval))
limit.bf <- min(log10.pval)
L.1 <- 0
L.2 <- 0

png(paste(manhattan.plot,".png",sep=""),width=1440,height=480,pointsize = 11)
par(mai=(c(0.65, 0.35, 0.35, 0.42) ))
plot(x,(log10.pval),ylim=c(limit.bf,max(log10.pval)),ylab="",axes=F,xlab="",cex=1.5,cex.lab=1.5,cex.axis=1.4)
axis(side=1, at=x.val, label=lab.chrom,cex.axis=1.4)
axis(side=2, cex.axis=1.5,line=-2.5)

for(i in 1: (length(color.length)-1)){
  L.1 <- L.1 + color.length[i]
  L.2 <- L.2 + color.length[i+1]
  points( x[L.1:L.2], log10.pval[ L.1 : L.2], col=color.choice[i]) 
}

abline(-log10(5e-8),0,col="grey",lty=1,cex=1.4)
abline(-log10(5e-6),0,col="grey",lty=2,cex=1.4)

dev.off()

####filter out MAF
dat <- dat[dat$freq > MAF & dat$freq < 1-MAF,]

####qqplot (filter out MAF)
png(paste(qqplot,"_MAF_",MAF,".png",sep=""))
lambda_logp(dat$Score.pval)
dev.off()

####manhattan plot (filter out MAF)
data.plot <- c()
for(i in 1:22){
  print(i)
  temp <- dat[which(dat$chr == i),]
  data.plot <- rbind(data.plot,temp[sort(temp$pos,index.return=T)[[2]],])
}

log10.pval <- -log10(data.plot$Score.pval+0.) ## what goes on y-axis
color.length <- c(1,unlist(table(data.plot$chr)))
color.choice <- c(rep(c("darkred", "red"), 11),"darkred")

L.1 <- 0
L.2 <- 0
x.val <- c()
for(i in 1: (length(color.length)-1)){
  L.1 <- L.1 + color.length[i]
  L.2 <- L.2 + color.length[i+1]
  x.val <- c(x.val,(L.1+L.2)/2)
}

lab.chrom <- as.character(c(1:22))
x<- c(1: length(log10.pval))
limit.bf <- min(log10.pval)
L.1 <- 0
L.2 <- 0

png(paste(manhattan.plot,"_MAF_",MAF,".png",sep=""),width=1440,height=480,pointsize = 11)
par(mai=(c(0.65, 0.35, 0.35, 0.42) ))
plot(x,(log10.pval),ylim=c(limit.bf,max(log10.pval)),ylab="",axes=F,xlab="",cex=1.5,cex.lab=1.5,cex.axis=1.4)
axis(side=1, at=x.val, label=lab.chrom,cex.axis=1.4)
axis(side=2, cex.axis=1.5,line=-2.5)

for(i in 1: (length(color.length)-1)){
  L.1 <- L.1 + color.length[i]
  L.2 <- L.2 + color.length[i+1]
  points( x[L.1:L.2], log10.pval[ L.1 : L.2], col=color.choice[i]) 
}

abline(-log10(5e-8),0,col="grey",lty=1,cex=1.4)
abline(-log10(5e-6),0,col="grey",lty=2,cex=1.4)

dev.off()
