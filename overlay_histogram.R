### Plot two overlaying histograms with density function with added color overlay

## Load and clean data

load(file="dat_compact.RData")
dat_compact2 <- na.omit(dat_compact)

# Age distribution based on ER status
er_neg <- dat_compact2[dat_compact2$er==0,]
er_pos <- dat_compact2[dat_compact2$er==1,]

hist(er_neg$age, col=rgb(0.1,0.1,0.1,0.5), probability = T, breaks=30,
xlab='age',
ylim = c(0, 0.06),
main = 'Histogram of age distribution depending on ER status')
abline(v=mean(er_neg$age), col=rgb(0.3,0.3,0.3), lwd=2, lty=2)
lines(density(er_neg$age), col=rgb(0.1,0.1,0.1), lwd=4, lty=4)

hist(er_pos$age, col=rgb(0,0,1,0.5), probability = T, breaks=30, add=T)
abline(v=mean(er_pos$age), col=rgb(0,0,0.5), lwd=2, lty=2)
lines(density(er_pos$age), col=rgb(0,0,0.9), lwd=4, lty=4)

legend("topright",
legend=c("ER Neg", "ER Pos"),
fill=c(rgb(0.1,0.1,0.1,0.5),
rgb(0,0,1,0.5)),
border=FALSE, bty="n",
cex=1.7)
dev.off()

# Age distribution based on OA status
oa_neg <- dat_compact2[dat_compact2$OncogeneAddiction_Pred==0,]
oa_pos <- dat_compact2[dat_compact2$OncogeneAddiction_Pred==1,]

hist(oa_neg$age, col=rgb(0.1,0.1,0.1,0.5), probability = T, breaks=30,
xlab='age',
ylim = c(0, 0.06),
main = 'Histogram of age distribution depending on Oncogene Addiction status')
abline(v=mean(oa_neg$age), col=rgb(0.3,0.3,0.3), lwd=2, lty=2)
lines(density(oa_neg$age), col=rgb(0.1,0.1,0.1), lwd=4, lty=4)

hist(oa_pos$age, col=rgb(0,0,1,0.5), probability = T, breaks=30, add=T)
abline(v=mean(oa_pos$age), col=rgb(0,0,0.5), lwd=2, lty=2)
lines(density(oa_pos$age), col=rgb(0,0,0.9), lwd=4, lty=4)

legend("topright",
legend=c("OA Neg", "OA Pos"),
fill=c(rgb(0.1,0.1,0.1,0.5),
rgb(0,0,1,0.5)),
border=FALSE, bty="n",
cex=1.7)
dev.off()

# T.dmfs distribution based on OA status
oa_neg <- dat_compact2[dat_compact2$OncogeneAddiction_Pred==0,]
oa_pos <- dat_compact2[dat_compact2$OncogeneAddiction_Pred==1,]

hist(oa_neg$t.dmfs, col=rgb(0.1,0.1,0.1,0.5), probability = T, breaks=50,
xlab='Time Distant Metastasis Free Survival [Months]',
ylim = c(0, 0.3),
main = 'Time of Distant Metastasis Free Survival depending on Predicted Oncogene Addiction Status')
abline(v=mean(oa_neg$t.dmfs), col=rgb(0.3,0.3,0.3), lwd=2, lty=2)
lines(density(oa_neg$t.dmfs), col=rgb(0.1,0.1,0.1), lwd=4, lty=4)

hist(oa_pos$t.dmfs, col=rgb(0,0,1,0.5), probability = T, breaks=50, add=T)
abline(v=mean(oa_pos$t.dmfs), col=rgb(0,0,0.5), lwd=2, lty=2)
lines(density(oa_pos$t.dmfs), col=rgb(0,0,0.9), lwd=4, lty=4)

legend("topright",
legend=c("OA Neg", "OA Pos"),
fill=c(rgb(0.1,0.1,0.1,0.5),
rgb(0,0,1,0.5)),
border=FALSE, bty="n",
cex=1.7)
dev.off()

plot(t.dmfs~age, data=oa_neg, pch=20, col='blue')
points(oa_pos$t.dmfs~oa_pos$age, pch=20, col='red')


####
## explore clustering
inds <- head(evar_ordr, length(evar)*.1)
d <- dist(x = t(y), method = 'euclidean')

par(mfrow = c(2, 2))
plot(hclust(d =  d, method = 'ward.D2'))
plot(hclust(d =  d, method = 'average'))
plot(hclust(d =  d, method = 'single'))
plot(hclust(d =  d, method = 'complete'))


########

