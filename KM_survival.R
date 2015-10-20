# Survival Curves - Oncogene Addiction
# adapted from http://rstudio-pubs-static.s3.amazonaws.com/5588_72eb65bfbe0a4cb7b655d2eee0751584.html

load(file="dat_compact.RData") # includes oncogene addiction prediction on 9 genes
dat_compact <- dat_compact[!is.na(dat_compact$OncogeneAddiction_Pred),]
head(dat_compact)

library(survival)

## Add survival object
dat_compact$SurvObj <- with(dat_compact, Surv(as.numeric(t.dmfs), e.dmfs == 1))

## Kaplan-Meier estimator. The "log-log" confidence interval.
km.as.one <- survfit(SurvObj ~ 1, data = dat_compact, conf.type = "log-log")
km.by.oa <- survfit(SurvObj ~ OncogeneAddiction_Pred, data = dat_compact, conf.type = "log-log")

## Plot
plot(km.as.one, mark.time=F)
plot(km.by.oa, mark.time=T,
     ylab = "Distant Metastasis Free",
     xlab='Time [Months]',
     col=c("red", "blue"))
legend("topright",
       legend=c("Oncogene Addicted", "Non-Addicted"),
       fill=c("blue", "red"),
       border=FALSE, bty="n", 
       #y.intersp = 0.7, 
       cex=0.8)

## Load rms package
library(rms)
objNpsurv <- npsurv(formula = Surv(as.numeric(t.dmfs),e.dmfs == 1) ~ OncogeneAddiction_Pred, data = dat_compact)

survplot(objNpsurv)

survplot(fit  = objNpsurv,
         conf = c("none","bands","bars")[1],
         xlab = "", ylab = "Distant Metastasis Free Probability",
         #ylim(0.4,1),
         label.curves = TRUE,                     # label curves directly
         ## label.curves = list(keys = "lines"),  # legend instead of direct label
         levels.only  = FALSE,                    # show only levels, no label
         abbrev.label = FALSE,                    # if label used, abbreviate
         ## fun = function(x) {1 - x},            # Cumulative probability plot         
         loglog   = FALSE,                        # log(-log Survival) plot
         logt     = FALSE,                        # log time
        # time.inc = 100,                          # time increment
        # dots     = TRUE,                         # dot grid
         n.risk   = TRUE,                         # show number at risk
         ## srt.n.risk = 0, sep.n.risk = 0.056, adj.n.risk = 1,
         ## y.n.risk = 0, cex.n.risk = 0.6
)


## Plot cumulative probability F(t) = 1 - S(t)
survplot(fit  = objNpsurv,
         conf = c("none","bands","bars")[1],
         xlab = "", ylab = "Cumulative Incidence",
         ## xlim(0,100),
         label.curves = TRUE,                     # label curves directly
         ## label.curves = list(keys = "lines"),  # legend instead of direct label
         levels.only  = FALSE,                    # show only levels, no label
         abbrev.label = FALSE,                    # if label used, abbreviate
         fun = function(x) {1 - x},             # Cumulative probability plot         
         loglog   = FALSE,                        # log(-log Survival) plot
         logt     = FALSE,                        # log time
         #time.inc = 100,                          # time increment
         dots     = FALSE,                        # dot grid
         n.risk   = TRUE,                         # show number at risk
         ## srt.n.risk = 0, sep.n.risk = 0.056, adj.n.risk = 1,
         ## y.n.risk = 0, cex.n.risk = 0.6
)

