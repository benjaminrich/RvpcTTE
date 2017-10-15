#' ---
#' title: "R-Based VPC for Time-To-Event Models with Nonlinear Hazard Functions"
#' author: "Benjamin Rich and Samer Mouksassi"
#' date:   "`r Sys.Date()`"
#' output: html_document
#' ---

#+ echo=FALSE
knitr::opts_chunk$set(fig.path="./", fig.width=4.5, fig.height=4.5, warning=FALSE, message=FALSE)
#+

library(survival)

set.seed(4093)

# Function to randomly sample a time-to-event from a survival function
# specified at discrete time points
# 
# Arguments
# - time        : a numeric vector of increasing times
# - survfn      : a numeric vector giving the value of the survival function at times t
# - interpolate : type of interpolation
#
# Returns       : A data.frame containing the event time and status
sample.tte <- function(time, survfn, interpolate=c("linear", "log-linear", "none")) {
    interpolate <- match.arg(interpolate)
    cdf <- 1 - survfn
    u <- runif(1)
    if (u > max(cdf)) {
        stime <- max(time)
        status <- 0  # Censored
    } else {
        if (interpolate == "linear") {
            fn <- approxfun(cdf, time, method="linear", rule=2)
            stime <- fn(u)
        } else if (interpolate == "log-linear") {
            fn <- approxfun(log(1 - cdf), time, method="linear", rule=2)
            stime <- fn(log(1 - u))
        } else {
            fn <- approxfun(cdf, time, method="constant", rule=2, f=1)
            stime <- fn(u)
        }
        status <- 1  # Event
    }
    data.frame(time=stime, status=status)
}

# Number of simulation repetitions for the VPC
nrep <- 100

# Times at which to compute the survival
tgrid <- seq(0, 150, len=100)

# This matrix will store the survival curves for each repetition
vpcsurv <- matrix(numeric(0), nrep, length(tgrid))


# Data: this is observed data (used to fit both models).
obsdat <- read.csv("SIMDAT100.csv")
obsdat <- obsdat[obsdat$mdv==0,]

#' Note: the true hazard used to generate these data is:
#' $$ h(t) = \exp( -10 - 0.006 t + 1.5 \log(t + 1) ) $$


# Model 1: this is output from fitting the correctly specified model
output1 <- read.csv("table_fit1.csv")

# Model 2: this is output from fitting the miss-specified model
output2 <- read.csv("table_fit2.csv")



#' ### Start of code to generate the VPC figures
#+ vpc-tte

xlab <- "Time (h)"
ylab <- "Probability of Survival"
xlim <- c(0, 150)

x <- survfit(Surv(time, status) ~ 1, data=obsdat)
plot(x, conf.int=FALSE, mark.time=FALSE, lty=0, axes=F, xlim=xlim)
mtext(xlab, 1, 3, las=0)
mtext(ylab, 2, 4, las=0)
axis(1, at=pretty(xlim, 5), las=1)
axis(2, at=seq(0, 1, by=.25), las=1)

for (rep in seq_len(nrep)) {
    vpctte <- by(output1, factor(output1$id),
        function(x) sample.tte(x$time, x$surv))
    vpctte <- do.call(rbind, vpctte)
    x.vpc <- survfit(Surv(time, status) ~ 1, data=vpctte)
    vpcsurv[rep,] <- summary(x.vpc, times=tgrid, extend=T)$surv
}
# Compute median and 90% PI over simulation repetitions
PI <- apply(vpcsurv, 2, quantile, probs=c(0.05, 0.5, 0.95), na.rm=T)
polygon(c(tgrid, rev(tgrid)), c(PI[1,], rev(PI[3,])), border=NA, col="gray90")
lines(tgrid, PI[2,], col="gray70", lwd=3.5)
lines(x, conf.int=F, lwd=2, col="darkcyan")
mtext("CORRECT MODEL", 3, 1, font=2, cex=1.1)


x <- survfit(Surv(time, status) ~ 1, data=obsdat)
plot(x, conf.int=FALSE, mark.time=FALSE, lty=0, axes=F, xlim=xlim)
mtext(xlab, 1, 3, las=0)
mtext(ylab, 2, 4, las=0)
axis(1, at=pretty(xlim, 5), las=1)
axis(2, at=seq(0, 1, by=.25), las=1)

for (rep in seq_len(nrep)) {
    vpctte <- by(output2, factor(output2$id),
        function(x) sample.tte(x$time, x$surv))
    vpctte <- do.call(rbind, vpctte)
    x.vpc <- survfit(Surv(time, status) ~ 1, data=vpctte)
    vpcsurv[rep,] <- summary(x.vpc, times=tgrid, extend=T)$surv
}
# Compute median and 90% PI over simulation repetitions
PI <- apply(vpcsurv, 2, quantile, probs=c(0.05, 0.5, 0.95), na.rm=T)
polygon(c(tgrid, rev(tgrid)), c(PI[1,], rev(PI[3,])), border=NA, col="gray90")
lines(tgrid, PI[2,], col="gray70", lwd=3.5)
lines(x, conf.int=F, lwd=2, col="darkcyan")
mtext("CONSTANT HAZARD MODEL", 3, 1, font=2, cex=1.1)

#+
#' #### R session information

#+ results='markup'
sessionInfo()

