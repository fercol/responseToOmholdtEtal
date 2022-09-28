# ============================== CODE METADATA =================================
# FILE NAME: ResponseToOmholdtEtal.R
# AUTHORS: Fernando Colchero
# DATE CREATED: 2022-09-28
# DATE MODIFIED: 
# DESCRIPTION: Reproducible code for Response to comment on “Slow and 
#              negligible senescence among testudines challenges evolutionary 
#              theories of senescence”. 
# DETAILS: This code reproduces the least squares fit of the model propose in
#          (2) and (3), to survival curves estimated in (1) and compares the 
#          resulting hazard rates to those estimated in (1).
# REFERENCES: 
#         1. R. da Silva, D.A. Conde, A. Baudisch, F. Colchero. Slow and
#            negligible senescence among testudines challenges evolutionary
#            theories of senescence. Science 376, 1466–1470 (2022).
#         2. S.W. Omholt, O. Omholt, T.B.L. Kirkwood. Comment on “Slow and
#            negligible senescence among testudines challenges evolutionary
#            theories of senescence” Science (under review)
#         3. S.W. Omholt, T.B.L. Kirkwood. Aging as a consequence of
#            selection to reduce the environmental risk of dying. PNAS 
#            118, e2102088118 (2021).
# ================================ CODE START ==================================
# ========================= #
# ==== GENERAL SETUP: =====
# ========================= #
# set working directory (change accordingly):
setwd("~/FERNANDO/PROJECTS/1.ACTIVE/TestudinesSenescence/analysisTestudinesSenes/code/responseToOmholdtEtal/")

# Load survival and mortality (i.e. hazard rate) curves in (1) 
# for T. scripta and T. graeca:
load("demoList.RData")

# ==================== #
# ==== FUNCTIONS: ====
# ==================== #
# Hazard based on (2) and (3):
hxOK <- function(x, mu0, alpha, kappa) {
  c <- mu0 - alpha
  a <- alpha * (1 + kappa)
  b <- log(1 + kappa)
  return(c + a * exp(b * x))
}

# Resulting survival function from (2) and (3):
SxOK <- function(x, mu0, alpha, kappa) {
  c <- mu0 - alpha
  a <- alpha * (1 + kappa)
  b <- log(1 + kappa)
  return(exp(-c * x + a/b * (1 - exp(b * x))))
}

# Least squares function to fit model in (2) and (3) to survival curves:
lsqOKpars <- function(theta) {
  mu0 <- abs(theta[1])
  alpha <- abs(theta[2])
  kappa <- theta[3]
  Sxfit <- SxOK(xref, mu0, alpha, kappa)
  return(sum((Sxfit - Sxref)^2))
}

# =================== #
# ==== ANALYSES: ====
# =================== #
# Number of species:
nsp <- nrow(modelTab)

# Output of fitted values:
fitList <- list()
paramList <- list()
for (isp in 1:nsp) {
  # species:
  species <- colnames(demoList$surv)[isp + 1]
  
  # Extract survival and mortality curves from results in (1):
  idref <- which(!is.na(demoList$surv[[species]]))
  Sxref <- demoList$surv[[species]][idref]
  muxref <- demoList$mort[[species]][idref]
  xref <- demoList$surv$Age[idref]

  # initial parameters for :
  kappa <- 0.1
  mu0 <- min(muxref)
  alpha <- 0.01
  inThe <- c(mu0, alpha, kappa)
  
  # Recursive LSQ fitting (to ensure minimum LSQ):
  for (ii in 1:20) {
    outlsq <- optim(inThe, lsqOKpars)
    inThe <- outlsq$par
    inThe[1:2] <- abs(inThe[1:2])
  }
  
  # Extract hazard rate parameters:
  theFit <- outlsq$par
  theFit[1:2] <- abs(theFit[1:2])
  names(theFit) <- c("mu0", "alpha", "kappa")
  
  # Caluculate resulting survival and hazard rate:
  Sxfit <- SxOK(xref, theFit[1], theFit[2], theFit[3])
  muxfit <- hxOK(xref, theFit[1], theFit[2], theFit[3])
  
  # Fill up output list:
  fitList[[species]] <- data.frame(x = xref, Sxref = Sxref, muxref = muxref,
                                  Sxfit = Sxfit, muxfit = muxfit)
  paramList[[species]] <- theFit
}

# Print parameter output:
print(paramList)

# ======================= #
# ==== PLOT RESULTS: ====
# ======================= #
# Species list (without underscore)
spList <- gsub("_", " ", colnames(demoList$surv)[-1])

# X and Y axes limits:
ylimS <- c(0, 1)
ylimM <- c(0, 0.2)

# Line specifications:
lwdref <- 4
lwdfit <- 2
ltyref <- 2

# Colors:
colref <- 'dark green'
colfit <- 'orange'

# Matrix and dimensions for layout matrix:
laymat <- rbind(c(2, 4, 3, 5), c(2, 6, 3, 7), 1)
widths <- c(0.1, 1, 0.1, 1)
heights <- c(0.7, 0.7, 0.2)

# Produce plot:
layout(mat = laymat, widths = widths, heights = heights)

# Xlabel:
par(mar = c(0, 1, 0, 1))
plot(c(0, 1), c(0, 1), col = NA, axes = FALSE, xlab = "", ylab = "")
text(0.5, 0.85, "Age", cex = 2, xpd = NA)
legend('bottomright', legend = c("da Silva et al. (2)", 
                                 "Fitted based on Omholt et al. (1)"), 
       col = c(colref, colfit), pch = NA, lwd = c(lwdref, lwdfit), 
       lty = c(ltyref, 1), bty = 'n', seg.len = 8)

# Y axis labels:
for (ii in 1:2) {
  par(mar = c(1, 0, 1, 0))
  plot(c(0, 1), c(0, 1), col = NA, axes = FALSE, xlab = "", ylab = "")
  text(0.5, 0.5, c("Survival", "Hazard rate")[ii], cex = 2, xpd = NA, srt = 90)
  
}

par(mar = c(1, 1, 1, 1))
# Survival:
ilet <- 0
for (isp in 1:nsp) {
  xlim <- c(0, max(fitList[[isp]]$x))
  
  ilet <- ilet + 1
  plot(xlim, ylimS + diff(ylimS) * c(0, 0.1), col = NA, xlab = "", ylab = "", 
       axes = FALSE)
  lines(fitList[[isp]]$x + 0.5, fitList[[isp]]$Sxref, col = colref, 
        lwd = lwdref, lty = ltyref)
  lines(fitList[[isp]]$x + 0.5, fitList[[isp]]$Sxfit, col = colfit, 
        lwd = lwdfit)
  lines(xlim + diff(xlim) * c(-0.05, 0.05), rep(ylimS[1], 2), lwd = 2)
  lines(rep(xlim[1], 2), ylimS + diff(ylimS) * c(-0.05, 0.05), lwd = 2)
  
  text(xlim[1] + diff(xlim) * 0.2, ylimS[2] + diff(ylimS) * 0.05, 
       labels = LETTERS[ilet], font = 2, cex = 2)
  text(xlim[2] - diff(xlim) * 0.05, ylimS[2] + diff(ylimS) * 0.05, 
       labels = spList[isp], font = 3, cex = 1.5, adj = 1)
  
  # Mortality:
  ilet <- ilet + 1
  plot(xlim, ylimM + diff(ylimM) * c(0, 0.1), col = NA, xlab = "", ylab = "", 
       axes = FALSE)
  lines(fitList[[isp]]$x + 0.5, fitList[[isp]]$muxref, col = colref, 
        lwd = lwdref, lty = ltyref)
  lines(fitList[[isp]]$x + 0.5, fitList[[isp]]$muxfit, col = colfit, 
        lwd = lwdfit)
  lines(xlim + diff(xlim) * c(-0.05, 0.05), rep(ylimM[1], 2), lwd = 2)
  lines(rep(xlim[1], 2), ylimM + diff(ylimM) * c(-0.05, 0.05), lwd = 2)
  
  text(xlim[1] + diff(xlim) * 0.2, ylimM[2] + diff(ylimM) * 0.05, 
       labels = LETTERS[ilet], font = 2, cex = 2)
  
}
