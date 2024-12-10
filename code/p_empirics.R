# import functions and data: depdendencies: dplyr, tidyr, stringr, readr
source("f_estimator_mde.R") # dependencies: pracma, calculus
vSymbols = 
  c("AEX", "AORD", "BFX", "BSESN", "BVSP", "DJI", "FCHI", "FTSE", "GDAXI",
    "HSI", "IBEX", "IXIC", "KS11", "KSE", "MXX", "N225", "NSEI", "RUT", "SPX",
    "SSEC", "SSMI", "S5E")
vCountry = c("EU", "AS", "EU", "AS", "AM", "AM", "EU", "EU", "EU",
            "AS", "EU", "AM", "AS", "AS", "AM", "AS", "AS", "AM",
            "AM", "AS", "EU", "EU")
vOrd = 
  data.frame(symb = vSymbols,
             country = vCountry) |>
  dplyr::arrange(country) |>
  dplyr::pull(symb) |>
  rev() |>
  stringr::str_replace("STOXX50E", "S5E")

tTS =
  read.csv("db_rv_230822.csv") |>
  dplyr::rename(date = X, symbol = Symbol) |>
  dplyr::filter(rv5 != 0) |>
  dplyr::mutate(symbol = stringr::str_remove(symbol, "."),
         lrv5 = log(100 * sqrt(252 * rv5)),
         date = as.Date(date)) |>
  dplyr::select(date, symbol, lrv5) |>
  tidyr::spread(symbol, -date) |>
  dplyr::rename(S5E = STOXX50E) |>
  dplyr::select(vSymbols)

# Perform estimates (or read them)
# vParamsExa = f_mde(mX = tTS, delta = 1/252, type = "exa") # estimate (2 hours required)
# vParamsAsy = f_mde(mX = tTS, delta = 1/252, type = "asy") # estimate (1.5 hours required)
cFile = "mde_exa_l52050_AEX_AORD_BFX_BSESN_BVSP_DJI_FCHI_FTSE_GDAXI_HSI_IBEX_IXIC_KS11_KSE_MXX_N225_NSEI_RUT_SPX_SSEC_SSMI_STOXX50E.csv"
tExa = readr::read_csv(cFile); tAsy = readr::read_csv(stringr::str_replace(cFile, "exa", "asy")) # read
vParamsExa = tExa$val; names(vParamsExa) = gsub("STOXX50E", "S5E", tExa$par) # read
vParamsAsy = tAsy$val; names(vParamsAsy) = gsub("STOXX50E", "S5E", tAsy$par) # read

# Convert to tables
vMu = apply(tTS, 2, mean, na.rm = TRUE)
vAlpha = vParamsExa[paste0("a_", vSymbols)]
vH = vParamsExa[paste0("H_", vSymbols)]
vNu = vParamsExa[paste0("nu_", vSymbols)]
data.frame(
  symbol = vOrd,
  mu = vMu[vOrd],
  alpha = vAlpha[vOrd],
  H = vH[vOrd],
  nu = vNu[vOrd],
  row.names = NULL) |>
  dplyr::mutate(across(where(is.numeric), ~ round(., 3))) 

mRho <- mEta <- mCoherence <- matrix(NA, nrow = length(vSymbols), ncol = length(vSymbols)) 
rownames(mRho) <- colnames(mRho) <- rownames(mEta) <- colnames(mEta) <- vSymbols
mCombn = combn(x = 1 : length(vSymbols), m = 2)
vNamesBiv = apply(X = combn(x = vSymbols, m = 2), MARGIN = 2, FUN = paste0, collapse = "/")
for(i in 1 : ncol(mCombn)){
  iI = mCombn[1, i]
  iJ = mCombn[2, i]
  s1 = unlist(strsplit(vNamesBiv[i], split = "/"))[1]
  s2 = unlist(strsplit(vNamesBiv[i], split = "/"))[2]
  mRho[iI, iJ] = vParamsExa[paste0("rho_", s1, "/", s2)]
  mEta[iI, iJ] = vParamsExa[paste0("eta_", s1, "/", s2)]
  mCoherence[iI, iJ] =  f_coherence(mRho[iI, iJ], mEta[iI, iJ], vH[iI], vH[iJ])
}
mRho[lower.tri(mRho)] = t(mRho)[lower.tri(mRho)]
mEta[lower.tri(mEta)] = t(mEta)[lower.tri(mEta)]
diag(mRho) = 1; diag(mEta) = 0

print(round(mRho[vOrd, vOrd], 3))
print(round(mEta[vOrd, vOrd], 3))

# Goodness of cross-covariance fit (dependency: latex2exp)
iLM = 50
vS1 = c("FCHI", "KS11")
vS2 = c("FTSE", "N225")
vS = c("FCHI", "FTSE", "N225", "KS11")

par(mfrow = c(2, 1), mar = c(2.7, 3, 0.5, 0), mgp = c(1.8, 0.8, 0))
for(i in 1 : 2){
  s1 = vS1[i]; s2 = vS2[i]
    vTheo = f_ccf_theo(vParamsExa[paste0("a_", s1)], vParamsExa[paste0("a_", s2)],
                     vParamsExa[paste0("nu_", s1)], vParamsExa[paste0("nu_", s2)],
                     vParamsExa[paste0("H_", s1)], vParamsExa[paste0("H_", s2)],
                     vParamsExa[paste0("rho_", s1, "/", s2)],  
                     vParamsExa[paste0("eta_", s1, "/", s2)], 
                     iLM)
    mX = cbind(tTS[, s1], tTS[, s2])
    cTitle = paste0(s1, "/", s2)
  vEmp = sapply(X = - iLM : iLM, FUN = ccfl, x = mX[, 1], y = mX[, 2])
  plot(x = - iLM : iLM, y = vEmp, type = "h", 
       ylim = c(0, max(vEmp, vTheo)),
       ylab = latex2exp::TeX("$\\gamma_{i,j}(k)$"),
       xlab = "k")
  lines(x = - iLM : iLM, y = vTheo, col = "red")
  mtext(cTitle, side = 3, line = - 1, adj = 0.01, cex = 0.92, font = 2)
}
par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1), mgp = c(3, 1, 0))

par(mfrow = c(2, 2), mar = c(2.7, 3, 0.5, 0), mgp = c(1.8, 0.8, 0))
vS = c("FCHI", "FTSE", "N225", "KS11")
for(i in 1 : 4){
  s1 = vS[i]
  vTheo = 
    f_acf_theo(vParamsExa[paste0("a_", s1)], vParamsExa[paste0("nu_", s1)],
             vParamsExa[paste0("H_", s1)], lag.max = iLM)
  vEmp = sapply(X = 0 : iLM, FUN = ccfl, x = tTS[, s1], y = tTS[, s1])
  plot(x = 0 : iLM, y = vEmp, type = "h", 
       ylim = c(0, max(vEmp, vTheo)),
       ylab = latex2exp::TeX("$\\gamma_{i,i}(k)$"),
       xlab = "k")
  lines(x = 0 : iLM, y = vTheo, col = "red")
  mtext(s1, side = 3, line = - 1, adj = 0.99, cex = 0.92, font = 2)
}
par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1), mgp = c(3, 1, 0))

par(mfrow = c(2, 1), mar = c(2.7, 3, 0.5, 0), mgp = c(1.8, 0.8, 0))
for(i in 1 : 2){
  s1 = vS1[i]; s2 = vS2[i]
  vX = (0 : iLM) ^ (vParamsExa[paste0("H_", s1)] + vParamsExa[paste0("H_", s2)])
  vEmp1 = sapply(X = 0 : iLM, FUN = ccfl, x = tTS[, s1], y = tTS[, s2])
  vEmp2 = sapply(X = 0 : iLM, FUN = ccfl, x = tTS[, s2], y = tTS[, s1])
  oFitted1 = lm(vEmp1 ~ vX)$fitted
  oFitted2 = lm(vEmp2 ~ vX)$fitted
  plot(x = c(- rev(vX), vX[- 1]), y = c(rev(vEmp2), vEmp1[- 1]),
       type = "h",
       ylab = latex2exp::TeX("$\\gamma_{i,j}(k)$"),
       xlab = latex2exp::TeX("$\\sign{(k)}|k|^{H_i+H_j}$"))
  lines(x = c(- rev(vX), vX[- 1]),
        y = c(rep(NA, iLM), oFitted1), col = "red")
  lines(x = c(- rev(vX), vX[- 1]),
        y = c(rev(oFitted2), rep(NA, iLM)), col = "red")
  mtext(paste0(s1, "-", s2), side = 3, line = - 1, adj = 0.01, cex = 0.92, font = 2)
}
par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1), mgp = c(3, 1, 0))

par(mfrow = c(2, 2), mar = c(2.7, 3, 0.5, 0), mgp = c(1.8, 0.8, 0))
for(i in 1 : 4){
  s1 = vS[i]
  vX = (0 : iLM) ^ (2 * vParamsExa[paste0("H_", s1)])
  vEmp = sapply(X = 0 : iLM, FUN = ccfl, x = tTS[, s1], y = tTS[, s1])
  oFitted = lm(vEmp ~ vX)$fitted
  plot(x = vX, y = vEmp, type = "h", 
       ylim = c(0, max(vEmp, vTheo)),
       ylab = latex2exp::TeX("$\\gamma_{i,i}(k)$"),
       xlab = latex2exp::TeX("$\\sign{(k)}k^{(2H_i)}$"))
  lines(x = vX, y = oFitted, col = "red")
  mtext(s1, side = 3, line = - 1, adj = 0.99, cex = 0.92, font = 2)
}
par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1), mgp = c(3, 1, 0))

par(mfrow = c(2, 1), mar = c(2.7, 3, 0.5, 0), mgp = c(1.8, 0.8, 0))
for(i in 1 : 2){
  s1 = vS1[i]; s2 = vS2[i]
    vTheo = 
      f_ccf_asy_theo(vParamsAsy[paste0("nu_", s1)], vParamsAsy[paste0("nu_", s2)],
                   vParamsAsy[paste0("H_", s1)], vParamsAsy[paste0("H_", s2)],
                   vParamsAsy[paste0("rho_", s1, "/", s2)],  
                   vParamsAsy[paste0("eta_", s1, "/", s2)], 
                   vParamsAsy[paste0("cov_", s1, "/", s2)],
                   iLM)
    mX = cbind(tTS[, s1], tTS[, s2])
    cTitle = paste0(s1, "/", s2)
  vEmp = sapply(X = - iLM : iLM, FUN = ccfl, x = mX[, 1], y = mX[, 2])
  plot(x = - iLM : iLM, y = vEmp, type = "h", 
       ylim = c(0, max(vEmp, vTheo)),
       ylab = latex2exp::TeX("$\\gamma_{i,j}^{a}(k)$"),
       xlab = "k")
  lines(x = - iLM : iLM, y = vTheo, col = "red")
  mtext(cTitle, side = 3, line = - 1, adj = 0.01, cex = 0.92, font = 2)
}
par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1), mgp = c(3, 1, 0))

par(mfrow = c(2, 2), mar = c(2.7, 3, 0.5, 0), mgp = c(1.8, 0.8, 0))
for(i in 1 : 4){
  s1 = vS[i]
  vTheo = 
    f_acf_asy_theo(vParamsAsy[paste0("nu_", s1)], vParamsAsy[paste0("H_", s1)],
                 vParamsAsy[paste0("var_", s1)], lag.max = iLM)
  vEmp = sapply(X = 0 : iLM, FUN = ccfl, x = tTS[, s1], y = tTS[, s1])
  plot(x = 0 : iLM, y = vEmp, type = "h", 
       ylim = c(0, max(vEmp, vTheo)),
       ylab = latex2exp::TeX("$\\gamma^a_{i,i}(k)$"),
       xlab = "k")
  lines(x = 0 : iLM, y = vTheo, col = "red")
  mtext(s1, side = 3, line = - 1, adj = 0.99, cex = 0.92, font = 2)
}
par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1), mgp = c(3, 1, 0))

# rho-based graph (dependency: igraph)
library(igraph)
tRho = 
  tExa |>
  dplyr::filter(stringr::str_detect(par, "rho")) |>
  dplyr::mutate(par_clean = stringr::str_remove(par, "rho_")) |>
  tidyr::separate(par_clean, into = c("V1", "V2"), sep = "/")  |>
  dplyr::select(V1, V2, val) |>
  dplyr::mutate(V1 = ifelse(V1 == "STOXX50E", "S5E", V1)) |>
  dplyr::mutate(V2 = ifelse(V2 == "STOXX50E", "S5E", V2)) |>
  dplyr::rename(rho = val)
data <- graph_from_data_frame(d = tRho, directed = FALSE)
g = simplify(data)
min_rho = min(tRho$rho)
max_rho = max(tRho$rho)
E(g)$weight <- (tRho$rho - min_rho + 0.001)/(max_rho - min_rho + 0.001)
vColor <- dplyr::case_when(
  vCountry == "EU" ~ "blue4",
  vCountry == "AM" ~ "red4",
  vCountry == "AS" ~ "grey"
)
V(g)$color = vColor
V(g)$label.cex = 0.65
set.seed(105) 
plot(g, 
     layout = layout_with_fr(g), 
     edge.length = E(g)$weight,
     edge.width = E(g)$weight,
     vertex.label.color = "white",
     vertex.size = 22,                  # Increases the size of the vertices
     vertex.frame.width = 0,
     vertex.color = V(g)$color)
legend("topleft", 
       legend = c("EU", "AM", "AS"), 
       col = c("blue4", "red4", "grey"),
       pch = rep(15, 3))

# spillovers (dependency: ggplot2)
source("f_spillovers.R")
# vParamsCau = f_mde(tTS, type = "cau") # estimate (2 hours requires)
tCau = readr::read_csv(stringr::str_replace(cFile, "exa", "cau")) # read
# vParamsCau = f_mde(tTS, type = "cau") # estimate (2 hours requires)
vParamsCau = tCau$val; names(vParamsCau) = gsub("STOXX50E", "S5E", tCau$par)
vH = vParamsCau[paste0("H_", vSymbols)]
for(i in 1 : ncol(mCombn)){
  iI = mCombn[1, i]
  iJ = mCombn[2, i]
  s1 = unlist(strsplit(vNamesBiv[i], split = "/"))[1]
  s2 = unlist(strsplit(vNamesBiv[i], split = "/"))[2]
  mRho[iI, iJ] = vParamsCau[paste0("rho_", s1, "/", s2)]
}
mRho[lower.tri(mRho)] = t(mRho)[lower.tri(mRho)]
f_spillovers(vH[paste0("H_", vOrd)], mRho[vOrd, vOrd], plot = TRUE)