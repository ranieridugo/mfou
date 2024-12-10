f_ttheta = function(vH, mRho){
  fGij = 
    function(i, j){
      sqrt(beta(vH[i] + 0.5, vH[i] + 0.5) * beta(vH[j] + 0.5, vH[j] + 0.5)) /
        (
          beta(vH[i] + 0.5, vH[j] + 0.5) * sqrt(sin(pi * vH[i]) * sin(pi * vH[j]))
        ) *
        sin(pi * (vH[i] + vH[j])) /
        (cos(pi * vH[i]) + cos(pi * vH[j])) * mRho[i, j]
    }
  
  
  iN = length(vH)
  mTheta = matrix(NA, nrow = iN, ncol = iN)
  for(i in 1 : iN){
    for(j in 1 : iN){ 
      mTheta[i, j] =
        (fGij(i, j) ^ 2 / sqrt(fGij(j, j))) /
        sum(
          sapply(X = 1 : iN,
                 FUN = function(x) (fGij(i = i, j = x) ^ 2 / sqrt(fGij(x, x))))
        )
    }
  }
  return(mTheta)
}
f_spillovers = function(vH, mRho, plot = FALSE){
  mTheta = f_ttheta(vH, mRho)
  dimnames(mTheta) = dimnames(mRho)
  
  dTot = sum(mTheta[upper.tri(mTheta)] + mTheta[lower.tri(mTheta)]) / sum(mTheta)
  iN = sum(mTheta)
  diag(mTheta) = 0
  vReceived = rowSums(mTheta) / iN * 100 # received by others 
  vTransmitted = colSums(mTheta) / iN * 100 # transmitted to others
  vNet = (colSums(mTheta) - rowSums(mTheta)) / iN * 100 # net
  mPair = (t(mTheta) - mTheta) / iN # net pairwise transmitted
  if(plot == TRUE){
    # Plot directional
    par(mfrow = c(3, 1), 
        mar = c(3, 3.5, 2, 0),   # Same margin for all plots, with minimal bottom space
        oma = c(2, 0.5, 0, 0),
        mgp = c(2.5, 1, 0))   # Outer margin for a shared x-axis label
    bpos= barplot(vReceived, 
                  col = "red4", 
                  ylab = "Received by others", 
                  cex.lab = 1.4,  
                  xaxt = "n",                # Hide x-axis for top plots
                  border = NA)
    bpos = barplot(vTransmitted, 
                   col = "blue4", 
                   ylab = "Transmitted to others", 
                   cex.lab = 1.4,  
                   xaxt = "n",                # Hide x-axis for middle plot
                   border = NA)
    net_colors <- ifelse(vNet > 0, "blue4", "red4")
    bpos = barplot(vNet, 
                   col = net_colors, 
                   ylab = "Net transmitted", 
                   cex.lab = 1.4,  
                   xlab = "",                 # No xlab here; we'll use oma for global x-axis
                   names.arg = names(vNet), 
                   las = 2,                   # Rotate x-axis labels for better readability
                   cex.names = 1.4,
                   border = NA)
    par(mfrow = c(1, 1))
    # Plot pairwise
    plot = 
      mPair |> 
      reshape2::melt() |>
      ggplot2::ggplot(ggplot2::aes(Var2, Var1, fill = value * 100)) + # x cols (v2), y rows (v1)
      ggplot2::geom_tile() +
      ggplot2::scale_fill_gradient2(low = "blue4", mid = "white", high = "red4",
                           midpoint = 0, name = "Net") + 
      ggplot2::theme_minimal() +
      ggplot2::labs(x = "To", y = "From") +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(size = 12 * 1.1, angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = ggplot2::element_text(size = 12 * 1.1, vjust = 0.5),
        axis.title.x = ggplot2::element_text(size = 12),
        axis.title.y = ggplot2::element_text(size = 12)
      )
    print(plot)
  }
  
  return(
    list(
      total = dTot,
      received = vReceived,
      transmitted = vTransmitted,
      net = vNet,
      pairwise = mPair
    )
  )
}
