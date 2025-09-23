#' Eccentricity spectrum from an orbital solution
#'
#' Transparent R-implementation of trigonometric expansion of eccentricity.
#' Based on (e, pi) expansion (e.g., `data(La88)`).
#' Ordered by decreasing frequency.
#'
#' @param spectrum_obj `discreteSpectrum` object for trigonometric development
#' @param verbose logical, default FALSE. If TRUE, prints progress messages
#' @return `discreteSpectrum` object for eccentricity development
#' @references Berger A. & M. F. Loutre (1990), Origine des frequences des
#' elements astronomiques intervenant dans le calcul de l'insolation,
#' Bulletin de la Classe des sciences, (1) 45-106.
#' @importFrom gtools combinations
#' @importFrom dplyr bind_rows
#' @importFrom stats aggregate
#' @note Spectrum of m(1+A/2-A^2/8+A/16) with A = e^2/m^2 - 1, m = sqrt(sum(epi$Amp^2))
#' errors are thus of the order of 5/128A^4, which perhaps surprisingly is still of the order
#' of 0.015 / 0.020 in eccentricity.
#' Errors are seem to be always positive (5/128*A^4 may reach 2.0, as A may reach 2.0, but not easy to prove
#' because the next term, in A^5, can be large as well.
#' Using higher order terms helps a bit for capturing the 'low eccentricity' bit, but it actually makes things
#' much worse above 0.04, so that the third order is actually the right deal.
-
#' @export eccentricity_spectrum
eccentricity_spectrum <- function(spectrum_obj, verbose = FALSE) {
  requireNamespace("gtools", quietly = TRUE)

  if (verbose) message("Initializing tables...")

  epi <- as.data.frame(spectrum_obj)
  freq <- epi[["Freq"]]
  amp <- epi[["Amp"]]
  phase <- epi[["Phases"]]

  combine_terms <- function(a, b, omega1, omega2, phi1, phi2, factor = 1) {
    amplitude <- factor * outer(a, b, "*")
    omega <- outer(omega1, omega2, "+")
    phi <- outer(phi1, phi2, "+")
    amplitude <- amplitude[upper.tri(amplitude)]
    omega <- omega[upper.tri(omega)]
    phi <- phi[upper.tri(phi)]
    data.frame(Amp = amplitude, Freq = omega, Phases = phi)
  }

  order_idx <- order(freq, decreasing = TRUE)
  freq <- freq[order_idx]
  amp <- amp[order_idx]
  phase <- phase[order_idx]

  m_norm <- sqrt(sum(amp^2))
  a_norm <- amp / m_norm

  if (verbose) message("Group 1: order 1 terms")

  terms <- combine_terms(a_norm, a_norm, freq, -freq, phase, -phase)
  b_k <- terms$Amp
  g_k <- terms$Freq
  p_k <- terms$Phases

  ek1 <- b_k
  ek2 <- 0.375 * b_k^3
  sum_bk2 <- rep(sum(b_k^2), length(b_k))
  sum_bk2 <- sum_bk2 - b_k^2
  ek3 <- 0.75 * b_k * sum_bk2
  ek <- ek1 + ek2 + ek3

  group <- list()
  group[[1]] <- data.frame(Amp = 1 - 0.25 * sum(b_k^2), Freq = 0, Phases = 0)
  group[[2]] <- data.frame(Amp = ek, Freq = g_k, Phases = p_k)

  order_dec <- order(b_k, decreasing = TRUE)
  b_k <- b_k[order_dec]
  g_k <- g_k[order_dec]
  p_k <- p_k[order_dec]

  group[[3]] <- data.frame(Amp = -0.25 * b_k^2, Freq = 2 * g_k, Phases = 2 * p_k)
  group[[4]] <- combine_terms(b_k, b_k, g_k, g_k, p_k, p_k, -0.5)
  group[[5]] <- combine_terms(b_k, b_k, g_k, -g_k, p_k, -p_k, -0.5)
  group[[6]] <- data.frame(Amp = 0.125 * b_k^3, Freq = 3 * g_k, Phases = 3 * p_k)
  group[[7]] <- combine_terms(b_k^2, b_k, 2 * g_k, g_k, 2 * p_k, p_k, 0.375)
  group[[8]] <- combine_terms(b_k, b_k^2, g_k, 2 * g_k, p_k, 2 * p_k, 0.375)
  group[[9]] <- combine_terms(b_k^2, b_k, 2 * g_k, -g_k, 2 * p_k, -p_k, 0.375)
  group[[10]] <- combine_terms(b_k, b_k^2, -g_k, 2 * g_k, -p_k, 2 * p_k, 0.375)

  max_terms <- 121
  comb_idx <- gtools::combinations(min(max_terms, length(b_k)), 3)
  bklm <- 0.75 * apply(comb_idx, 1, function(v) b_k[v[1]] * b_k[v[2]] * b_k[v[3]])
  gklm1 <- apply(comb_idx, 1, function(v) g_k[v[1]] + g_k[v[2]] + g_k[v[3]])
  gklm2 <- apply(comb_idx, 1, function(v) g_k[v[1]] + g_k[v[2]] - g_k[v[3]])
  gklm3 <- apply(comb_idx, 1, function(v) g_k[v[1]] - g_k[v[2]] + g_k[v[3]])
  gklm4 <- apply(comb_idx, 1, function(v) g_k[v[1]] - g_k[v[2]] - g_k[v[3]])
  pklm1 <- apply(comb_idx, 1, function(v) p_k[v[1]] + p_k[v[2]] + p_k[v[3]])
  pklm2 <- apply(comb_idx, 1, function(v) p_k[v[1]] + p_k[v[2]] - p_k[v[3]])
  pklm3 <- apply(comb_idx, 1, function(v) p_k[v[1]] - p_k[v[2]] + p_k[v[3]])
  pklm4 <- apply(comb_idx, 1, function(v) p_k[v[1]] - p_k[v[2]] - p_k[v[3]])

  group[[11]] <- data.frame(
    Amp = c(bklm, bklm, bklm, bklm),
    Freq = c(gklm1, gklm2, gklm3, gklm4),
    Phases = c(pklm1, pklm2, pklm3, pklm4)
  )

  if (verbose) message("Final development... Filtering terms... Generating partial reconstructions...")

  times <- seq(0, 5000, 2) * 1000
  X <- t(sapply(group, function(dev) {
    sapply(times, function(t)
      sum(m_norm * dev$Amp * cos(dev$Freq * t + dev$Phases)))
  }))

  esinpi <- sapply(times, function(t) sum(amp * sin(freq * t + phase)))
  ecospi <- sapply(times, function(t) sum(amp * cos(freq * t + phase)))
  Y <- sqrt(esinpi^2 + ecospi^2)

  XX <- colSums(X)
  Xts <- ts(XX, deltat = times[2] - times[1])
  Yts <- ts(Y, deltat = times[2] - times[1])

  A <- (Yts / m_norm)^2 - 1
  # the belows lines were used to check accuracy. They are kept for reference.
  # Yts2 <- m_norm * (1 + A / 2 - A^2 / 8 + A^3 / 16)
  # Yts3 <- m_norm * (1 + A / 2 - A^2 / 8 + A^3 / 16 - 5 * A^4 / 128)
  # Yts4 <- m_norm * (1 + A / 2 - A^2 / 8 + A^3 / 16 - 5 * A^4 / 128 + 7 * A^5 / 256)

  if (verbose) message("Removing duplicates...")

  development <- dplyr::bind_rows(group)
  dev1 <- data.frame(Amp = development$Amp, Freq = development$Freq)
  dev2 <- data.frame(Phases = development$Phases, Freq = development$Freq)

  agg1 <- stats::aggregate(. ~ Freq, data = dev1, FUN = sum)
  agg2 <- stats::aggregate(. ~ Freq, data = dev2, FUN = function(x) x[1])
  aff <- cbind(agg1, Phases = agg2$Phases)

  if (verbose) message("Final list...")

  neg_idx <- which(aff$Amp < 0)
  aff$Amp[neg_idx] <- -aff$Amp[neg_idx]
  aff$Freq[neg_idx] <- -aff$Freq[neg_idx]
  aff$Phases[neg_idx] <- pi - aff$Phases[neg_idx]

  order_out <- order(aff$Amp, decreasing = TRUE)
  out <- list(
    Freq = aff$Freq[order_out],
    Amp = aff$Amp[order_out] * m_norm,
    Phases = aff$Phases[order_out]
  )
  class(out) <- "discreteSpectrum"
  return(out)
}
