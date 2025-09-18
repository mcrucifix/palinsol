# require(devtools)
# load_all()

#' Caloric insolation
#'
#' Computes caloric summer insolation for a given astronomical configuration
#' and latitude.
#'
#' The caloric summer is a notion introduced by M. Milankovitch. It is defined
#' as the halve of the tropical year during for which daily mean insolation are
#' greater than all days of the other halves. The algorithm is an original
#' algorithm by M. Crucifix, but consistent with earlier definitions and
#' algorithms by A. Berger (see examples). Do not confuse this Berger (1978)
#' reference with the Berger (1978), J. Atm. Sci. of the astronomical solution.
#'
#' @param orbit Output from a solution, such as \code{ber78}, \code{ber90} or
#' \code{la04}
#' @param lat latitude
#' @param ... Other arguments passed to Insol
#' @return Time-integrated insolation in kJ/m2 during the caloric summer.
#' @author Michel Crucifix, U. catholique de Louvain, Belgium.
#' @references Berger (1978) Long-term variations of caloric insolation
#' resulting from the earth's orbital elements, Quaternary Research, 9, 139 -
#' 167.
#' @examples
#'
#' ## reproduces Table 2 of Berger 1978
#' lat <- seq(90, 0, -10) * pi/180. ## angles in radians.
#' orbit_1 = ber78(0)
#' orbit_2 = orbit_1
#' orbit_2 ['eps'] = orbit_2['eps'] + 1*pi/180.
#
#' Table2 <-  sapply(lat, function(x) c(lat = x * 180/pi,
#'                           calins(orbit_2, lat=x, S0=1365) / (4.18 * 1e1)
#'                         - calins(orbit_1, lat=x, S0=1365) / (4.18 * 1e1) ) )
#' data.frame(t(Table2))
#' # there are still some differences, of the order of 0.3 %, that are probably related to
#' # the slightly different methods.
#' # 41.8 is the factor from cal/cm2 to  kJ/m2
#'
#' @export calins
calins <- function(orbit, lat, degrees = FALSE, ...) {
  if (degrees == TRUE) {
    lat <- lat * pi / 180
  }

  # supply difference of insolation between a given time (mean longitude + 80.5) and
  # 6 months later

  diff_insol <- function(orbit, day, lat, ...) {
    l1 <- day2l(orbit, day)
    l2 <- day2l(orbit, day + 180.)
    Insol(orbit = orbit, long = l2, lat = lat, ...) -
        Insol(orbit = orbit, long = l1, lat = lat, ...)
  }

  # this first guess will find the location of the maximu insolation
  # this depends on the longitude of the perihelion,
  first_guess <- function(orbit, lat) {
    modvarpi <- (orbit["varpi"] + pi / 2) %% (2 * pi)
    ifelse(modvarpi > pi,
      (170 - 90 * tanh(1.7 * lat) + 360) %% 360,
      (-10 + 90 * tanh(1.7 * lat) + 360) %% 360
    )
  }

  fg <- first_guess(lat)
  # we use the uniroot function to find the root
  root_result <- stats::uniroot(
    function(day) {
      diff_insol(orbit, day, lat)
    },
    c(fg - 60, fg + 60),
    extendInt = "yes"
  )

  if (root_result$root) {
    d1 <- root_result$root
    d2 <- root_result$root + 180.
    sol <- Insol_d1d2(orbit, d1 = d1, d2 = d2, lat = lat, avg = FALSE)
  } else {
    sol <- NULL
  }

  return(sol)
}

# TO DO : check the notion of caloric equator
