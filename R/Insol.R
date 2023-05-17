# Copyright (c) 2012 Michel Crucifix <michel.crucifix@uclouvain.be>

# Permission is hereby granted, free of charge, to any person obtaining
# a copy of this software and associated documentation files (the
# "Software"), to deal in the Software without restriction, including
# without limitation the rights to use, copy, modify, merge, publish,
# distribute, sublicense, and/or sell copies of the Software, and to
# permit persons to whom the Software is furnished to do so, subject

# the following conditions:

# The above copyright notice and this permission notice shall be
# incluudedin all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND INFRINGEMENT
# IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR

# CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
# TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
# SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

# When using this package for actual applications, always
# cite the authors of the original insolation solutions 
# Berger, Loutre and/or Laskar, see details in man pages


# ------------------------------------------------------------------
# R Code developed for R version 2.15.2 (2012-10-26) -- "Trick or Treat"
# ------------------------------------------------------------------ 

 SIDERAL_YEAR = 365.25636 
 TROPIC_YEAR  = 365.24219876
 YEAR = TROPIC_YEAR 
 ## delibarate difference with 
 ## Berger et al.(2010), who consider the SIDERAL_YEAR

 
 # Returns daily mean incoming solar insolation after Berger (1978) 

 # eps  : obliquity (radians)
 # varpi: longitude of perihelion (radians)
 # e : eccentricity 
 # long : true solar longitude 
 #        (radians; = pi/2 for summer solstice )
 #                     pi   for autumn equinox  )
 #                     3* pi/2 for winter solstice )
 #                     0 for spring equinox )
 # lat : latitude
 # orbit is a list containing eps, varpi and e (see Ber78)
 # S0 : Total solar irradiance (default : 1365 W/m^2)
 # returns : daily mean incoming solar radiation at the top of the atmosphere
 #           in the same unit as S0. 
 # H : if NULL (default): compute daily mean insolation
 # else: hour angle (in radians) at which insolation is being computed


#' Computes incoming solar radiation (insolation)
#' 
#' Computes incoming solar radiation (insolation) for a given astronomical
#' configuration, true solar longitude and latitude
#' 
#' True solar longitude is measured in radians: \tabular{ll}{ pi/2 \tab for
#' June solstice\cr pi \tab for September equinox\cr 3 * pi/2 \tab for December
#' solstice\cr 0 \tab for Spring equinox\cr } It may be obtained for a given
#' day in the year using the function \code{day2l}.
#' 
#' @param orbit Output from a solution, such as \code{ber78}, \code{ber90} or
#' \code{la04}
#' @param long true solar longitude
#' @param lat latitude
#' @param S0 Total solar irradiance
#' @param H Sun hour angle, in radians
#' @return Daily-mean insolation (assuming fixed astronomical parameters during
#' a true solar day) if `H` is null.  Otherwise, insolation at specified hour
#' angle (H = 0 at noon, H = pi at midnight).
#' @author Michel Crucifix, U. catholique de Louvain, Belgium.
#' @references Berger, A. L. (1978).  Long-term variations of daily insolation
#' and Quaternary climatic changes, J. Atmos. Sci., 35, 2362-2367.
#' @keywords misc
#' @examples
#' 
#' 
#' ## make a little wrapper, with all default values
#' 
#' insolation <- function(times, astrosol=ber78,...)
#'   sapply(times, function(tt) Insol(orbit=astrosol(tt)))
#' 
#' tts <- seq(from = -400e3, to = 0, by = 1e3)
#' isl <- insolation(tts, ber78)
#' plot(tts, isl, typ='l')
#' 
#' 
#' @export Insol
 Insol <- function (orbit,long=pi/2, lat=65*pi/180,S0=1365, H=NULL)
 {
  # orbit is a list that contains the following variables: varpi, eps, and ecc
  varpi <- NULL
  eps <- NULL
  ecc <- NULL
  for (i in names(orbit)) assign(i,orbit[[i]])
  nu <- long - varpi
  rho <- (1-ecc^2)/(1+ecc*cos(nu))
  sindelta <- sin(eps)*sin(long)
  cosdelta <- sqrt(1-sindelta^2)
  sinlatsindelta <- sin(lat)*sindelta
  coslatcosdelta <- cos(lat)*cosdelta
  if (is.null(H))
  {
    cosH0 <- pmin(pmax(-1,-sinlatsindelta/coslatcosdelta),1)
    sinH0 <- sqrt(1-cosH0^2)
    H0 <- acos(cosH0)
    insol <- S0/(pi*rho^2)*(H0*sinlatsindelta+coslatcosdelta*sinH0)
  }
  else 
  {
    insol <- max(0, S0/(rho^2)*(sinlatsindelta+coslatcosdelta*cos(H)))
  }
  return(as.numeric(insol))
  }
 
 ## time increment corresponding a tsl increment
 .dtdnu <- function (orbit,long=pi/2) {
  nu <- long - orbit['varpi']
  ecc <- orbit['ecc']
  xec <- ecc*ecc
  rho <- (1-xec)/(1+ecc*cos(nu))
  .dtdnu <- rho^2/sqrt(1.-xec)
 }
 # Provides an insolation times series
 # astrosol = astronomical solution (defaults to Ber78)
 # times = times in yr epoch 1950.0
 # ...   = any argument passed to Insol
 #InsolWrapper <- function(times=times,astrosol=ber78,...)
 #  sapply (times,function(tt)  Insol(astrosol(tt),...) )

 ## caloric_insolation
 ## integrated insolation over the 180 days receiving above median insolation


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
#' 
#' T <-  sapply(lat, function(x) c(lat = x * 180/pi, 
#'                           calins(orbit_2, lat=x, S0=1365) / (4.18 * 1e1)
#'                         - calins(orbit_1, lat=x, S0=1365) / (4.18 * 1e1) ) )
#' data.frame(t(T))
#' # there are still some differences, of the order of 0.3 %, that are probably related to
#' # the slightly different methods. 
#' # 41.8 is the factor from cal/cm2 to  kJ/m2
#' 
#' @export calins
 calins <- function (orbit,lat=65*pi/180,...) {
    ins   <- sapply(seq(1:360)*pi/180, function(x) Insol(orbit,long=x, lat=lat,...))
    dt    <- sapply(seq(1:360)*pi/180, function(x) .dtdnu (orbit,long=x))
   
    is    <- sort(ins,decreasing=TRUE,index.return=TRUE)$ix
    cs    <- cumsum(dt[is])
    ### BUG FIX ON 20.05.2021
    ## is    <- which(cs <= 180)           ## the 180 days whose cumulative length is half total
    i_halfyear   <- is[which(cs <= 180)]           ## the 180 days whose cumulative length is half total
                                        ## year length, picking days by decreasing order of 
                                        ## insolation. 
    XCORR = 86.4 *  YEAR / 360
    sum(ins[i_halfyear]*dt[i_halfyear])* XCORR    ## result in kJ
   }


 ## integrated insolation over the 360 days receiving insolation above a threshold


#' Integrated insolation for all days exceeding a threshold
#' 
#' Integrated insolation over the part during which daily-mean insolation
#' exceeds a threshold, expressed in W/m2
#' 
#' Algorithm is by M. Crucifix, but the idea of thresholded insolation is due
#' to Huybers and Tziperman (2008), reference below.
#' 
#' @param lat latitude
#' @param orbit Output from a solution, such as \code{ber78}, \code{ber90} or
#' \code{la04}
#' @param threshold threshold insolation ,in W/m2
#' @param ... other arguments to be passed to Insol
#' @return Time-integrated insolation in kJ/m2 . The quantity is calculated by
#' brute-force integration with a 1-degree time-step in true solar longitude
#' and this can be quite slow if long series are to be calculated.
#' @author Michel Crucifix, U. catholique de Louvain, Belgium.
#' @references P. Huybers and E. Tziperman (2008), Integrated summer insolation
#' forcing and 40,000-year glacial cycles: The perspective from an
#' ice-sheet/energy-balance model, Paleoceanography, 23.
#' @export thrins
 thrins <- function (lat=65*pi/180,orbit,threshold=400,...)
   {
    ins   <- sapply(seq(1:360)*pi/180, function(x) Insol(orbit,long=x, lat=lat,...))
    dt    <- sapply(seq(1:360)*pi/180, function(x) .dtdnu (orbit,long=x, ...))
    
    is    <- which(ins >= threshold)
    XCORR = 86.4 *  YEAR / 360
    sum(ins[is]*dt[is])* XCORR   ## result in kJ
   }



 ## time-integrated between two true solar longitude bounds


#' Time-integrated insolation
#' 
#' Computes time-integrated incoming solar radiation (Insol) either between
#' given true solar longitudes (\code{Insol_l1l2}) or days of year
#' (\code{Insol_d1d2}) for a given orbit and latitude
#' 
#' All angles input measured in radians.
#' 
#' Note that in contrast to Berger (2010) we consider the tropic year as the
#' reference, rather than the sideral year, which partly explains some of the
#' small differences with the original publication
#' 
#' @aliases Insol_l1l2 
#' @param orbit Output from a solution, such as \code{ber78}, \code{ber90} or
#' \code{la04}
#' @param lat latitude
#' @param l1 lower true solar longitude bound of the time-integral
#' @param l2 upper true solar longitude bound of the time-integral
#' @param d1 lower calendar day (360-day year) of the time-integral
#' @param d2 upper calendar day (360-day year) of the time-integral
#' @param avg performs a time-average.
#' @param ell uses elliptic integrals for the calculation (much faster)
#' @param ... other arguments to be passed to \code{Insol}
#' @return Time-integrated insolation in kJ/m2 if \code{avg=TRUE}, else
#' time-average in W/m2
#' @author Michel Crucifix, U. catholique de Louvain, Belgium.
#' @references Berger, A., Loutre, M.F. and Yin Q. (2010), Total irradiation
#' during any time interval of the year using elliptic integrals, Quaternary
#' Science Reviews, 29, 1968 - 1982, doi:10.1016/j.quascirev.2010.05.007
#' @examples
#' 
#' ## reproduces Table 1a of Berger et al. 2010:
#' lat <- seq(85, -85, -10) * pi/180. ## angles in radians. 
#' orbit=c(eps=  23.446 * pi/180., ecc= 0.016724, varpi= (102.04 - 180)* pi/180. )
#' T <-  sapply(lat, function(x) c(lat = x * 180/pi, 
#'         m1 =  Insol_l1l2(orbit, 0, 70 * pi/180, lat=x, ell= TRUE, S0=1368) / 1e3,
#'         m2 =  Insol_l1l2(orbit, 0, 70 * pi/180, lat=x, ell=FALSE, S0=1368) / 1e3) ) 
#' data.frame(t(T))
#'  ## reproduces Table 1b of Berger et al. 2010:
#' lat <- c(85, 55, 0, -55, -85) * pi/180. ## angles in radians. 
#' T <-  sapply(lat, function(x) c(lat = x * 180/pi, 
#'          m1 =  Insol_l1l2(orbit, 30 * pi/180. , 75 * pi/180, 
#'                lat=x, ell= TRUE, S0=1368) / 1e3,
#'          m2 =  Insol_l1l2(orbit, 30 * pi/180. , 75 * pi/180, 
#'                lat=x, ell=FALSE, S0=1368) / 1e3) ) 
#'  ## reproduces Table 2a of Berger et al. 2010:
#' lat <- seq(85, -85, -10) * pi/180. ## angles in radians. 
#' 
#' ## 21 march in a 360-d year. By definition : day 80 = 21 march at 12u
#' d1 = 79.5 
#' d2 = 79.5 + (10 + 30 + 30 ) * 360/365.2425 ## 30th May in a 360-d year
#' 
#' T <-  sapply(lat, function(x) c(lat = x * 180/pi, 
#'         m1 =  Insol_d1d2(orbit, d1,d2, lat=x, ell= TRUE, S0=1368) / 1e3,
#'         m2 =  Insol_d1d2(orbit, d1,d2, lat=x, ell= FALSE, S0=1368) / 1e3))
#'                           
#' ## I did not quite get the same results as on the table 
#' ## on this one; probably a matter of calendar
#' ## note : the authors in fact used S0=1368 (pers. comm.) 
#' ## 1366 in the paper is a misprint
#' 
#' data.frame(t(T))
#' 
#' ## annual mean insolation at 65N North, as a function of the longitude of the perihelion
#' ## (expected to be invariant)
#'
#' varpis <- seq(0,360,10)*pi/180.  
#' sapply(varpis, function(varpi)
#'    {   orbit=c(eps=  23.446 * pi/180., ecc= 0.016724, varpi= varpi )
#'        amean <- Insol_l1l2 (orbit, lat=65*pi/180., avg=TRUE)
#'        return(amean)
#'    })

#' @export Insol_l1l2

Insol_l1l2 <- function (orbit,l1=0,l2=2*pi,lat=65*pi/180,avg=FALSE,ell=TRUE,...)
   {
    # parameters: orbit : supplied by orbit calculator; e.g. : ber78 or ber90
    # l1 and l2 : longitudes bonds in radians. Defaults to annual average.
    # discretize longitude intreval in N intervals
    # avg : supplies an average insolation
    # ell : defaults to TRUE, use elliptic integrals for calculation (much faster)
    #       rather than trapeze rule integral.  Currently incompatible
    #       with avg=TRUE (this can be fixed later)

     if (ell &&  requireNamespace("gsl", quietly = TRUE))
      ## use elliptic integrals if required and available
      {
          if (l1 < l2) 
          { 
          INT = W(lat, orbit['eps'], orbit['ecc'],l2,...) - 
                W(lat, orbit['eps'], orbit['ecc'], l1,...)} 
          else
          {
           INT = W(lat, orbit['eps'], orbit['ecc'], 2*pi,...) -
                 W(lat, orbit['eps'], orbit['ecc'], l1,...) +
                 W(lat, orbit['eps'], orbit['ecc'], l2,...)
           }
           if (avg) { 
                DT <- l2day(orbit,l2) - l2day(orbit,l1)
                if (DT <= 0) DT = DT+360.
                XCORR = 86.4 *  YEAR / 360
                INT = INT / (DT*XCORR) ## result in W/m2
                }
      }
      else
      ## integration using trapeze rule
      {
        Dl = ((l2-l1) %% (2*pi))
        if (Dl == 0) Dl=2*pi
        N =  1*ceiling(Dl* 180/pi)
        dl = Dl/N
        L  = l1+(0:N)*dl
        ins   <- sapply(L, function(x) Insol(orbit=orbit,long=x, lat=lat,...))
        dt    <- sapply(L, function(x) .dtdnu (orbit=orbit,long=x)) * 180./pi
   

        is    <- ins*dt
        XCORR = 86.4 *  YEAR / 360
        INT = (sum(is[2:N]) +  0.5 *is[1] +  0.5 *is[N+1]) * dl * XCORR  ## result in kJ
        
        if (avg) {
         DT =  (sum(dt[2:N]) + 0.5 * dt[1] + 0.5 * dt[N+1]) * dl * XCORR
         INT = INT / DT ## result in W/m2
         }
      }
      as.numeric(INT)
    }



#' Converts calendar day into true solar longitude and vice-versa
#' 
#' Converts calendar day into true solar longitude for a given astronomical
#' configuration and vice-versa
#' 
#' The 360-d calendar is a conventional calendar, for which day 80 is the day
#' of NH spring equinoxe. The tropic year, which in reality is 365.24219876 *
#' 86400 seconds was the practical reference to define the Gregorian Calendar
#' since this is the time needed to go through all the seasons. More discussion
#' of calendars and conversions in Berger et al. (2010) appendix D.
#' 
#' The \code{day2l} and \code{l2day} is based on algoritms given in Berger
#' (1978), but which can be traced back to expansions of the mean and true
#' anomaly by Brouwer and Clemente (1961), pp. 65 and 77 (see code for further
#' details).
#' 
#' @aliases day2l 
#' @param orbit Output from a solution, such as \code{ber78}, \code{ber90} or
#' \code{la04}
#' @param l true solar longitude, in radians
#' @param day calendar day, in a 360-d year
#' @return day of year (360-d cal.) or true solar longitude (in radians).
#' @author Michel Crucifix, U. catholique de Louvain, Belgium.
#' @references Brouwer D. and G. M. Clemence, (1961), Methods of celestial
#' mechanics, Academic Press, New York.
#' 
#' Berger, (1978) Long-term variations of daily insolation and Quaternary
#' climatic changes, J. Atmos. Sci., 35, 2362-2367 1978,
#' doi:10.1175/1520-0469(1978)035<2362:LTVODI>2.0.CO;2
#' 
#' Berger, A. Loutre, M.F. and Yin Q. (2010), Total irradiation during any time
#' interval of the year using elliptic integrals, Quaternary Science Reviews,
#' 29, 1968 - 1982, doi:10.1016/j.quascirev.2010.05.007
#' @examples
#' 
#' ## date of perihelion throughout today
#' orbit=c(eps=0.409214, ecc=0.01672393, varpi=4.92251)
#' date_of_perihelion(orbit)
#' ## date of winter solstice)
#' l2day(orbit, 270*pi/180.)
#' 
#' @export day2l
day2l  <- function (orbit,day)
   { 
    ## converts day to longitude.
    ## source : Berger 78, from Brower and Clemence.
    ## day using a 360-d calendar
    ## here :  day of spring equinox is in fact conventionally 80
    
    ecc = orbit['ecc']
    varpi= orbit['varpi']
    # definitions for speed-up
    xee= ecc*ecc
    xec = xee * ecc
    xse= sqrt(1.-xee)
    # true anomaly of vernal equinox point 
    xlp= - varpi 
    # mean anomaly  of the vernal equinox point
    lm =  xlp - 2.*((ecc/2 + xec/8)*(1+xse)*sin(xlp) - 
          xee/4*(0.5+xse)*sin(2*xlp) +
          xec/8*(1/3+xse)*sin(3*xlp) )

    # mean anomaly  of the point considered
    M = (day- 80) * pi/ 180. + lm 

    # TRUE anomaly of the point considered
    V = M + (2*ecc-xec/4.)*sin(M) +
          5/4*xee*sin(2*M) +
          13/12*xec*sin(3*M)

    # TRUE longitude of the point considered

    L = ( V + varpi ) %% (2*pi)
    L
   }

#' @describeIn day2l Converts true solar longitude into calendar day
#' @export
l2day <- function (orbit,l)
    ## source :  Brouwer and Clemence
   { 
    ecc = orbit['ecc']
    varpi= orbit['varpi']
    # definitions for speed-up
    xee= ecc*ecc
    xec = xee * ecc
    xse= sqrt(1.-xee)
    # true anomaly
    xlp= - varpi 
    # mean anomaly  of the vernal equinox point
    lm =  xlp - 2.*((ecc/2 + xec/8)*(1+xse)*sin(xlp) - 
          xee/4*(0.5+xse)*sin(2*xlp) +
          xec/8*(1/3+xse)*sin(3*xlp) )

    # true anomaly of the point considered
    V = (l + xlp) %% (2*pi)

    # mean anomaly  of the point considered
    M =  V - 2.*((ecc/2 + xec/8)*(1+xse)*sin(V) - 
          xee/4*(0.5+xse)*sin(2*V) +
          xec/8*(1/3+xse)*sin(3*V) )

    # anomaly in deg. elapsed between vernal equinox point and point

    DAY = ( 80 + (M - lm)  * 360.0/(2*pi) ) %% (360.0)
    DAY

   }


#' @describeIn day2l Returns date of perihelion
#' @export 
date_of_perihelion <- function(orbit)
 {
    ecc = orbit['ecc']
    varpi= orbit['varpi']
    # definitions for speed-up
    xee= ecc*ecc
    xec = xee * ecc
    xse= sqrt(1.-xee)
    # true anomaly
    xlp= - varpi 
    # mean anomaly  of the vernal equinox point
    lm =  xlp - 2.*((ecc/2 + xec/8)*(1+xse)*sin(xlp) - 
          xee/4*(0.5+xse)*sin(2*xlp) +
          xec/8*(1/3+xse)*sin(3*xlp) )

    # mean anomaly  of the point considered
    M =  0.

    # anomaly in deg. elapsed between vernal equinox point and perihelion passage

    DAY = ( 80 + (M - lm)  * 360.0/(2*pi) ) %% (360.0)
    names(DAY) <- 'day'
    DAY

   }


#' @describeIn Insol_l1l2 Mean and integrated insolation over an interval bounded by calendar days
#' @aliases Insol_d1d2
#' @export Insol_d1d2

Insol_d1d2 <- function (orbit,d1,d2,lat=65*pi/180,avg=FALSE,...)
## as insol but given days rather than longitudes

   {
      l1 = day2l(orbit,d1)
      l2 = day2l(orbit,d2)
      Insol_l1l2(orbit,lat=lat,l1,l2,avg=avg,...)
    }  


# if (isTrue(getOption('debug')) && interactive())
# { t <- seq(-1e6,0,by=1e3)
#   F <- InsolWrapper(t)



#' Milankovitch graph for a given astronomical configuration
#' 
#' Computes the distrubition in latitude and longitude of incoming solar
#' radiation, known as a Milankovitch graph, with possibility of plotting with
#' a dedicated plot function
#' 
#' 
#' @param orbit Output from a solution, such as \code{ber78}, \code{ber90} or
#' \code{la04}
#' @param S0 Total solar irradiance
#' @param lat latitudes, passed as an array
#' @param long true solar longitudes, passed as an array
#' @param deg If true : the axes of the Milankovitch object are expressed in
#' degrees.  Inputs are always in radians
#' @return A object of Milankovitch class, which may be plotted using the
#' regular plot function
#' @note The polar night option may not be bullet-proof for exotic obliquities
#' @author Michel Crucifix, U. catholique de Louvain, Belgium.
#' @references Berger, A. L. (1978).  Long-term variations of daily insolation
#' and Quaternary climatic changes, J. Atmos. Sci., 35, 2362-2367.
#' @keywords misc
#' @examples
#' 
#' orbit <- c(eps=0.409214, ecc=0.01672393, varpi=4.92251)
#' M <- Milankovitch(orbit)
#' plot(M, plot=contour)
#' plot(M, plot=contour, month=FALSE)
#' 
#' @export Milankovitch
Milankovitch <- function(orbit, S0=1365, lat=seq(-pi/2, pi/2, l=73), long=seq(0, 2*pi, l=145), deg=TRUE)

{
# returns a 'Milankovitch graph' of incoming solar irradiance
CONVERT=ifelse(deg, 180/pi, 1)

long = long[-length(long)]

M <- outer(long,lat,Vectorize(Insol,c("long","lat")),orbit=orbit,  S0=S0)
attr(M,  "class") <- 'Milankovitch'
attr(M,  "lat") <- lat * CONVERT
attr(M,  "long") <- long * CONVERT
attr(M,  "deg") <- deg
attr(M,  "polar_nights") <- .polar_night_curves(orbit)
M
}



#' plot Milankovitch graph
#' 
#' plot Milankovitch object
#' 
#' 
#' @param x Milankovitch object
#' @param months if true : x-axis of the plot indicates months conventionnally
#' defined with the true solar longitude; x-axis is simply the true solar
#' longitude otherwise
#' @param polar_night if true : the polar night zone will be hashed
#' @param col main color for positive contours (negative contours are in red)
#' @param ... Other arguments passed to plotting function
#' @param plot_function function used to plot the matrix. Typically
#' \code{contour} or \code{image} but may also be \code{image.plot} if using
#' the \code{fields} package.
#' @author Michel Crucifix, U. catholique de Louvain, Belgium.
#' @keywords misc
#' @importFrom graphics axis contour polygon
#' @export
plot.Milankovitch <- function(x, months=TRUE, polar_night=TRUE, plot_function=contour, col="black",...)
{
   long = attr(x, "long")
   if( ! attr(x,  "deg")) 
    {stop ("plot Milankovitch only when longitude in degrees")
    }
if (months)
 { 
   Col  = c(which(long >= (360-80)) , which(long < (360-80)))
   Col = c(Col, Col[1])
   Long = c(long, long[1]+360)
   MM = x[ Col,]
   if (length(which(x < -1))) ## this is a diff plot
   {
     levels <- pretty(range(x), 8)
     cols <- rep("red", length( levels) )
     cols[which(levels < 0)] <- "blue"
     cols[which(levels == 0)] <- "black"
     plot_function(Long, attr(x, "lat"), MM, axes=FALSE,xlab="Month",ylab="Latitude",xaxs="i",yaxs="i",col=cols, levels=levels,...)


     polar_night <- FALSE
   } else {
   plot_function(Long, attr(x, "lat"), MM, axes=FALSE,xlab="Month",ylab="Latitude",xaxs="i",yaxs="i",col=col,...)
   }
    if(!exists("legend.only"))  {
   axis (1, at=seq(0,12)*30, labels=rep("",13))
   axis (3, at=seq(0,12)*30, labels=rep("",13))
   axis (1, at=seq(0,11)*30+15, labels=c("J","F","M","A","M","J","J","A","S","O","N","D"), tick=FALSE)}
   if (polar_night)
   {
   polar_night_curves <- attr(x, ".polar_nights")
   if (!is.null (polar_night_curves))
   {
   lapply(polar_night_curves, function(p) {
      polygon(p[[1]], p[[2]], density=12, col=col)
      polygon(p[[1]], p[[2]], density=12, angle=-45, border=NA, col=col) })
   }
   }
 }
 else { 
 plot_function(attr(x, "long"), attr(x, "lat"), x, axes=FALSE, 
      xlab='True solar Longitude',ylab='Latitude',xaxs='i',yaxs='i',...)
   if(!exists("legend.only")) {
       axis (1, at=seq(0,360,30))
       axis (3, at=seq(0,360,30))
    }
}
 if(!exists("legend.only")) {
 axis(2, at = seq(-90,90,30), labels = c('90S','60S', '30S','Eq.','30N','60N','90N'))
 axis(4, at = seq(-90,90,30), labels = rep('',7))
 }
}

#' Annual Mean insolation
#'
#' Compute annual mean insolation for a given orbit and solar constant
#' @param  orbit Output from a solution, such as \code{ber78}, \code{ber90} or 
#' \code{la04}
#' @param  S0   Total solar irradiance
#' @param  lats  list of latitudes at which annual mean is computed
#' @param  degree true if latitudes are provided in degrees
#' @return an object of class "AnnualMean" with annual mean insolations
#' @export AnnualMean
AnnualMean <- function(orbit, S0 = 1365, lats = seq(-90, 90, 1), degree = TRUE) {
  if (degree) {
    conv <- pi/180; rconv <- 1} 
  else {
    conv <- 1. ; rconv <- 180/pi}
  # colats <- cos(lats)
  A <- as.numeric(sapply(conv*lats, function(l) Insol_l1l2(orbit, lat = l, avg = TRUE)))
  attr(A,  "class") <- "AnnualMean"
  attr(A,  "lats")  <- lats * rconv
  return(A)
}


#' @importFrom graphics lines
#' @export 
plot.AnnualMean <- function(x, vertical = TRUE, yaxs='i',add=FALSE,...) {
  lats <- attr(x, "lats")
  if (!add)
  {
  plot(as.numeric(x), lats, type = "l", axes = FALSE, xlab = "Annual Mean [W/m2]", ylab = "Latitude", 
       yaxs = yaxs, ...)
  axis(1)
  axis(2, at = seq(-90, 90, 30), labels = c("90S", "60S", "30S", "Eq.", "30N", "60N", "90N"))
  } else  {
  lines(as.numeric(x), lats, type = "l",...)
  }
}

#' Begin and end of the polar night
#' 
#' Provides the true solar longitude (in degrees) of the beginning and end of
#' the polar night for a given latitude and orbit
#' 
#' 
#' @param lat latitude
#' @param orbit Output from a solution, such as \code{ber78}, \code{ber90} or
#' \code{la04}
#' @return Either a message about the absence of polar night (for specified
#' reasons), or the true solar longitude, in degrees, of the beginning and end
#' of the polar night.
#' @author Michel Crucifix, U. catholique de Louvain, Belgium.
#' @references any standard text book of spherical astronomy
#' @examples
#' 
#' 
#' current_orbit <- la04(0)
#' 
#' # polar night at the equator ? 
#' polar_night (0, current_orbit)
#' 
#' # polar night at 80 N ? 
#' polar_night (80*pi/180, current_orbit)
#' 
#' # polar nights expressed as day of year 
#' l2day(current_orbit, polar_night (80*pi/180, current_orbit))
#' 
#' 
#' @export polar_night

polar_night <- function(lat, orbit) {
  if (is.null(orbit["eps"])) stop("provided orbit does not include obliquity")
  eps <- orbit["eps"]
  if (eps == 0) return("no polar night/day on a zero-obliquity planet")
  lam <- asin(cos(abs(lat)) / sin(eps))
  if (abs(lat) < abs(pi / 2 - eps))
    return("no polar night / day at this latitude")
  pn <- c(2 * pi - lam, lam) * 360 / (2 * pi)
  # if latitude and obliquity of opposite sign...
  if (lat * eps < 0)  pn <- rev(pn)
  return(pn)
}
     



.polar_night_curves <- function(orbit) {
   pi2 <- pi / 2
   eps <- orbit['eps']
   if (eps < 0.001) return(NULL)  # no polar night
   lam <- function(phi) asin ( - cos (phi) / sin (eps) ) 
   sections <- list( seq(pi/2, pi/2 - eps, -0.005),
                     seq(-pi/2, -pi/2+eps-0.0001, 0.005), 
                     seq(-pi/2+eps-0.0001, -pi/2, -0.005))

   curves <- list( 
      list(c(sections[[1]]),
             c(pi-lam(sections[[1]]))),
      list(c(sections[[1]]),
             (c(lam(sections[[1]])))),
      list(c(sections[[2]],sections[[3]]),
            - c(lam(sections[[2]]), pi-lam(sections[[3]]))))
   curves_degree <- lapply ( curves,  function(l) {
             list(l[[1]] * 180 / pi,
             (l[[2]] * 180 / pi + 360 + 80 ) %% 360 )} ) 

   i <- which(curves_degree [[2]][[2]] > 180 )
   # will not work for exotic obliquities; this needs more work
   if (length(i)) {
   polygon_degree <- list(
     list ( 
                    c( curves_degree[[1]][[2]], 
                      curves_degree[[2]][[2]][(i):length(curves_degree[[2]][[2]])],
                      360, 
                      360, 
                      curves_degree[[1]][[2]][1] ),
                    c( curves_degree[[1]][[1]], 
                      curves_degree[[2]][[1]][(i):length(curves_degree[[2]][[1]])],
                      curves_degree[[2]][[1]][length(curves_degree[[2]][[1]])],
                      90 , 90 ))  ,
     list ( 
                    c( curves_degree[[2]][[2]][1:(i-1)], 0 , 0, 
                      curves_degree[[2]][[2]][1] ),
                    c( curves_degree[[2]][[1]][1:(i-1)], 
                       curves_degree[[2]][[1]][(i-1)], 
                      90 , 90 ))  ,
     list ( 
                    c( curves_degree[[3]][[2]][1] , curves_degree[[3]][[2]], 
                      curves_degree[[3]][[2]][length(curves_degree[[3]][[2]])],  
                      curves_degree[[3]][[2]][1] ),
                    c( -90, curves_degree[[3]][[1]], -90 , -90 ) ))
   } else  { 
   polygon_degree <- list(
     list ( 
                    c( curves_degree[[1]][[2]], 
                      360, 
                      360, 
                      curves_degree[[1]][[2]][1] ),
                    c( curves_degree[[1]][[1]], 
                      curves_degree[[1]][[1]][length(curves_degree[[1]][[1]])],
                      90 , 90 ))  ,
     list ( 
                    c( curves_degree[[2]][[2]], 0 , 0, 
                      curves_degree[[2]][[2]][1] ),
                    c( curves_degree[[2]][[1]], 
                       curves_degree[[2]][[1]][length(curves_degree[[2]][[1]])], 
                      90 , 90 ))  ,
     list ( 
                    c( curves_degree[[3]][[2]][1] , curves_degree[[3]][[2]], 
                      curves_degree[[3]][[2]][length(curves_degree[[3]][[2]])],  
                      curves_degree[[3]][[2]][1] ),
                    c( -90, curves_degree[[3]][[1]], -90 , -90 ) ))
   }
   polygon_degree
}
