# Code currently taken from gtseries. Will be merged at some point

#' Develop a spectrum into a time series (generic)
#' @param  arg: input class
#' @param  times: if supplied, times of the decomposition
#' @param  start: if supplied, overrides time and will generate a time series with start and deltat, which must then
#'         be supplied as well
#' @param  deltat : see start. 
#' @note   place holder for type-specific develop functions
#' @export develop
#' @return nothing
develop <- function(M, start=NULL, end=NULL, deltat=NULL, times = NULL, ...){
     UseMethod("develop") 
}

#' Cis function
#'
#' The mathematic function defined by cisx = cos x + i sin x
#'
#' @param x : input
#' @return complex number cos x + i sin x
#' @export 
cis <- function(x) exp(1i*x)

#' Discrete spectrum reconstruction
#' @importFrom stats deltat start
#' @param  M : discreteSpectrum object 
#' @param  times: if supplied, times of the decomposition
#' @param  start: if supplied, overrides time and will generate a time series with start and deltat, which must then
#'         be supplied as well
#' @param  deltat : see start. 
#' @param  sum : TRUE if user wants to sum components %in% the reconstruction
#' @param  trendshift : TRUE if user wants to account for trend and shift encoded in the discreteSpectrum object
#' @param  dfunction is the trigonometrical function. Classically one of 'cos', 'sin', or 'cis'
#' @note   if none if times, start and deltat are supplied, will reconstruct based on the attribute `xdata`
#'         which must then be present. If no `xdata` is availble, return an error. 
#' @return list of reconstructed components if sum=FALSE,  full
#'         reconstructed time series otherwise
#' @method develop discreteSpectrum
#' @export
develop.discreteSpectrum  <- function(M, start = NULL, end = NULL, deltat = NULL, times = NULL,  dfunction = cos, maxfreq = NULL, sum=TRUE, trendshift = TRUE){
 if (!("discreteSpectrum" %in% class(M))) stop ("object is not a discreteSpectrum decomposition")

 timesIsATseries = FALSE
 if (is.ts(times)){
   start = start(times)
   deltat= deltat(times)
   end = start + length(times)*deltat
   timesIsATseries <- TRUE
 }
 if (!is.null(start)){
   if (is.null(deltat) || is.null(end)) stop ("if you supply start, you must also supply deltat and end");
   n <- (end-start) %/% deltat
   times <- ts(start + seq(0, n) * deltat, start=start, deltat=deltat)
   timesIsATseries <-  TRUE
 }

 if (is.null(times)){
   if (is.null(attr(M,"data"))) stop ("if you do not supply any time argument (times, or (start, end, deltat)), then object must have a valid data attribute")
 xdata <- attr(M,"data")
 start <- stats::start(xdata)[1]
 deltat <- stats::deltat(xdata)
 times <- (seq(length(xdata))-1) * deltat + start
 timesIsATseries <- TRUE
 }

 nfreq <- attr(M,"nfreq")
 if (is.null(nfreq)) nfreq <- length(M$Amp)
 if (!is.null(maxfreq)) nfreq <- min(nfreq, maxfreq) 

 if (timesIsATseries){
   reconstructed <- lapply(seq(nfreq), function(i) ts( M$Amp[i] * dfunction(M$Freq[i] * times + M$Phase[i]), start=start, deltat=deltat) )
 } else {
 reconstructed <- sapply(seq(nfreq), function(i) M$Amp[i] * dfunction(M$Freq[i] * times + M$Phase[i]) )
 }

if ( sum ) { 
   shift <- attr(M, "shift"); if (is.null (shift)) shift = 0
   trend <- attr(M, "trend"); if (is.null (trend)) trend = 0
   if (!(trendshift)) { shift = 0; trend = 0 }
   if (timesIsATseries) 
     reconstructed <- Reduce('+', reconstructed)  + trend * times + shift 
   else
     reconstructed <- apply(reconstructed, 1 , sum)  + trend * times + shift
 } else if (timesIsATseries) { 
   reconstructed <- lapply(reconstructed, function(x) ts(x,start=start, deltat=deltat)) 
} 
 return(reconstructed)
}

#' Discrete spectrum class
#' 
#' Class designed to contain amplitudes, frequencies
#' and phases of a trigonometric decomposition
#'
#' @rdname discreteSpectrum
#' @export
as.data.frame.discreteSpectrum <- function(x) {data.frame(Freq=x$Freq, Amp=x$Amp, Phases=x$Phases)}


#' Plot function for descrete spectrum
#'
#' Plot a discrete spectrum object
#'
#' @rdname discreteSpectrum
#' @param a `discreteSpectrum` object, typically the output of a `mfft` call. 
#' @param labels to be set above the frequency peaks. Can be the output of `attributeTone`
#' @param periods if TRUE will add a lower axis with period labels
#' @importFrom graphics mtext axis plot text
#' @export
plot.discreteSpectrum <- function (M,periods=FALSE,labels=NULL,...){
#   O <- order(M$Freq)
  plot(abs(M$Freq), abs(M$Amp),'h',ylab="Amplitudes", xlab="",  ...)
  if (periods) {
    frequencies <- pretty(range(M$Freq/(2*pi)))
    plabels <- as.character(1/frequencies)
    if (0 %in% frequencies) plabels[which(frequencies == 0)] = "∞"
    axis(1, line=3, at=2*pi*frequencies, labels=plabels)
    mtext("Rate", 1, 2)
    mtext("Period", 1, 4)
  } else {
    mtext("Rate", 1, 3)
  }
  # points(abs(M$Freq), abs(M$Amp),'p',...)
  if (!is.null(labels)) {
    yshift <- 0.05*diff(range(M$Amp))
    text(M$Freq, M$Amp, labels, srt=90, , adj=-0.4)
  }
}

#' @rdname discreteSpectrum
#' @export
lines.discreteSpectrum <- function (M,...){
#   O <- order(M$Freq)
  lines(abs(M$Freq), abs(M$Amp),'h',...)
  points(abs(M$Freq), abs(M$Amp),'p',...)
}


#' @rdname discreteSpectrum
#' @export
print.discreteSpectrum <- function (M,...){
  N <- nrow(as.data.frame(M))
  to_print <- seq(min(10,N))
  print.data.frame(cbind(as.data.frame(M)[to_print, ], Period=2*pi/M$Freq[to_print]))
  if (N > 10) print(sprintf("... + %d other rows \n", N-10))
}



