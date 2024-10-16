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
develop <- function(M, times=NULL, start=NULL, end=NULL, deltat=NULL,...){
     UseMethod("develop") 
}


#' @export cis
cis <- function(x) exp(1i*x)

#' Discrete spectrum reconstruction
#' @importFrom stats deltat start
#' @param  M : discreteSpectrum object 
#' @param  times: if supplied, times of the decomposition
#' @param  start: if supplied, overrides time and will generate a time series with start and deltat, which must then
#'         be supplied as well
#' @param  deltat : see start. 
#' @param  sum : TRUE if user wants to sum components %in% the reconstruction
#' @param  dfunction is the trigonometrical function. Classically one of 'cos', 'sin', or 'cis'
#' @note   if none if times, start and deltat are supplied, will reconstruct based on the attribute `xdata`
#'         which must then be present. If no `xdata` is availble, return an error. 
#' @return list of reconstructed components if sum=FALSE,  full
#'         reconstructed time series otherwise
#' @method develop discreteSpectrum
#' @export
develop.discreteSpectrum  <- function(M, times=NULL, start=NULL, end=NULL, deltat=NULL, dfunction = cos, sum=TRUE){
 if (!("discreteSpectrum" %in% class(M))) stop ("object is not a discreteSpectrum decomposition")

 timesIsATseries = FALSE
 if (!is.null(start)){
   if (is.null(deltat) || is.null(end)) stop ("if you supply start, you must also supply deltat and end");
   n <- (end-start) %*% deltat
   times <- start + seq(0, n) * deltat
   timesIsATseries = TRUE
 }

 if (is.null(times)){
   if (is.null(attr(M,"data"))) stop ("if you do not supply any time argument (times, or (start, end, deltat)), then object must have a valid data attribute")
 xdata <- attr(M,"data")
 start <- stats::start(xdata)
 deltat <- stats::deltat(deltat)
 times <- (seq(length(xdata))-1) * stats::deltat(xdata) + stats::start(xdata)[1]
 timesIsATseries = TRUE
 }

 nfreq <- attr(M,"nfreq")
 if (is.null(nfreq)) nfreq <- length(M$Amp)
 if (timesIsATseries){
   reconstructed <- lapply(seq(nfreq), function(i) ts( M$Amp[i] * dfunction(M$Freq[i] * times + M$Phase[i]), start=start, deltat=deltat) )} 
 else {
   reconstructed <- lapply(seq(nfreq), function(i) M$Amp[i] * dfunction(M$Freq[i] * times + M$Phase[i]))
 }

 if ( sum ) reconstructed <- apply(simplify2array(reconstructed), 1 , sum)
 if (timesIsATseries) reconstructed <- ts(reconstructed, start=start, deltat=deltat)
 return(reconstructed)
}

#' @rdname discreteSpectrum
#' @export
as.data.frame.discreteSpectrum <- function(x) {data.frame(Freq=x$Freq, Amp=x$Amp, Phases=x$Phases)}


#' @rdname discreteSpectrum
#' @param a `discreteSpectrum` object, typically the output of a `mfft` call. 
#' @param labels to be set above the frequency peaks. Can be the output of `attributeTone`
#' @param periods if TRUE will add a lower axis with period labels
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
  print.data.frame(cbind(as.data.frame(M), Period=2*pi/M$Freq))
}



