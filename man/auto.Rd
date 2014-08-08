\name{auto}

\alias{auto}

\docType{data}

\title{Statistics of 1979 automobile models}

\description{
  The data give the following statistics for 74 automobiles in the
  1979 model year as sold in the US.
}

\usage{data(auto)}

\format{
  A data frame with 74 observations on the following 14 variables.
  \describe{
    \item{\code{Model}}{Make and model of car.}
    \item{\code{Origin}}{a factor with levels \code{A} \code{E} \code{J}}
    \item{\code{Price}}{Price in dollars.}
    \item{\code{MPG}}{Miles per gallon.}
    \item{\code{Rep78}}{Repair record for 1978 on 1 (worst) to 5 (best) scale.}
    \item{\code{Rep77}}{Repair record for 1978 on 1 to 5 scale.}
    \item{\code{Hroom}}{Headroom in inches.}
    \item{\code{Rseat}}{Rear seat clearance in inches.}
    \item{\code{Trunk}}{Trunk volume in cubic feet.}
    \item{\code{Weight}}{Weight in pounds.}
    \item{\code{Length}}{Length in inches.}
    \item{\code{Turn}}{Turning diameter in feet.}
    \item{\code{Displa}}{Engine displacement in cubic inches.}
    \item{\code{Gratio}}{Gear ratio for high gear.}
  }
}

\source{
  This data frame was created from
  \url{http://euclid.psych.yorku.ca/ftp/sas/sssg/data/auto.sas}.
}

\references{
  Originally published in
  Chambers, Cleveland, Kleiner, and Tukey, \emph{Graphical Methods for Data Analysis},
  1983, pages 352-355.

  The data is from various  sources, primarily \emph{Consumer Reports}, April,
  1979, and the United States government EPA statistics on fuel consumption.
}

\examples{

corrgram(auto[, -c(1:2)])

}


