\name{baseball}

\alias{baseball}

\docType{data}

\title{Baseball Hitter's Data}

\description{
The data are for 322 Major Leaque Baseball regular and substitute hitters in 1986.
}

\usage{data(baseball)}

\format{
  A data frame with 322 observations on the following 22 variables.
  \describe{
    \item{\code{Name}}{The hitter/player's name}
    \item{\code{League}}{Player's league (American/National) at the beginning of 1987}
    \item{\code{Team}}{Player's team at the beginning of 1987}
    \item{\code{Position}}{Player's position in 1986: \code{1B} (first
      base), \code{2B} (second base), \code{3B} (third base), \code{C }
      (catcher),  \code{OF} (outfield), \code{DH} (designated hitter),
      \code{SS} (short stop), \code{UT} (utility).}
    \item{\code{Atbat}}{Number of times at bat in 1986}
    \item{\code{Hits}}{Number of hits in 1986}
    \item{\code{Homer}}{Number of home runs in 1986}
    \item{\code{Runs}}{Number of runs in 1986}
    \item{\code{RBI}}{Runs batted in during 1986}
    \item{\code{Walks}}{Number of walks in 1986}
    \item{\code{Years}}{Number of years in the major leagues}
    \item{\code{Atbatc}}{Number of times at bat in his career}
    \item{\code{Hitsc}}{Number of hits in career}
    \item{\code{Homerc}}{Number of home runs in career}
    \item{\code{Runsc}}{Number of runs in career}
    \item{\code{RBIc}}{Number of Runs Batted In in career}
    \item{\code{Walksc}}{Number of walks in career}
    \item{\code{Putouts}}{Number of putouts in 1986}
    \item{\code{Assists}}{Number of assists in 1986}
    \item{\code{Errors}}{Number of errors in 1986}
    \item{\code{Salary}}{Annual salary (in thousands) on opening day 1987}
    \item{\code{logSal}}{Log of salary}
  }
}

\details{
  The levels of the player's positions have been collapsed to fewer
  levels for a simpler analysis.  See the original data for the full
  list of positions.
}

\source{
  The version of the data used to create this data was found at
  \url{http://euclid.psych.yorku.ca/ftp/sas/sssg/data/baseball.sas}.

  Michael Friendly analyzed the data in: \bold{Corrgrams: Exploratory
  Displays for Correlation Matrices}, \emph{The American Statistician},
  Nov 2002, Vol 56.
}

\references{
  The data was originally published for the 1988 ASA Statistical Graphics and
  Computing Data Exposition:
  \url{http://lib.stat.cmu.edu/data-expo/1988.html}.

  The salary data were taken from Sports Illustrated, April 20, 1987.
  The salary of any player not included in that article is listed as an
  NA.  The 1986 and career statistics were taken from The 1987 Baseball
  Encyclopedia Update published by Collier Books, Macmillan Publishing
  Company, New York.
}

\examples{

vars2 <- c("Assists","Atbat","Errors","Hits","Homer","logSal",
           "Putouts","RBI","Runs","Walks","Years")
corrgram(baseball[,vars2],
         lower.panel=panel.shade, upper.panel=panel.pie)

}
