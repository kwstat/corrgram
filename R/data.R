

#' Statistics of 1979 automobile models
#' 
#' The data give the following statistics for 74 automobiles in the 1979 model
#' year as sold in the US.
#' 
#' 
#' @name auto
#' @docType data
#' 
#' @format A data frame with 74 observations on the following 14 variables.
#' \describe{ 
#' \item{\code{Model}}{Make and model of car.}
#' \item{\code{Origin}}{a factor with levels A,E,J}
#' \item{\code{Price}}{Price in dollars.} 
#' \item{\code{MPG}}{Miles per gallon.} 
#' \item{\code{Rep78}}{Repair record for 1978 on 1 (worst) to 5 (best) scale.} 
#' \item{\code{Rep77}}{Repair record for 1978 on 1 to 5 scale.}
#' \item{\code{Hroom}}{Headroom in inches.} 
#' \item{\code{Rseat}}{Rear seat clearance in inches.} 
#' \item{\code{Trunk}}{Trunk volume in cubic feet.}
#' \item{\code{Weight}}{Weight in pounds.} 
#' \item{\code{Length}}{Length in inches.} 
#' \item{\code{Turn}}{Turning diameter in feet.}
#' \item{\code{Displa}}{Engine displacement in cubic inches.}
#' \item{\code{Gratio}}{Gear ratio for high gear.} 
#' }
#' 
#' @source This data frame was created from
#' \url{http://euclid.psych.yorku.ca/ftp/sas/sssg/data/auto.sas}.
#' 
#' @references Originally published in Chambers, Cleveland, Kleiner, and Tukey,
#' \emph{Graphical Methods for Data Analysis}, 1983, pages 352-355.
#' 
#' The data is from various sources, primarily \emph{Consumer Reports}, April,
#' 1979, and the United States government EPA statistics on fuel consumption.
#' 
#' @examples
#' 
#' corrgram(auto[, -c(1:2)])
#' 
NULL




#' Baseball Hitter's Data
#' 
#' The data are for 322 Major Leaque Baseball regular and substitute hitters in
#' 1986.
#' 
#' The levels of the player's positions have been collapsed to fewer levels for
#' a simpler analysis.  See the original data for the full list of positions.
#' 
#' @name baseball
#' @docType data
#' 
#' @format A data frame with 322 observations on the following 22 variables.
#' \describe{ 
#' \item{\code{Name}}{The hitter/player's name}
#' \item{\code{League}}{Player's league (American/National) at the beginning
#' of 1987} 
#' \item{\code{Team}}{Player's team at the beginning of 1987}
#' \item{\code{Position}}{Player's position in 1986: 1B=First base, 
#' 2B=Second base, 3B=Third base, C=Catcher, OF=Outfild, DH=Designated hitter,
#' SS=Short stop, UT=Utility} 
#' \item{\code{Atbat}}{Number of times at bat in 1986}
#' \item{\code{Hits}}{Number of hits in 1986} 
#' \item{\code{Homer}}{Number of home runs in 1986} 
#' \item{\code{Runs}}{Number of runs in 1986}
#' \item{\code{RBI}}{Runs batted in during 1986} 
#' \item{\code{Walks}}{Number of walks in 1986} 
#' \item{\code{Years}}{Number of years in the major leagues}
#' \item{\code{Atbatc}}{Number of times at bat in his career}
#' \item{\code{Hitsc}}{Number of hits in career} 
#' \item{\code{Homerc}}{Number of home runs in career} 
#' \item{\code{Runsc}}{Number of runs in career}
#' \item{\code{RBIc}}{Number of Runs Batted In in career}
#' \item{\code{Walksc}}{Number of walks in career}
#' \item{\code{Putouts}}{Number of putouts in 1986}
#' \item{\code{Assists}}{Number of assists in 1986}
#' \item{\code{Errors}}{Number of errors in 1986} 
#' \item{\code{Salary}}{Annual salary (in thousands) on opening day 1987} 
#' \item{\code{logSal}}{Log of salary} }
#' 
#' @source 
#' The data was originally published for the 1988 ASA Statistical
#' Graphics and Computing Data Exposition:
#' http://lib.stat.cmu.edu/data-expo/1988.html.
#' 
#' The version of the data used to create this data was found at
#' \url{http://euclid.psych.yorku.ca/ftp/sas/sssg/data/baseball.sas}.

#' The salary data were taken from Sports Illustrated, April 20, 1987.  The
#' salary of any player not included in that article is listed as an NA.  The
#' 1986 and career statistics were taken from The 1987 Baseball Encyclopedia
#' Update published by Collier Books, Macmillan Publishing Company, New York.
#' 
#' @references 
#' Michael Friendly (2002). Corrgrams: Exploratory Displays for Correlation 
#' Matrices, \emph{The American Statistician}, Vol 56.
#' 
#' @examples
#' 
#' vars2 <- c("Assists","Atbat","Errors","Hits","Homer","logSal",
#'            "Putouts","RBI","Runs","Walks","Years")
#' corrgram(baseball[,vars2],
#'          lower.panel=panel.shade, upper.panel=panel.pie)
#' 
NULL





#' Voting correlations
#' 
#' Voting correlations
#' 
#' These are the correlations of traits, where each trait is measured for 17 
#' developed countries (Europe, US, Japan, Australia, New Zealand).
#' 
#' @name vote
#' @docType data
#' @format A 12x12 matrix.
#' 
#' @source 
#' Torben Iversen and David Soskice (2006). Electoral institutions and
#' the politics of coalitions: Why some democracies redistribute more than
#' others.  \emph{American Political Science Review}, 100, 165-81.  Table A2.
#' 
#' @references 
#' Using Graphs Instead of Tables.
#' http://tables2graphs.com/doku.php?id=03_descriptive_statistics
#' 
#' @examples
#' 
#' corrgram(vote, order=TRUE)
#' 
NULL



