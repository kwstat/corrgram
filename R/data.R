# data.R

#' Statistics of 1979 automobile models
#' 
#' Statistics for 74 automobiles in the 1979 model year as sold in the US.
#' 
#' The data is from various sources, primarily \emph{Consumer Reports}, April,
#' 1979, and the United States government EPA statistics on fuel consumption.
#' 
#' @format
#' A data frame with 74 observations on the following 14 variables.
#' \describe{ 
#' \item{Model}{Make and model of car.}
#' \item{Origin}{a factor with levels A,E,J}
#' \item{Price}{Price in dollars.} 
#' \item{MPG}{Miles per gallon.} 
#' \item{Rep78}{Repair record for 1978 on 1 (worst) to 5 (best) scale.} 
#' \item{Rep77}{Repair record for 1978 on 1 to 5 scale.}
#' \item{Hroom}{Headroom in inches.} 
#' \item{Rseat}{Rear seat clearance in inches.} 
#' \item{Trunk}{Trunk volume in cubic feet.}
#' \item{Weight}{Weight in pounds.} 
#' \item{Length}{Length in inches.} 
#' \item{Turn}{Turning diameter in feet.}
#' \item{Displa}{Engine displacement in cubic inches.}
#' \item{Gratio}{Gear ratio for high gear.} 
#' }
#' 
#' @source 
#' This data frame was created from
#' http://euclid.psych.yorku.ca/ftp/sas/sssg/data/auto.sas
#' 
#' @references 
#' Originally published in Chambers, Cleveland, Kleiner, and Tukey,
#' \emph{Graphical Methods for Data Analysis}, 1983, pages 352-355.
#' 
#' @examples
#' 
#' corrgram(auto[, -c(1:2)])
#' 
"auto"




#' Baseball Hitter's Data
#' 
#' Data are for 322 Major Leaque Baseball regular and substitute hitters in
#' 1986.
#' 
#' The levels of the player's positions have been collapsed to fewer levels for
#' a simpler analysis.  See the original data for the full list of positions.
#'  
#' The salary data were taken from Sports Illustrated, April 20, 1987.  The
#' salary of any player not included in that article is listed as an NA.  The
#' 1986 and career statistics were taken from The 1987 Baseball Encyclopedia
#' Update published by Collier Books, Macmillan Publishing Company, New York.
#' 
#' @format A data frame with 322 observations on the following 22 variables.
#' \describe{ 
#' \item{Name}{The hitter/player's name}
#' \item{League}{Player's league (American/National) at the beginning
#' of 1987} 
#' \item{Team}{Player's team at the beginning of 1987}
#' \item{Position}{Player's position in 1986: 1B=First base, 
#' 2B=Second base, 3B=Third base, C=Catcher, OF=Outfild, DH=Designated hitter,
#' SS=Short stop, UT=Utility} 
#' \item{Atbat}{Number of times at bat in 1986}
#' \item{Hits}{Number of hits in 1986} 
#' \item{Homer}{Number of home runs in 1986} 
#' \item{Runs}{Number of runs in 1986}
#' \item{RBI}{Runs batted in during 1986} 
#' \item{Walks}{Number of walks in 1986} 
#' \item{Years}{Number of years in the major leagues}
#' \item{Atbatc}{Number of times at bat in his career}
#' \item{Hitsc}{Number of hits in career} 
#' \item{Homerc}{Number of home runs in career} 
#' \item{Runsc}{Number of runs in career}
#' \item{RBIc}{Number of Runs Batted In in career}
#' \item{Walksc}{Number of walks in career}
#' \item{Putouts}{Number of putouts in 1986}
#' \item{Assists}{Number of assists in 1986}
#' \item{Errors}{Number of errors in 1986} 
#' \item{Salary}{Annual salary (in thousands) on opening day 1987} 
#' \item{logSal}{Log of salary} }
#' 
#' @source 
#' The data was originally published for the 1988 ASA Statistical
#' Graphics and Computing Data Exposition:
#' http://lib.stat.cmu.edu/data-expo/1988.html.
#' 
#' The version of the data used to create this data was found at
#' http://euclid.psych.yorku.ca/ftp/sas/sssg/data/baseball.sas
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
"baseball"





#' Voting correlations
#' 
#' Voting correlations
#' 
#' These are the correlations of traits, where each trait is measured for 17 
#' developed countries (Europe, US, Japan, Australia, New Zealand).
#' 
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
"vote"



