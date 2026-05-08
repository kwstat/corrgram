# Baseball Hitter's Data

Data are for 322 Major Leaque Baseball regular and substitute hitters in
1986.

## Usage

``` r
baseball
```

## Format

A data frame with 322 observations on the following 22 variables.

- Name:

  The hitter/player's name

- League:

  Player's league (American/National) at the beginning of 1987

- Team:

  Player's team at the beginning of 1987

- Position:

  Player's position in 1986: 1B=First base, 2B=Second base, 3B=Third
  base, C=Catcher, OF=Outfild, DH=Designated hitter, SS=Short stop,
  UT=Utility

- Atbat:

  Number of times at bat in 1986

- Hits:

  Number of hits in 1986

- Homer:

  Number of home runs in 1986

- Runs:

  Number of runs in 1986

- RBI:

  Runs batted in during 1986

- Walks:

  Number of walks in 1986

- Years:

  Number of years in the major leagues

- Atbatc:

  Number of times at bat in his career

- Hitsc:

  Number of hits in career

- Homerc:

  Number of home runs in career

- Runsc:

  Number of runs in career

- RBIc:

  Number of Runs Batted In in career

- Walksc:

  Number of walks in career

- Putouts:

  Number of putouts in 1986

- Assists:

  Number of assists in 1986

- Errors:

  Number of errors in 1986

- Salary:

  Annual salary (in thousands) on opening day 1987

- logSal:

  Log of salary

## Source

The data was originally published for the 1988 ASA Statistical Graphics
and Computing Data Exposition:
http://lib.stat.cmu.edu/data-expo/1988.html.

The version of the data used to create this data was found at
http://euclid.psych.yorku.ca/ftp/sas/sssg/data/baseball.sas

## Details

The levels of the player's positions have been collapsed to fewer levels
for a simpler analysis. See the original data for the full list of
positions.

The salary data were taken from Sports Illustrated, April 20, 1987. The
salary of any player not included in that article is listed as an NA.
The 1986 and career statistics were taken from The 1987 Baseball
Encyclopedia Update published by Collier Books, Macmillan Publishing
Company, New York.

## References

Michael Friendly (2002). Corrgrams: Exploratory Displays for Correlation
Matrices, *The American Statistician*, Vol 56.

## Examples

``` r

vars2 <- c("Assists","Atbat","Errors","Hits","Homer","logSal",
           "Putouts","RBI","Runs","Walks","Years")
corrgram(baseball[,vars2],
         lower.panel=panel.shade, upper.panel=panel.pie)

```
