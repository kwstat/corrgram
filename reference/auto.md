# Statistics of 1979 automobile models

Statistics for 74 automobiles in the 1979 model year as sold in the US.

## Usage

``` r
auto
```

## Format

A data frame with 74 observations on the following 14 variables.

- Model:

  Make and model of car.

- Origin:

  a factor with levels A,E,J

- Price:

  Price in dollars.

- MPG:

  Miles per gallon.

- Rep78:

  Repair record for 1978 on 1 (worst) to 5 (best) scale.

- Rep77:

  Repair record for 1978 on 1 to 5 scale.

- Hroom:

  Headroom in inches.

- Rseat:

  Rear seat clearance in inches.

- Trunk:

  Trunk volume in cubic feet.

- Weight:

  Weight in pounds.

- Length:

  Length in inches.

- Turn:

  Turning diameter in feet.

- Displa:

  Engine displacement in cubic inches.

- Gratio:

  Gear ratio for high gear.

## Source

This data frame was created from
http://euclid.psych.yorku.ca/ftp/sas/sssg/data/auto.sas

## Details

The data is from various sources, primarily *Consumer Reports*, April,
1979, and the United States government EPA statistics on fuel
consumption.

## References

Originally published in Chambers, Cleveland, Kleiner, and Tukey,
*Graphical Methods for Data Analysis*, 1983, pages 352-355.

## Examples

``` r

corrgram(auto[, -c(1:2)])

```
