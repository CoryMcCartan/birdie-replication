# Replication Materials for "Estimating Racial Disparities when Race is Not Observed"

**Paper link:**

### Code

To run the replication, run the R files in this folder in numerical order:
```r
lapply(sort(Sys.glob("*.R")), source)
```

### Data

The geocoded L2 voter file is proprietary and confidential and is not included
in this repository.
However, a version without the proprietary information is included.
Moreover, the estimates from the data are available in
`data-out/nc_ests_party.rda`, which the replication code in `02_nc_fig.R` will 
read from in lieu of full results being available from `01_nc_fit.R`.

The file `data-out/surn_cov.rds` contains the surname classification based on
1930 census data which is used as part of the sensitivity diagnostic.

The IRS tax data studied in the application section, and the analysis code,
is confidential. However, the text output from the code is included as part of
`06_irs_sum.R`, which also reproduces the relevant figures in the paper.
