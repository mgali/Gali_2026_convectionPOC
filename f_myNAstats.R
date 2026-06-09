# Wraper function for a few simple custom stats

nancounts <- function(x) {sum(!(is.na(x) | is.nan(x)))} # Define function to count non-NA and non-NaN values
nansum <- function(x) {sum(x, na.rm = T)}
nanmean <- function(x) {mean(x, na.rm = T)}
nansd <- function(x) {sd(x, na.rm = T)}
nanskewness <- function(x) {skewness(x, na.rm = T)}
nanmedian <- function(x) {median(x, na.rm = T)}
nanq75 <- function(x) {quantile(x, 0.75, na.rm = T)}
na_0 <- function(x) {x[x==0] <- NA; return(x)}
no0mean <- function(x) {mean(x[x!=0], na.rm = T)}
no0sd <- function(x) {sd(x[x!=0], na.rm = T)}
