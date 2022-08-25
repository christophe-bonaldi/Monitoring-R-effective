## ---
## author:"[Juliette PAIREAU](juliette.paireau@pasteur.fr)"
##        "[Anne FOUILLET] (ANNE.FOUILLET@santepubliquefrance.fr)"
##        "[Christophe BONALDI](Christophe.BONALDI@santepubliquefrance.fr)"
## date:  2020-11-16
## Last modification: 2020-
## R version 4.1.1 (2021-08-10)
## ---
##  Description:
## ---
##  Functions and tools
## ---


# Moving average 7 days -----
ma <- function(x, n = 7){stats::filter(x, rep(1 / n, n), sides = 2)}

#  Function to compute R -----
#' @parameters
#'   - dat : a data.frame necessarily containing the number
#'     of cases (@argument : var.cases) by day (@argument : dday)
#'     and by geographical level (@arguments : zone, var.zone)
#'   - dday : date of the day (the default id "d")
#'   - var.zone : character string specifying the column's
#'     name for geographic level
#'   - zone : character string specifying the geographical
#'     level to be used for estimates
#'   - var.cases : character string specifying the column
#'     for incidence time serie (non-negative)
#'   - first_d : character string specifying the first date
#'     to begin estimates into in format "yyyy-mm-dd"
#'   - method : Estimated instantaneous reproduction number
#'     method. The default is Cori with "parametric_si", others
#'     methods not implemented yet
#'   - window : moving time window in day
#'   - SI_mean : average serial interval
#'   - SI_sd : serial sd
#' @authors Juliette Paireau, Christophe Bonaldi
#' @return compute_R return a data.frame contaning the estimated
#'  instantaneous reproduction number for each day, std and quantiles
#' @call EpiEstim::estimate_R, ma
compute_R <- function( dat, dday = "d", var.zone, zone, var.cases,
                       first_d = "2020-03-15",
                       method = "Cori", window, SI_mean, SI_sd){

  require('EpiEstim')
  dat <- as.data.frame(dat)
  tmp <- dat[
    which(dat[, var.zone] == zone & !is.na(dat[, var.cases])),
    c(var.cases, dday)
    ]
  names(tmp)<-c("cases","d")

  if (sum(tmp$cases != 0) == 0){
    return(NULL)
  } else {

    # number of incident cases within the last window
    cas_sum <- sum(tmp[(nrow(tmp)-window+1):nrow(tmp), "cases"], na.rm = T)

    # smooth time-series :
    #' find the smoothing parameter that best reproduces a 7-day moving
    #' average (without loosing the last 3 days)
    smoothed_ma <- ma(tmp$cases)
    cor <- 0
    for (k in seq(0.2, 0.7, by=0.1)){
      smoothed_ss <- smooth.spline(x = tmp$d, y = tmp$cases, spar = k)$y
      cor.tmp <- cor.test(smoothed_ma, smoothed_ss)$estimate
      if (cor.tmp > cor){
        cor <- cor.tmp
        sparam <- k
      }
    }
    tmp$smoothed <- smooth.spline(x = tmp$d, y = tmp$cases, spar = sparam)$y
    tmp$smoothed[tmp$smoothed < 0] <- 0

    if (method == "Cori"){
      incid <- tmp[, c("d", "smoothed")]
      names(incid)<-c("dates","I")
      # beginning estimates at first_d date
      t_start <- seq(which(
                       incid$dates == (as.Date(first_d) - window + 1)),
                       nrow(incid)-(window-1)
                     )
      t_end <- t_start + window - 1

      resR <- estimate_R(
                incid = incid,
                method = "parametric_si",
                config = make_config(
                           list(
                             mean_si = SI_mean,
                             std_si = SI_sd,
                             t_start = t_start,
                             t_end = t_end,
                             mean_prior = 1,
                             std_prior=2)
                           )
                )
      resR$R$dates <- incid$d[resR$R$t_end[1]:resR$R$t_end[nrow(resR$R)]]
      resR$R$var <- var.cases
      resR$R[,var.zone] <- zone
      resR$R$cas_sum <- cas_sum
      resR$R$window <- window
      resR$R$SI <- SI_mean
    } else {
      stop("method not implemented yet")
    }
    return(resR$R)
  }
}


#  Function to smooth epidemic curve -----
#' @parameters
#'   - dt : a data.frame necessarily containing the number
#'     of cases (@argument : cases) by day (column name = "d")
#'   - cases : names of column containing the number of cases per day
#' @authors Juliette Paireau, Christophe Bonaldi
#' @return smoothed number of cases
fun_smooth <- function(dt, cases)
{
  dt <- as.data.frame(dt)
  smoothed_ma <- ma(dt[, cases])
  cor <- 0
  for (k in seq(0.2, 0.7, by = 0.1)){
    smoothed_ss <- smooth.spline(x = dt[, "d"],
                                 y = dt[, cases], spar = k)$y
    cor.tmp <- cor.test(smoothed_ma,  smoothed_ss)$estimate
    if (cor.tmp > cor){
      cor <- cor.tmp
      sparam <- k
    }
  }
  smoothed <- smooth.spline(x = dt[ , "d"], 
                            y = dt[ , cases], spar = sparam)
  smoothed$y[smoothed$y < 0] <- 0
  smoothed$y
}

#' Function to plot raw and smoothed epidemic curve -------------
#' @parameters
#'   - dt : a data.frame necessarily containing the number
#'     of cases (necessary column's name : cases) by day (necessary 
#'     column's name : d) and the indicator for raw or smoothed data 
#'     (necessary column's name : grp)
#' @authors Christophe Bonaldi
#' @return epi_curve return a plot of raw and smoothed number
#' number of cases
epi_curve <- function(dt) 
{
  dt %>% 
    ggplot(aes(x = d, y = cases, group = grp)) +
    geom_line(aes(color = grp, size = grp)) +
    scale_color_manual(values=c('grey40','blue')) +
    scale_size_manual(values=c(.7, 1)) +  
    xlab("") +
    ylab("Daily numbers") +
    scale_x_date(date_breaks = "2 months", date_labels = "%Y %b") +
    theme_light() +
    theme(
      axis.text.x = element_text(angle = 30, hjust = 1, size = 8 ),
      axis.title=element_text(size = 8),
      axis.title.x = element_blank(),
	  legend.position = c(.15, .95),
      legend.title = element_blank()
    ) 
}    


# Function to plot R trajectory ----------------
#' @parameters
#'   - res_dt : a data.frame resulting of a call to compute_R
#'     necessarily containing the R-effective estimates (column name :
#'     `Mean(R)`), their lower and upper confidence limits (column names :
#'     `Quantile.0.025(R)` and `Quantile.0.975(R)`, by date (names : dates)
#'     and by geographical level (zone, var.zone) and dataset (var)
#'   - colour : colors for curves
#' @authors Christophe Bonaldi
#' @return fun_plot_R return a plot of time serie(s) estimates of
#'   instantaneous reproduction number for each names_var componant
fun_plot_R <- function(res_dt, colour)
  {
    bq <- res_dt %>% ggplot(aes(x = dates, y = `Mean(R)`, group = var)) +
    geom_ribbon(
      aes(ymin = `Quantile.0.025(R)`, ymax = `Quantile.0.975(R)`, fill = var,  group = var),
      alpha = 0.3
    ) +
    scale_fill_manual(values = colour) +
    geom_line(aes(color = var)) +
    scale_color_manual(values = colour) +
    xlab("") +
    ylab(as.expression(bquote(R[t]))) +
    scale_x_date(date_breaks = "2 months", date_labels = "%Y %b") +
    geom_hline(yintercept = 1) +
    theme_light() +
    theme(
      axis.text.x = element_text(angle = 30, hjust = 1, size = 8 ),
      axis.title=element_text(size = 11),
      axis.title.x = element_blank(),
      legend.position = c(.15, .85),
      legend.title = element_blank()
    )
  
    bq
}



# Function to produce doubling time from the Reproduction number  -------
#' @parameter
#' Reff : Reproduction number estimates
#' @authors Juliette Paireau
get_doubling_time <- function(Reff) {
  g1 <- 1. / 4.
  g2 <- 1. / 1.
  g3 <- 1. / 3.
  res <- (log(2.) / optim(c(0.1), to_minimize, R = Reff,
                         g1 = g1, g2 = g2, g3 = g3, method = "Brent",
                         lower = 0, upper = 1)$par[1])
  return(ifelse(Reff <= 1, NA, res))
}

# Function to minimize
to_minimize <- function(x, R, g1, g2, g3) {
  c2 <- g3 / (g2 + g3)
  c3 <- g2 / (g2 + g3)
  tmp <- R - (1. + x / g1) * (1. + x / g2) * (1. / (c2 + c3 * g3 / (g3 + x)))
  return(tmp ^ 2)
}
