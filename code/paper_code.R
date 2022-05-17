## ---
## author: "[Christophe BONALDI](Christophe.BONALDI@santepuliquefrance.fr)"
## date:  2020-11-16
## source:
## Last modification: 2022-03-09
## R version 4.1.1 (2021-08-10)
## ---
##  Description: 
##  R script for the paper : "Near-real time monitoring of the reproductive 
##  number of COVID-19 in France: estimates compared from 3 datasets"
## ---

## Useful library ----

require( dplyr )
require( tidyr )
require( ggplot2 )
require ( stringr )
require( purrr )

## Input data : ----- 
#'  a data.frame necessarily containing the number
#'  of cases by day and location identity (even if single)
#'  Data frame dt_cas : number of positive tests, emergency 
#'  admissions and hospitalizations for Covid-19 in
#'  Metropolitan France between May 13, 2020 and March 12, 2022


rm( list=ls() )
dt_cas <- readRDS(file = "./data/ds_mt.rds")
str(dt_cas)


#' source useful functions 
source(file = "./fun/Nrt_R_eff_fun.R", encoding = "utf-8")


## Raw and smoothed epidemic curve ----

#' Use function fun_smooth with data
#' according to each information systems
#' Smoothed number of cases are in ss_name_of_SI column 

formals(fun_smooth)
body(fun_smooth)

cas_sys <- c("pos_SIDEP", "pass_OSCOUR", "hosp_SIVIC")

for (i in seq_along(cas_sys) ) {
nm_sys <- unlist(strsplit(cas_sys[i],"_"))[2]
dt_cas[ , paste0("ss_",nm_sys)] <-
  fun_smooth(dt = dt_cas, cases = cas_sys[i])
}


## Plot of the epidemic curves ----
Sys.setlocale("LC_TIME", "English")

#' Transform the data set in format long :
#' raw and smoothed numbers of cases are given in the 
#' column 'cases', 'grp' contain the levels "smoothed" or
#' "raw". These column's names are necessary for 
#' epi_curve function. 

dt_long  <- pivot_longer(
  select(dt_cas, -contains("name")),
  !d,
  names_to = "dat",
  values_to = "cases"
  ) %>%
  mutate(
    grp = ifelse(str_sub(dat,1,2) == "ss", "smoothed", "raw"),
    si = str_split(dat,"_", simplify = T)[,2]
    )

#' A little function to plot lock-down period
lock_down <- function (p, dat1, dat2) {
  #@ dat1/2 : "AAAA-MM-JJ"
  annotate(geom = "rect",
           xmin = as.Date(dat1), 
           xmax = as.Date(dat2),
           ymin = 0, ymax = Inf,
           fill = "blue",
           alpha = 1/5)
}

args_1 <- list(dat1 = "2020-10-29", dat2 = "2020-12-14") #2nd lockdown
args_2 <- list(dat1 = "2021-04-03", dat2 = "2021-05-02") #3rd lockdown


#' Plot of raw and epidemic curves
#' Exemple With SI-DEP numbers only
epi_curve( dt = dplyr::filter(dt_long, si == "SIDEP"))


#' The code to view the 3 data sets, with a logarithmic scale for the y-axis 
#' for a better view and the lockdown periods are added.  
X11(width = 6, height = 7.5, pointsize = 9)
epi_curve( dt = dt_long ) +
  scale_y_continuous(trans='log10', labels = scales::comma) +
  facet_grid(
    factor(si, levels = c("SIDEP","OSCOUR","SIVIC")) ~ .,
    scales = "free_y") +
  do.call(lock_down, args_1) +
  do.call(lock_down, args_2)

#' Warning messages just indicate that the infinite value for the ymax argument
#' (lockdown function above) is transformed into an infinite value by the 
#' logarithm; which is the expected effect! 


## Rt estimates, Metropolitan France ----

#' Require EpiEstim package installed

list_var <- c("hosp_SIVIC", "pass_OSCOUR", "pos_SIDEP")

res3s <- NULL
for(i in seq_along(list_var)) {
  out <- compute_R (dat = dt_cas, dday = "d", var.zone = "reg_name",
                    zone = "metropole",
                    var.cases = list_var[i],
                    first_d = "2020-05-20",
                    window = 7, SI_mean = 7, SI_sd = 5.2)
  res3s <- rbind(res3s, out)
}

tail(res3s)

#' mutate var column in proper character format
dt_rt <- filter(res3s, dates > as.Date("2020-05-26"))  %>% 
  mutate(var = factor(
    var,
    levels = c("pos_SIDEP", "pass_OSCOUR", "hosp_SIVIC"),
    labels = c("SIDEP", "OSCOUR", "SIVIC")))

X11(width = 7, height = 4, pointsize = 9)
fun_plot_R(dt_rt, colour =  c('blue', 'orange', 'red'))


# Supplementary analysis : Infra-national level and Overseas territories ----

#' Load data :

dt_rg <- readRDS(file = "./data/ds_rg.rds")

#' Long format
dt_rglg  <- pivot_longer(
  select(dt_rg, -contains("name")),
  cols = hosp_SIVIC:pos_SIDEP,
  names_to = "dat",
  values_to = "cases"
) %>%
  mutate(
    grp = "raw",
    si = str_split(dat,"_", simplify = T)[,2]
  ) %>% 
  filter(dat != "rea_SIVIC") %>% 
  filter(!(region  == 2 &  si == "OSCOUR"))


#' Smooth the time series of number by area ('region') and dataset ('si')  
dt_nest <- dt_rglg %>%
  group_by(region, si) %>% 
  nest() %>% 
  mutate( ss = purrr::map(data, ~ {fun_smooth(dt = ., cases = "cases")}) ) %>% 
  unnest(cols = c(data, ss)) %>% 
  select ( -cases ) %>% 
  mutate ( grp = "smoothed" ) %>% 
  rename( cases = ss )

dt_rglg <- bind_rows(dt_rglg, dt_nest)

#' Epidemic curves 

# The list of ID for area is in the dt_rg data frame, column 'region' :

dt_rg %>% select(region, name_long) %>% distinct()

#' Metropolitan territory

gm <- dt_rglg %>%
  filter(!(region %in% c(1, 2, 3, 4, 6))) %>% 
  group_by(region) %>% 
  nest() %>% 
  mutate(r_plot =  purrr::map(.x = data, .f = ~ {
    epi_curve(.) +
      scale_y_continuous(trans='log10', labels = scales::comma) +
      facet_grid(
        factor(si, levels = c("SIDEP","OSCOUR","SIVIC")) ~ .,
        scales = "free_y") +
      do.call(lock_down, args_1) +
      do.call(lock_down, args_2) 
      }
    )
  )

 # Exemple : Ile de France ID is 11 ('REGION' column)
 X11(width = 6, height = 7.5, pointsize = 9, title = "Ile de France")
 gm %>% filter(region == 11) %>% pull(r_plot) 
  

#' Overseas territory 

 gom <- dt_rglg %>%
   filter(region %in% c(1, 2, 3, 4, 6)) %>%
   group_by(region) %>%
   nest() %>%
   mutate(r_plot =  map(
     .x = data,
     .f = ~  epi_curve(.) +
       scale_y_continuous(trans='log10', labels = scales::comma) +
       facet_grid(factor(si, levels = c(
         "SIDEP", "OSCOUR", "SIVIC"
       )) ~ .,
       scales = "free_y")
     )
   )
 
 # Exemple : La Réunion island - ID area = 4
 X11(width = 6, height = 7.5, pointsize = 9, title = "La Réunion")
 gom %>% filter(region == 4) %>% pull(r_plot)

 
#' Infra-national Rt estimates 

#' Rt estimates function

fun_comp_R <- function(dt, zone){
  res3s <- NULL
  for(i in seq_along(list_var)) { 
    out <- compute_R (dat = dt, dday = "d", var.zone = "reg_name",
                      zone = zone,
                      var.cases = list_var[i],
                      first_d = "2020-05-20",
                      window = 7, SI_mean = 7, SI_sd = 5.2)
    res3s <- rbind(res3s, out)
  }  
  res3s
}

R_reg <- dt_rg %>%
  group_by(region) %>% 
  nest() %>% 
  mutate(r_res3 =  map(
    .x = data, 
    .f = ~  fun_comp_R(., zone = unique(.x$reg_name))  
  )
)

#' Some warning messages state:
#' "You are estimating R too early in the epidemic
#' to obtain the desired posterior CV". This means 
#' that the estimate of R is based on a small number
#' of observed cases and is not reliable for certain
#' time periods. 

#' Plot of temporal series for estimated R

dt_prepa <- function(res_dt)
{
  res_dt <- filter(res_dt, dates > as.Date("2020-05-26"))  %>% 
    mutate(var = factor(
      var,
      levels = c("pos_SIDEP", "pass_OSCOUR", "hosp_SIVIC"),
      labels = c("SIDEP", "OSCOUR", "SIVIC")))
}

R_reg <- R_reg %>%
  mutate(r_res3 =  map(
    .x = r_res3, 
    .f = ~  dt_prepa(.)  
  )
)


# Metropolitan territories : 


qm <- R_reg %>%
  filter(!(region %in% c(1, 2, 3, 4, 6))) %>% 
  mutate(r_plot =  map(.x = r_res3, .f = ~ {
    fun_plot_R(., colour = c('blue', 'orange', 'red')) +
      do.call(lock_down, args_1) +
      do.call(lock_down, args_2) 
     }
    )
  )


X11(width = 7, height = 4.5, pointsize = 9, title = "Ile de France")
qm %>% filter(region == 11) %>% pull(r_plot)

# Overseas territories : 

qom <- R_reg %>%
  filter(region %in% c(1, 2, 3, 4, 6)) %>% 
  mutate(r_plot =  
    map(.x = r_res3, .f = ~ fun_plot_R(., colour = c('blue', 'orange', 'red')))
    )

X11(width = 7, height = 4.5, pointsize = 9, title = "La Réunion")
qom %>% filter(region == 4) %>% pull(r_plot)


