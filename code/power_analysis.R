#### Power analysis of model relating species turnover at BBS routes to environmental change

library(tidyverse)
library(sf)
library(purrr)
library(pwr)
library(ggplot2)

### Read in data

## BBS

routes <- read.csv("/Volumes/hurlbertlab/Databases/BBS/2017/bbs_routes_20170712.csv")
weather <- read.csv("/Volumes/hurlbertlab/Databases/BBS/2017/bbs_weather_20170712.csv")
counts <- read.csv("/Volumes/hurlbertlab/Databases/BBS/2017/bbs_counts_20170712.csv")
species <- read.csv("/Volumes/hurlbertlab/Databases/BBS/2017/bbs_species_20170712.csv")

routes$stateroute <- routes$statenum*1000 + routes$route
weather$stateroute <-weather$statenum*1000 + weather$route
RT1 <- subset(weather, runtype == 1, select = c("stateroute", "year"))
RT1.routes <- merge(RT1, routes[ , c("countrynum", "statenum", "stateroute", "latitude", "longitude","bcr")], by = "stateroute", all.x = TRUE)
counts$stateroute <- counts$statenum*1000 + counts$route

# Diurnal land birds, no birds of prey
species_list <- species %>%  
  filter(aou > 2880) %>%
  filter(aou < 3650 | aou > 3810) %>%
  filter(aou < 3900 | aou > 3910) %>%
  filter(aou < 4160 | aou > 4210) %>%
  filter(aou != 7010) %>%
  filter(sporder != "Accipitriformes", 
         sporder != "Falconiformes", 
         sporder != "Anseriformes",
         sporder != "Cathartiformes")

counts.subs <- counts %>%
  inner_join(RT1.routes, by = c("statenum", "stateroute", "year")) %>% # Remove RT = 0 route-years
  filter(aou %in% species_list$aou) %>% # Only diurnal land bird species, no birds of prey
  filter(rpid == 101) %>%
  filter(year >= 1990, year < 2017)

## BCRs

bcrs <- read_sf("/Volumes/hurlbertlab/DiCecco/bcr_terrestrial_shape/BCR_terrestrial_master.shp") %>%
  filter(COUNTRY == "USA" | COUNTRY == "CANADA") %>%
  st_crop(c(xmin = -178, ymin = 18.9, xmax = -53, ymax = 60)) %>%
  filter(PROVINCE_S != "ALASKA" & PROVINCE_S != "HAWAIIAN ISLANDS" & PROVINCE_S != "NUNAVUT" & PROVINCE_S != "NORTHWEST TERRITORIES" & PROVINCE_S != "YUKON") %>%
  filter(WATER == 3)

## ENV change at BBS routes

env_change <- read.csv("data/bbs_route_env_change.csv", stringsAsFactors = F) %>%
  rename(cover_change = "deltaProp",
         edge_change = "deltaED")

## ENV change at BCRs

bcr_change <- env_change %>%
  left_join(dplyr::select(routes, stateroute, longitude, latitude)) %>%
  st_as_sf(coords = c("longitude", "latitude")) %>%
  st_set_crs(st_crs(bcrs)) %>%
  st_intersection(bcrs) %>%
  group_by(BCR) %>%
  summarize(mean_tmax = mean(tmax, na.rm = T),
            mean_tmin = mean(tmin, na.rm = T),
            mean_cover = mean(cover_change, na.rm = T),
            mean_edge = mean(edge_change, na.rm = T)) %>%
  st_set_geometry(NULL)

### Format data - species turnover at 2 scales - route, BCR

temporal_turnover <- function(specieslist1, specieslist2) {
  gamma <- length(unique(c(specieslist1, specieslist2)))
  intersec <- length(specieslist1[specieslist1 %in% specieslist2])
  jaccard <- (gamma - intersec)/gamma
  success <- gamma - intersec
  fail <- intersec
  res <- data.frame(jaccard = jaccard, success = success, fail = fail)
  return(res)
}

route_turnover <- counts.subs %>%
    filter(stateroute %in% env_change$stateroute) %>%
    group_by(stateroute) %>%
    nest() %>%
    mutate(jaccard = purrr::map(data, ~{
      df <- .
      
      df_early <- df %>%
        filter(year <= 1994) %>%
        group_by(aou) %>%
        filter(n_distinct(year) > 1) %>%
        ungroup()
      sp_early <- unique(df_early$aou)
      
      df_late <- df %>%
        filter(year >= 2013) %>%
        filter(n_distinct(year) > 1) %>%
        ungroup()
      sp_late <- unique(df_late$aou)
      
      temporal_turnover(sp_early, sp_late)
    })) %>%
    dplyr::select(-data) %>%
  unnest(cols = c(jaccard)) %>%
  left_join(env_change)

bcr_turnover <- counts.subs %>%
  filter(stateroute %in% env_change$stateroute) %>%
  group_by(bcr) %>%
  nest() %>%
  mutate(jaccard = purrr::map(data, ~{
    df <- .
    
    df_early <- df %>%
      filter(year <= 1994) %>%
      group_by(aou) %>%
      filter(n_distinct(year) > 1) %>%
      ungroup()
    sp_early <- unique(df_early$aou)
    
    df_late <- df %>%
      filter(year >= 2013) %>%
      filter(n_distinct(year) > 1) %>%
      ungroup()
    sp_late <- unique(df_late$aou)
    
    temporal_turnover(sp_early, sp_late)
  })) %>%
  dplyr::select(-data) %>%
  unnest(cols = c(jaccard)) %>%
  left_join(bcr_change, by = c("bcr" = "BCR"))
  
### Build model for power analysis: species turnover at a route over time ~ env change at 3 scales

e <- min(route_turnover$jaccard, na.rm = T)/2

logit <- function(x) {
  logit <- log(x/(1 - x))
  return(logit)
}

logit_e <- function(x) {
  logit <- log((x- e)/(1 - (x - e)))
  return(logit)
}

route_turnover_model <- route_turnover %>%
  filter(!is.na(jaccard)) %>%
  mutate(logit_jac = case_when(
                               jaccard == 1 ~ logit_e(jaccard),
                               jaccard < 1 ~ logit(jaccard)))

e <- min(bcr_turnover$jaccard, na.rm = T)/2

bcr_turnover_model <- bcr_turnover %>%
  filter(!is.na(jaccard)) %>%
  mutate(logit_jac = case_when(
    jaccard == 1 ~ logit_e(jaccard),
    jaccard < 1 ~ logit(jaccard)))

route_turnover_mod <- lm(logit_jac ~ tmax + tmin + cover_change + edge_change, data = route_turnover_model)
# route_turnover_mod <- glm(cbind(success,fail) ~ tmax + tmin + cover_change + edge_change, data = route_turnover, family = binomial(link = "logit"))

bcr_turnover_mod <- lm(logit_jac ~ mean_tmax + mean_tmin + mean_cover + mean_edge, data = bcr_turnover_model)
# bcr_turnover_mod <- glm(cbind(success,fail) ~ mean_tmax + mean_tmin + mean_cover + mean_edge, data = bcr_turnover, family = binomial(link = "logit"))

### Power analysis of current model

power_analysis <- function(model, alpha) {
  sum <- summary(model)
  rsq <- sum$r.squared
  f2 <- rsq/(1 - rsq)
  numdf <- sum$fstatistic[[2]]
  dendf <- sum$fstatistic[[3]]
  
  pwr.f2.test(numdf, dendf, f2, alpha)
}

power_analysis_output <- function(model, pvals) {
  map_df(pvals, ~{
    p <- .
    pwr <- power_analysis(model, p)
    flatten_dfc(pwr)
  })
}

pvals <- seq(0.0001, 0.2, by = 0.005)

route_power <- power_analysis_output(route_turnover_mod, pvals) %>%
  mutate(mod = "Stateroute")
bcr_power <- power_analysis_output(bcr_turnover_mod, pvals) %>%
  mutate(mod = "BCR")

power_results <- bind_rows(route_power, bcr_power)

theme_set(theme_classic(base_size = 15))

ggplot(power_results, aes(x = sig.level, y = power, color = mod)) + geom_line(cex = 1) +
  labs(x = "Significance level", y = "Power", color = "Scale of model") +
  scale_color_manual(values = c("springgreen3", "skyblue3"))
ggsave("figures/existing_sample_power_analysis.pdf")

### Power analysis: artificially decrease number of routes, how does power change

# What would power be for whole US model if we resample from routes to give same density in east as western

na_area_west <- bcrs %>%
  st_crop(c(xmin = -178, ymin = 18.9, xmax = -100, ymax = 60)) %>%
  st_union() %>%
  st_area()

na_area_east <- bcrs %>%
  st_crop(c(xmin = -100, ymin = 18.9, xmax = -50, ymax = 60)) %>%
  st_union() %>%
  st_area()

regional_power <- route_turnover_model %>%
  left_join(dplyr::select(routes, stateroute, latitude, longitude)) %>%
  mutate(region = case_when(longitude < -100 ~ "West",
                            longitude > -100 ~ "East")) %>%
  group_by(region) %>%
  summarize(nRoutes = n_distinct(stateroute))

west_density <- regional_power$nRoutes[[2]]/na_area_west
east_new_nRoutes <- as.numeric(round(west_density*na_area_east))

eastern_routes <- route_turnover_model %>%
  left_join(dplyr::select(routes, stateroute, latitude, longitude)) %>%
  mutate(region = case_when(longitude < -100 ~ "West",
                            longitude > -100 ~ "East")) %>%
  filter(region == "East")

new_eastern_routes <- eastern_routes[sample(nrow(eastern_routes), east_new_nRoutes), ]

subsampled_bbs <- route_turnover_model %>%
  left_join(dplyr::select(routes, stateroute, latitude, longitude)) %>%
  mutate(region = case_when(longitude < -100 ~ "West",
                            longitude > -100 ~ "East")) %>%
  filter(region == "West") %>%
  bind_rows(new_eastern_routes)

subsampled_mod <- lm(logit_jac ~ tmax + tmin + cover_change + edge_change, data = subsampled_bbs)

subsampled_power <- power_analysis_output(subsampled_mod, pvals) %>%
  mutate(mod = "Western US density") %>%
  bind_rows(filter(route_power, mod == "Stateroute"))

ggplot(subsampled_power, aes(x = sig.level, y = power, color = mod)) + geom_line(cex = 1) +
  labs(x = "Significance level", y = "Power", color = "Sample of model") +
  scale_color_manual(values = c("springgreen3", "skyblue3"))
ggsave("figures/power_subsampled_routes.pdf")
