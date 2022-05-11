# select powell center sites based on NEP and river type for eco memory analysis

library(tidyverse)
library(lubridate)

#Powell data
dd <- read_tsv('data/powell_data/site_data.tsv')
dat <- readRDS('data/powell_data/high_quality_daily_metabolism_with_SP_covariates.rds')


dat <- select(dat, site_id = site_name, name = long_name, date, year, DOY,
       GPP_raw, GPP, ER_raw, ER, NEP, temp.water, discharge, LAI, Stream_PAR,
       PAR, width, slope) %>%
    as_tibble()

# Autotrophic sites to start off with:
# nwis_10133800 East Canyon Creek: high elevation desert stream
# nwis_14211010 Clackamas River
# nwis_08447000 Pecos river

sites <- c('nwis_10133800', 'nwis_14211010', 'nwis_08447000')

dat <- dat %>% filter(site_id %in% sites) %>%
    group_by(site_id) %>%
    complete(date = seq(min(date), max(date), by = 'day'))%>%
    mutate(name = name[1],
           site_id = site_id[1],
           width = width[1],
           slope = slope[1],
           year = year(date),
           DOY = as.numeric(format(date, '%j')))

ggplot(dat, aes(date, GPP)) +
    geom_point() +
    geom_point(aes(y = ER)) +
    facet_wrap(.~site_id, scales = 'free')

dat$ecosystem = 'river'
write_csv(dat, 'data/powell_selected_sites.csv')
