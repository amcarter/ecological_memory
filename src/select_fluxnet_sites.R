# Combine fluxnet metadata files,
# select fluxnet sites based on biome type for eco memory analysis

library(tidyverse)

# Fluxnet Data
dd <- readRDS('data/fluxnet_data/fluxnet_site_info_filtered.rds') # From Bernhardt et al 2022
igbp <- read_csv('data/fluxnet_data/fluxnet_site_info_with_IGBP.csv') %>% # from the fluxnet site info page
    select(Site_ID, Elev, IGBP, MATemp_C, MAPrecip_mm)
key <- read_csv('data/fluxnet_data/IGBP_key.csv')

dd <- left_join(dd, igbp, by = 'Site_ID') %>% as_tibble()
write_csv(dd, 'data/fluxnet_site_info.csv')

list_dat <- readRDS('data/fluxnet_data/fluxnet_filtered_metabolism.rds')
dat <- tibble()
for(i in 1:length(list_dat)){
    tmp <- list_dat[[i]]
    tmp$Site_ID = names(list_dat)[i]
    dat <- bind_rows(dat, tmp)
}
dat$U_ID <- dat$U_ID[,1]


# Select Ecosystem Types:
dd %>% mutate(MAPrecip_mm = MAPrecip_mm/100) %>%
    pivot_longer(cols = c('MATemp_C', 'MAPrecip_mm'), names_to = 'variable',
                    names_pattern = '^MA([a-zA-Z_]+$)',
                    values_to = 'value') %>%
    group_by(IGBP, variable) %>%
    summarize(val = mean(value, na.rm = T)) %>%
    ggplot(aes(IGBP, val, fill = variable)) +
    geom_bar(stat = 'identity', position = 'dodge')

# Grassland:
dd %>% filter(IGBP == 'GRA') %>%
    arrange(-nyears)
# CH-Cha: is a 'heavily managed' grassland in Switzerland
# US-SRG: a semi-desert grassland
# NL-Hor: dutch grassland
sites <- c('CH-Cha', 'US-SRG', 'NL-Hor')

#Desiduous: and evergreen forests:
dd %>% filter(IGBP == 'DBF') %>%
    arrange(-nyears)
# US-HA1: Harvard Forest (temperate broadleaf)
dd %>% filter(IGBP == 'ENF') %>%
    arrange(-nyears)
# US-NR1: Niwot Ridge Forest (high elecation evergreen)
sites <- c(sites, 'US-Ha1', 'US-NR1')

#Wetland:
dd %>% filter(IGBP == 'WET') %>%
    arrange(-nyears)
# CZ-wet: sedgegrass marsh land
sites <- c(sites, 'CZ-wet')

#Cropland:
dd %>% filter(IGBP == 'CRO') %>%
    arrange(-nyears)
# FR-Gri: crop rotation site
sites <- c(sites, 'FR-Gri')

dat <- filter(dat, Site_ID %in% sites)
ggplot(dat, aes(DOY, GPP, col = log(Precip))) +
    geom_point() +
    geom_point(aes(y = ER)) +
    facet_wrap(.~Site_ID)
ggplot(dat, aes(Date, GPP, col = Temp)) +
    geom_point() +
    geom_point(aes(y = ER)) +
    facet_wrap(.~Site_ID, scales = 'free', ncol = 2)


dat <- rename(dat, site_id = Site_ID, date = Date, year = Year, temp_C = Temp,
              precip_mm = Precip, shortwave = SW) %>%
    select(-U_ID)

key <- data.frame(IGBP = unique(flux$IGBP),
                  ecosystem = c('grassland', 'wetland', 'cropland',
                                'deciduous forest', 'evergreen forest'))


dat <- select(dd, site_id = Site_ID, name = Name, IGBP) %>%
    right_join(dat, by = 'site_id') %>%
    right_join(key, by = 'IGBP')

write_csv(dat, 'data/fluxnet_selected_sites.csv')

