library(tidyverse)
library(furrr)
library(gnnlab)
library(mice)
library(brms)
library(tidybayes)
library(ggrepel)
library(showtext)

font_add_google("Noto Sans", "notosans")


# Load logger dataset ----

survey_dates = read_csv("~/Lab_Data/kawatea/period_info_230829.csv") |> select(location, start_date, end_date)

folder = "~/Lab_Data/kawatea/Oxygen/"
data = tibble(fnames = dir(folder,
                           recursive = T,
                           pattern = "[Cc][Ss][Vv]",
                           full.names = TRUE))

plan(multisession, workers = 16)

data = data |> 
  mutate(data = future_map(fnames, read_onset))
  
data = data %>% ungroup() |> 
  mutate(bnames = basename(fnames)) %>%
  filter(!str_detect(fnames, "mushima|sand|edge|calibration")) %>%
  separate(bnames, c("type", "id", "location", "depth", "survey", "x")) %>%
  select(location, depth, survey, data)

data = data %>% unnest(data) %>%
  mutate(mgl = ifelse(near(mgl,0), NA, mgl)) %>%
  drop_na()

# Check if there is a full dataset (144 points) for each date.

data = data %>%
  mutate(date = as_date(datetime)) %>%
  mutate(datetime = floor_date(datetime, "minutes")) %>%
  group_by(location, depth, date) %>%
  mutate(n = length(date)) %>%
  filter(near(n, 24*6))

tmp =  data %>% 
  mutate(datetime = floor_date(datetime, "day")) %>% 
  group_by(location, datetime) %>% 
  summarise(temperature = mean(temperature), .groups = "drop") |> 
  filter(str_detect(location, "arikawa"),
         datetime >= ymd_h("2018-01-01 0")) 

arikawa = tmp |> 
  mutate(datetime = floor_date(datetime, "month")) |> 
  group_by(datetime) |> 
  summarise(min = mean(temperature)) |> 
  mutate(year = year(datetime),
         month = month(datetime),
         datetime = floor_date(datetime, "year")) |> 
  filter(month %in% c(1,2,3)) |> 
  group_by(datetime) |> 
  summarise(min = mean(min))

# Data from:https://www.data.jma.go.jp/kaiyou/data/shindan/a_1/japan_warm/cfig/warm_area.html?area=C#title
URL = "http://www.data.jma.go.jp/gmd/kaiyou/data/shindan/a_1/japan_warm/cfig/data/areaC_SST.txt"
URL2 = "https://www.data.jma.go.jp/fcd/yoho/typhoon/statistics/generation/generation.csv"
typhoon = read_csv(URL2, locale = locale(encoding = "cp932"))

dset0 = read_csv(URL, skip = 2, n_max = 3)
dset = read_csv(URL, skip = 7, na = c("", "NA", "NoData"))

dset0 = 
  dset0 %>% dplyr::select(Annual:Autumn) %>% 
  slice(1) %>% pivot_longer(everything())

dset  = 
  dset %>% 
  pivot_longer(-Year) %>% 
  full_join(dset0, by = "name") %>% 
  mutate(temperature = value.x + value.y) %>% 
  dplyr::select(-value.y) 


dset2 = tibble(temperature = c(13,14,13.6,13.5,14,14.2),
               year = c(1962,1978,2000,2006,2007,2013),
               location = c("Narao", "Tama-no-ura", "Aoki", "Ojika", "Onihitode","Arikawa"))


dset = 
  dset %>% select(year = Year, name, temperature) %>% 
  pivot_wider(names_from = name, values_from = temperature) 


winter = brm(Winter ~ year, data = dset |> filter(!is.na(Winter)),chains = 4)
summary(winter)

pdata = dset |> expand(year = 1900:2024)

winter_pred = 
  pdata |> 
  add_predicted_draws(winter, ndraws = 4000) %>% 
  ungroup() %>% 
  group_by(year) %>% 
  mean_hdci(.prediction) 

winter_expect = pdata |> 
  add_epred_draws(winter, ndraws = 2000) %>% 
  ungroup() %>% 
  group_by(year) %>% 
  mean_hdci(.epred) 

get_variables(winter)

sumout= spread_draws(winter, b_year) |> 
  mutate(b_year = 100*b_year) |> mean_hdci()
sumout

(0.01*sumout[1] * 124 / 6) * 180 / pi

summary(winter)
library(ggforce)

showtext_auto()
library(ggtext)


set.seed(2023)
LABEL = sprintf("Mean winter temperatures increasing at a rate of **%1.1f (%1.1f to %1.1f)** °C / 100 years (mean and 95%% CI)", sumout[1], sumout[2], sumout[3])

RATIO = 7/16
p1 = ggplot() + 
  geom_point(aes(x = year, y = Winter, color = "Japan Meterological Agency (SST)"), data = dset)  +
  geom_ribbon(aes(x = year, ymin = .lower, ymax = .upper), alpha = 0.3, data = winter_pred) +
  # geom_ribbon(aes(x = year, ymin = .lower, ymax = .upper), alpha = 0.5, data = winter_expect) +
  geom_line(aes(x = year, y = .epred), color = "white", size = 2, data = winter_expect) +
  
  scale_color_viridis_d() +
  geom_point(aes(x = year, y = temperature, color = "Nagasaki region (reports)"), data = dset2, size = 2) +
  
  geom_point(aes(x = year, y = temperature, color = "Arikawa Bay (data)"), size = 3,
             data = arikawa %>% mutate(year = year(datetime), temperature = min)) +
  
  scale_x_continuous("Year", 
                     limits = c(1900, 2024),
                     breaks = c(seq(1900, 2000, by = 20), 2024)) +
  scale_y_continuous("Mean winter water temperature (°C)",
                     limits = c(10, 16),
                     breaks = seq(10, 16, by = 2),
                     labels = \(x) {str_c(x, "°C")}) +
  annotate("richtext", x= 1960.5, y = 14.5, label = LABEL,
           hjust = 0.5, vjust = 0, label.color = "grey30",
           angle = (atan(0.01*sumout[1] * 124 / 6 * RATIO)) * 180 / pi
           ) +
  theme(legend.position = c(1,0),
        legend.justification = c(1,0),
        legend.background = element_blank(),
        legend.direction = "horizontal",
        legend.key = element_rect(fill = "grey80"),
        legend.title = element_blank(), 
        legend.text = element_text(size = 10),
        text = element_text("notosans", size = 16), 
        panel.border = element_rect(color = "grey30", fill = NA),
        panel.background = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        axis.text = element_text(color = "grey30"),
        axis.ticks = element_line(color = "grey30")
        ); p1
pdfname = "2023temperature.pdf"
pngname = str_replace(pdfname, "pdf", "png")

width = 200
height = width * RATIO
save_plot(pdfname = pdfname, plot = p1, width = width, height = height, units = "mm")

