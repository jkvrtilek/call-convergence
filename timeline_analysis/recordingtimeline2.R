# CITE CHERULEAN PALETTE: https://gist.github.com/cherfychow/e9ae890fd16f4c86730748c067feee2b

library(tidyverse)

setwd("/Users/jkvrtilek/Desktop/OSU/PhD/GitHub/call-convergence/timeline_analysis")

# load first contact data from 2016
fc16 <- read_csv("/Users/jkvrtilek/Desktop/OSU/PhD/Calls/first_contact_events2016.csv")

# load recording date and capture site data
metadata <- read_csv("metadata.csv") %>% 
  select("date", "bat_ID", "capture.site") %>% 
  distinct() %>% 
  separate(date, into = c("year","month","day"), remove = FALSE)

# separate out 2016-17
rec_date16 <- metadata %>% 
  filter(year == "2016" | year == "2017") %>% 
  mutate(date = as.Date(date, "%Y-%m-%d"))

# get rid of non-Tole/Las Pavas bats
rec_date16$capture.site <- case_when(rec_date16$capture.site == "tole" ~ "Tole",
                                   rec_date16$capture.site == "gamboa" ~ "Never introduced",
                                   rec_date16$capture.site == "las.pavas" ~ "Las Pavas",
                                   rec_date16$capture.site == "chilibre" ~ "Never introduced",
                                   rec_date16$capture.site == "captive" ~ "Never introduced")
rec_date16$capture.site <- factor(rec_date16$capture.site, levels = c("Tole", "Las Pavas", "Never introduced"))
rec_date16 <- rec_date16 %>% 
  filter(capture.site != "Never introduced")

# get rid of bats with no post-intro recordings
postbats <- rec_date16 %>% 
  filter(date > "2016-09-01") %>% 
  select(bat_ID) %>% 
  distinct() %>% 
  filter(bat_ID != "dos") %>% 
  filter(bat_ID != "ss")
rec_date16 <- rec_date16 %>% 
  subset(bat_ID %in% postbats$bat_ID)

# make df of 2016 first contact dates
cont_date16 <- fc16 %>% 
  mutate(time = as.character(time)) %>% 
  separate(col = time, into = c("date", "time"), sep = " ") %>% 
  distinct(date) %>% 
  filter(!row_number() %in% c(3,4,6)) %>% 
  mutate(date = as.Date(date, "%Y-%m-%d"))
# Feb 7 is a typo for Jul 2
# Jul 7 and Aug 26 carry over from introductions on Jul 6 and Aug 25
# legitimate dates: 2016-07-02, 2016-07-06, 2016-08-25
smallcage <- cont_date16$date[1]
largecage <- cont_date16$date[3]

# make figure
rec_date16 %>% 
  mutate(capture.site = fct_reorder(capture.site, desc(capture.site))) %>% 
  ggplot(aes(x=date, y=bat_ID, group=bat_ID)) +
  facet_wrap(~capture.site, scales = "free_y", ncol = 1) +
  geom_point(aes(color = capture.site)) + 
  geom_line() +
  scale_x_date(date_labels = "%b %Y", breaks = "3 months") +
  scale_color_manual(values = c("#0f2d59","#3eb7c7")) +
  guides(color=guide_legend("Bat Origin")) +
  geom_vline(xintercept = as.numeric(smallcage), color = "red", linetype = "dashed") +
  geom_vline(xintercept = as.numeric(largecage), color = "red") +
  ylab("individual bats") +
  xlab("recording dates") +
  theme_bw() +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        #panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 20),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        strip.text.x = element_text(size = 12))

# make figure with bats reordered and no las pavas
reorder2 <- rec_date16 %>% 
  filter(capture.site == "Tole") %>% 
  mutate(bat_ID = factor(bat_ID, levels = c("cs","dcd","ccs","scs","ccc","rc","r"))) #tole

reorder2 %>% 
  ggplot(aes(x=date, y=bat_ID, group=bat_ID, color = "#331875")) +
  geom_point(size=3) + 
  geom_line(size=1) +
  scale_x_date(date_labels = "%b %Y", breaks = "3 months") +
  scale_color_manual(values = "#331875") +
  geom_vline(xintercept = as.numeric(smallcage), color = "red", linetype = "dashed") +
  geom_vline(xintercept = as.numeric(largecage), color = "red") +
  ylab("individual bats") +
  xlab("recording dates") +
  theme_bw() +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        #panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 20),
        axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24))
