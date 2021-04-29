library(readxl)
library(plyr)
library(lubridate)
library("scales")

# User modifiable variables
file <- 'genomeMetaData.xlsx'
start_date <- "2020-01-01"
lineage_important <- c("B.1.526", "B.1.526.1", "B.1.525", "P.2", "B.1.1.7", "P.1", "B.1.351", "B.1.427", "B.1.429", "B.1.311",  "B.1", "B.1.2", "B.1.243")

# R Script _____________________________________________________

# Read the excel file
df_original <- read_excel(file)

# Make a new dataframe with only the necessary columns
df1 <- df_original[,c("sample_date", "lineage", "rationale")]

# Remove missing dates from df1
df2 <- na.omit(df1)

# Make the date column a date object
df2$sample_date <- as.Date(df2$sample_date, "%Y-%m-%d")

# Make a dataframe of only the random samples
df3a <- df2[df2$rationale == 'Asymptomatic',]
df3b <- df2[df2$rationale == 'Hospitalized',]
df3c <- df2[df2$rationale == 'Random',]
df3d <- df2[df2$rationale == 'Random; Vaccine Breakthrough, S drop',]
df3e <- df2[df2$rationale == 'Vaccine breakthrough',]
df4 <- rbind(df3a, df3b, df3c)

# Summarize the data
df4$sample_date_week <- floor_date(df4$sample_date, "week") # Adds a column with date for the week only
df5 <- ddply(df4, c("sample_date_week", "lineage"), summarise, n = length(lineage)) # Makes a new dataframe with the counts of each lineage for each of the week categories
df6 <- na.omit(df5) # Remove any data that does not have an associated date

# Make a list of the unique dates and lineages
lineage_unique <- unique(df6[c("lineage")])
date_week_unique <- unique(df6[c('sample_date_week')])

# Make the columns of the graph with the repeating lineages for each date
lineage_repeat <- lineage_unique
for(d in 2:lengths(date_week_unique)){
  lineage_repeat <- Map(c, lineage_repeat, lineage_unique)
}

# Make the column of the graph wit the repeating date, where there are enough dates for each lineage.
date_repeat <- list(sample_date = date_week_unique[1,1])
for(d in 1:lengths(date_week_unique)){
  for(l in 1:lengths(lineage_unique)){
    date_repeat <-Map(c, date_repeat, date_week_unique[d,1])
  }
} 
lengths(date_repeat)

# Remove the initial repeated date from the creation of this new list
date_repeat <- date_repeat$sample_date[-2]
date_repeat <- list(sample_date = date_repeat)
lengths(date_repeat)

# Combine the repeated dates and lineages into a single data frame
df7 <- c(date_repeat, lineage_repeat)
df7 <- data.frame(df7)

# Make a list of zero counts to be temporary placeholders.  
count_zero <- rep(c(0),each=length(df7$sample_date))
count_zero <- list(counts = count_zero)
df8 <- c(date_repeat, lineage_repeat, count_zero)
df8 <- data.frame(df8)

# Loop through the data to determine which dates have counts for each lineage and replace zero placeholders with appropriate counts
short_length <- length(df6$sample_date_week)
long_length <- length(df8$sample_date)
for(s in 1:short_length){
  temp_short_date <- df6[s,1]
  temp_short_lineage <- df6[s,2]
  for(l in 1:long_length){
    temp_long_date <- df8[l,1]
    temp_long_lineage <- df8[l,2]
    if(temp_short_date == temp_long_date && temp_short_lineage == temp_long_lineage){
      df8[l,3] <- df6[s,3] 
    }
  }
}

# Load libraries required for graphing the data
library(ggplot2)
library(dplyr)
library(viridis)
library(hrbrthemes)

# Apply the date cut off to include reads only as old as the start_date specified
df9 <- df8 %>%
  select(sample_date, lineage, counts) %>%
  filter(sample_date >= as.Date(start_date))

# Assign data frame
data <- df9

# Convert variant counts to percentages
data2 <- data  %>%
  group_by(sample_date, lineage) %>%
  summarise(n = sum(counts)) %>%
  mutate(percentage =100 * n / sum(n),
         date=as.Date(sample_date, "%m/%d/%Y")) %>%
  ungroup() %>% as.data.frame()

# Determine which lineages to drop
lineage_unique2 <- c("")
for(lu in 1:lengths(lineage_unique)){
  lineage_unique2 <- c(lineage_unique2, lineage_unique[lu,1])
}
lineage_unique2 <- c(lineage_unique2[2:length(lineage_unique2)])
lineage_drop <- setdiff(lineage_unique2, lineage_important)
#lineage_drop <- c("NA")

# Lump variants into 'Other' category
df9 %>%
  as_tibble() %>%
  mutate(date = as.Date(sample_date, "%Y-%m-%d")) %>%
  # Lumped "Other" Variants chosen at random here
  mutate(Lineage = factor(lineage),
         Lineage = forcats::fct_other(f = Lineage, drop = lineage_drop, other_level = "Other")) %>%
  group_by(date, Lineage) %>%
  summarise(Lineage_count = sum(counts, na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(date) %>%
  mutate(Lineage_perc = 100 * Lineage_count / sum(Lineage_count)) %>%
  ungroup() %>%
  identity() -> data3

# Add specimen totals
data3 %>%
  group_by(date) %>%
  summarise(total_seqs = sum(Lineage_count, na.rm = TRUE)) %>%
  ggplot(data = ., aes(x = date, y = total_seqs)) +
  geom_point() +
  geom_line() +
  theme_ipsum() +
  labs(x = "Date", y = "Genomes Sequenced")

# # Plot lineages using line graph
# p1 <- ggplot(data3,aes(x=date, y=Lineage_perc, fill=Lineage)) +
#   geom_area(alpha=0.6 , size=.5, color="white") +
#   labs(fill = "", x = "Date", y = "Percentage") +
#   theme_ipsum() +
#   ggtitle("Philadelphia Random Sampling SARS-CoV-2 Variants") +
#   scale_x_date(date_labels = "%b %Y", breaks = breaks_pretty(6))
# 
# # Plot sequencing counts using line graph
# p2 <- data3 %>%
#   group_by(date) %>%
#   summarise(total_seqs = sum(Lineage_count, na.rm = TRUE)) %>%
#   ggplot(data = ., aes(x = date, y = total_seqs)) +
#   geom_point() +
#   geom_line() +
#   theme_ipsum() +
#   labs(x = "Date", y = "Sequence Total") +
#   scale_x_date(date_labels = "%b %Y", breaks = breaks_pretty(6))


# Plot lineages using bar graph
library("scales")
p1 <- ggplot(data3, aes(fill=Lineage, y=Lineage_perc, x=date)) +
   geom_bar(position="stack", stat="identity") +
   theme_ipsum() +
   ggtitle("Philadelphia Random Sampling SARS-CoV-2 Variants") +
   labs(x = "Date", y = "Percentage") +
   scale_x_date(date_labels = "%b %Y", breaks = breaks_pretty(6)) + 
   scale_fill_manual(values = c("yellow2", "palegreen", "lightseagreen", "skyblue", "steelblue1", "royalblue", "navy", "purple", "hotpink", "indianred", "red", "orange","honeydew4"))


# Plot sequencing counts using bar graph
p2 <- data3 %>%
  group_by(date) %>%
  summarise(total_seqs = sum(Lineage_count, na.rm = TRUE)) %>%
  ggplot(data = ., aes(x = date, y = total_seqs)) +
  geom_bar(position="stack", stat="identity") +
  theme_ipsum() +
  labs(x = "Date", y = "Genomes Sequenced") +
  scale_x_date(date_labels = "%b %Y", breaks = breaks_pretty(6))

# Combine lineages and sequence counts graphs
library(patchwork)
p1 + p2 + plot_layout(ncol = 1, heights = c(1,0.4), guides = 'collect') 

