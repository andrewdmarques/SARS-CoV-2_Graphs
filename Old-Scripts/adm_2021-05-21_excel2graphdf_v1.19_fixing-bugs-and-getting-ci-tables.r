library(readxl)
library(plyr)
library(dplyr)
library(lubridate)
library(scales)

# ______________________________________________________________________________
# ______________________________________________________________________________
# Define the graphing function
graph_sars <- function(df_original, lineage_important, title, rationale_list, graph_type, date_start, date_stop, date_label, time_bin, breaks_num, bar_width, scale){

# Make a new dataframe with only the necessary columns
if(rationale == 'y'){
  df1 <- df_original[,c("sample_date", "lineage", "rationale")]
}
if(rationale == 'n'){
  df1 <- df_original[,c("sample_date", "lineage")]
}  

# Condense B.1.526.1 and B.1.526.2 into the parent lineage B.1.526
df1[df1=='B.1.526.1']<- 'B.1.526'
df1[df1=='B.1.526.2']<- 'B.1.526'
df1[df1=='B.1.526.3']<- 'B.1.526'
df1[df1=='B.1.617.1']<- 'B.1.617'
df1[df1=='B.1.617.2']<- 'B.1.617'
df1[df1=='B.1.617.3']<- 'B.1.617'

# Remove missing dates from df1
df2 <- na.omit(df1)

# Make the date column a date object
df2$sample_date <- as.Date(df2$sample_date, "%Y-%m-%d")

# Make a data frame for the rationales selected
df3 <- data.frame(Date=as.Date(character()),
                  File=character(), 
                  User=character(), 
                  stringsAsFactors=FALSE)
if(rationale == 'y'){
  for(rl in 1:length(rationale_list)){
    df3 <- rbind(df3, df2[df2$rationale == rationale_list[rl],])
  }
  df4 <- df3
  if(rationale_list[1] == 'All'){
    df4 <- df2
  }
}
if(rationale == 'n'){
  df4 <- df2
}

# Summarize the data
df4$sample_date_week <- floor_date(df4$sample_date, time_bin) # Adds a column with date for the week only
df5 <- ddply(df4, c("sample_date_week", "lineage"), summarise, n = length(lineage)) # Makes a new dataframe with the counts of each lineage for each of the week categories
df6 <- na.omit(df5) # Remove any data that does not have an associated date

# Make a list of the unique dates and lineages
lineage_unique <- unique(df6[c("lineage")])
date_week_unique <- unique(df6[c('sample_date_week')])

# Add important lineages to unique lineages. This makes sure the color assigned to lineages is consistent.
for(li in 1:length(lineage_important)){
  lineage_unique <- lineage_unique %>% add_row(lineage = lineage_important[li])
}
lineage_unique <- unique(lineage_unique[c("lineage")])

# Make the columns of the graph with the repeating lineages for each date
lineage_repeat <- lineage_unique
for(d in 2:lengths(date_week_unique)){
  lineage_repeat <- Map(c, lineage_repeat, lineage_unique)
}

# Make the column of the graph with the repeating date, where there are enough dates for each lineage.
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

# Apply the date cut off to include reads only as old as the date_start specified
df9 <- df8 %>%
  select(sample_date, lineage, counts) %>%
  filter(sample_date >= as.Date(date_start)) %>%
  filter(sample_date <= as.Date(date_stop))

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
  # Lumped "Other"
  mutate(lineage = factor(lineage),
         lineage = forcats::fct_other(f = lineage, drop = lineage_drop, other_level = "Other")) %>%
  group_by(date, lineage) %>%
  summarise(lineage_count = sum(counts, na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(date) %>%
  mutate(lineage_perc = 100 * lineage_count / sum(lineage_count)) %>%
  ungroup() %>%
  identity() -> data1

# Add specimen totals
data1 %>%
  group_by(date) %>%
  summarise(total_seqs = sum(lineage_count, na.rm = TRUE))

# Add colloquial names
data2 <- data1
data2 <- data.frame(lapply(data2, function(x) {gsub("B.1.526", "B.1.526 (NY)", x)}))
data2 <- data.frame(lapply(data2, function(x) {gsub("B.1.1.7", "B.1.1.7 (UK)", x)}))
data2 <- data.frame(lapply(data2, function(x) {gsub("P.1", "P.1 (Brazil)", x)}))
data2 <- data.frame(lapply(data2, function(x) {gsub("B.1.351", "B.1.351 (South Africa)", x)}))
data2 <- data.frame(lapply(data2, function(x) {gsub("B.1.427", "B.1.427 (California)", x)}))
data2 <- data.frame(lapply(data2, function(x) {gsub("B.1.429", "B.1.429 (California)", x)}))
data2 <- data.frame(lapply(data2, function(x) {gsub("B.1.617", "B.1.617 (India)", x)}))
data2$date <- as.Date(data2$date)
data2$lineage <- as.factor(data2$lineage)
data2$lineage_count <- as.numeric(data2$lineage_count)
data2$lineage_perc <- as.numeric(data2$lineage_perc)
data2$lineage <- factor(data2$lineage, levels = c("B.1", "B.1.1.7 (UK)", "B.1.2", "B.1.243", "B.1.311", "B.1.351 (South Africa)", "B.1.427 (California)", "B.1.429 (California)", "B.1.525", "B.1.526 (NY)", "B.1.526.1 (NY)", "B.1.617 (India)", "P.1 (Brazil)", "P.2", "Other"))

# Change the values of the genome counts from 1 to 1.1892071 so that they appear as a value on the genome graph
data2$counts <- data2$lineage_count
data2$lineage_count[data2$lineage_count==1] <- 1.1892071

# Reorder the columns in data2
col_order <- c("date", "lineage", "lineage_count", "counts", "lineage_perc")
data2 <- data2[, col_order]

# ______________________________________________________________________________
# Calculate 95% Confidence Interval using Clopper-Pearson "Exact" confidence interval
# Make a list of zero counts to be temporary placeholders.  
group_zero <- rep(c(0),each=length(data2$counts))
group_count_zero <- list(group_count = group_zero)
ci_lo_zero <- list(lower_ci = group_zero)
ci_up_zero <- list(upper_ci = group_zero)
data3 <- c(data2, group_count_zero, ci_lo_zero, ci_up_zero)
data3 <- data.frame(data3)

# Add the number of counts for each group of dates as a new column
group_new <- TRUE
group_date <- unique(data3$date)
for(g in 1:length(group_date)){
  data3_temp <- data3[data3$date == group_date[g],]
  group_count <- sum(data3_temp$counts)
  data3_temp$group_count[1:length(data3_temp$group_count)] <- group_count
  if(group_new == FALSE){
    data4 <- rbind(data4,data3_temp)
  }
  if(group_new == TRUE){
    group_new <- FALSE
    data4 <- data3_temp
  }
}

# Calculate the upper and lower bounds of the confidence interval.
for(i in 1:length(data4$group_count)){
  ci_all <- binom.test(data4$counts[i], data4$group_count[i], conf.level = 0.95)
  ci_summary <- as.list(ci_all[4])
  ci_upper_lower <- data.frame(ci_summary)
  ci_upper <- ci_upper_lower[2,1]
  ci_lower <- ci_upper_lower[1,1]
  data4$lower_ci[i] <- ci_lower*100
  data4$upper_ci[i] <- ci_upper*100
}

# Remove Unnecessary columns.
data4<- subset(data4, select = -lineage_count)
data4<- subset(data4, select = -counts)
data4<- subset(data4, select = -group_count)

# Rename columns in the dataframe.
data4 <- rename(data4, c("percent_lineage" = "lineage_perc", "lower_95%_CI" = "lower_ci", "upper_95%_CI" = "upper_ci"))

# Save the file as a .csv.
write.csv(data4, gsub(' ', '', paste(gsub(' ', '-', title), '_confidence-interval.csv')), row.names = TRUE)

# ______________________________________________________________________________
# Graphing
# Load libraries required for graphing the data
library(ggplot2)
library(viridis)
library(hrbrthemes)

# Plot lineages using line graph
if(graph_type == 'line'){
  p1 <- ggplot(data2,aes(x=date, y=lineage_perc, fill=lineage)) +
    geom_area(alpha=0.6 , size=.5, color="white") +
    labs(fill = "", x = "Date", y = "Percentage") +
    theme_ipsum() +
    ggtitle("Philadelphia Random Sampling SARS-CoV-2 Variants") +
    scale_x_date(date_labels = date_label, breaks = breaks_pretty(4)) +
    guides(fill=guide_legend(title="lineage")) +
    scale_fill_manual(values = c("yellow2", "palegreen", "lightseagreen", "skyblue", "steelblue1", "royalblue", "navy", "purple", "hotpink", "indianred", "red", "orange","lightgoldenrod1", "honeydew4"))
}

# Plot lineages using bar graph
if(graph_type == 'bar'){
  p1 <- ggplot(data2, aes(fill=lineage, y=lineage_perc, x=date)) +
    geom_bar(position="stack", stat="identity") +
    theme_ipsum() +
    ggtitle(title) +
    labs(x = "Date", y = "Percentage") +
    scale_x_date(date_labels = date_label, breaks = breaks_pretty(breaks_num)) + 
    scale_fill_manual(values = c("yellow2", "palegreen", "lightseagreen", "skyblue", "steelblue1", "royalblue", "navy", "purple", "hotpink", "indianred", "red", "orange","lightgoldenrod1", "honeydew4"))
}

# Plot sequencing counts using bar graph
p2 <- data2 %>%
  group_by(date) %>%
  summarise(total_seqs = sum(lineage_count, na.rm = TRUE)) %>%
  ggplot(data = ., aes(x = date, y = total_seqs)) +
  geom_bar(position="stack", stat="identity", width = bar_width) +
  theme_ipsum() +
  labs(x = "Date", y = "Genomes") +
  scale_x_date(date_labels = date_label, breaks = breaks_pretty(breaks_num)) +
  if(scale == 'log'){
    scale_y_continuous(trans = 'log2')
  }

# Combine lineages and sequence counts graphs
library(patchwork)
p1 + p2 + plot_layout(ncol = 1, heights = c(1,0.4), guides = 'collect') 
}

# ______________________________________________________________________________
# ______________________________________________________________________________
# User modifiable variables

# General parameters
lineage_important <- c("B.1.526", "B.1.525", "P.2", "B.1.1.7", "P.1", "B.1.351", "B.1.427", "B.1.429", "B.1.311",  "B.1", "B.1.2", "B.1.243", "B.1.617")

# Graph 1: Philadelphia All Samples
title <- 'Delaware Valley SARS-CoV-2 Variants'
file <- '2021-05-21_genomeMetaData.xlsx'
df_original <- read_excel(file)
rationale <- 'y'
rationale_list <- list('All')
graph_type <- 'bar'
date_start <- '2020-01-01'
date_stop <- '2021-05-22'
date_label <- "%b %Y"
time_bin <- 'week'
breaks_num <- 6
bar_width <- 5
scale <- 'log'
g1 <- graph_sars(df_original, lineage_important, title, rationale_list, graph_type, date_start, date_stop, date_label, time_bin, breaks_num, bar_width, scale)
g1

# Graph 2: Washington DC
title <- 'Washington DC SARS-CoV-2 Variants'
file <- 'wdc_2021-05-13_all_lineage_report_dated.csv'
df_original <- read.csv(file)
# date_in_file <- "%m/%"
rationale <- 'n'
rationale_list <- list('All')
graph_type <- 'bar'
date_start <- '2020-01-01'
date_stop <- '2021-06-27'
date_label <- "%b %Y"
time_bin <- 'week'
breaks_num <- 6
bar_width <- 5
scale <- 'log'
g2 <- graph_sars(df_original, lineage_important, title, rationale_list, graph_type, date_start, date_stop, date_label, time_bin, breaks_num, bar_width, scale)
g2

# Graph 3: NYC
title <- 'NYC SARS-CoV-2 Variants'
file <- 'nyc_repeat_all_lineage_report_dated.csv'
df_original <- read.csv(file)
# date_in_file <- "%m/%"
rationale <- 'n'
rationale_list <- list('All')
graph_type <- 'bar'
date_start <- '2020-01-01'
date_stop <- '2021-06-27'
date_label <- "%b %Y"
time_bin <- 'week'
breaks_num <- 6
bar_width <- 5
scale <- 'log'
g4_repeat <- graph_sars(df_original, lineage_important, title, rationale_list, graph_type, date_start, date_stop, date_label, time_bin, breaks_num, bar_width, scale)
g4_repeat
