library(readxl)
library(plyr)
library(lubridate)

# User modifiable variables
file <- 'genomeMetaData.xlsx'
start_date <- "2021-02-28"

# R Script _____________________________________________________
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
df4 <- rbind(df3a, df3b, df3c)
# Summarize the data
df4$sample_date_week <- floor_date(df4$sample_date, "week") # Adds a column with date for the week only
df5 <- ddply(df4, c("sample_date_week", "lineage"), summarise, n = length(lineage)) # Makes a new dataframe with the counts of each lineage for each of the week categories
df6 <- na.omit(df5) # Remove any data that does not have an associated date
# Make a list of the unique dates and lineages
unique_lineage <- unique(df6[c("lineage")])
unique_date_week <- unique(df6[c('sample_date_week')])
# Make the columns of the graph with the repeating lineages for each date
repeat_lineage <- unique_lineage
for(d in 2:lengths(unique_date_week)){
  repeat_lineage <- Map(c, repeat_lineage, unique_lineage)
}
# Make the column of the graph wit the repeating date, where there are enough dates for each lineage.
repeat_date <- list(sample_date = unique_date_week[1,1])
for(d in 1:lengths(unique_date_week)){
  for(l in 1:lengths(unique_lineage)){
    repeat_date <-Map(c, repeat_date, unique_date_week[d,1])
  }
} 
lengths(repeat_date)

# Remove the initial repeated date from the creation of this new list
repeat_date <- repeat_date$sample_date[-2]
repeat_date <- list(sample_date = repeat_date)
lengths(repeat_date)

# Combine the repeated dates and lineages into a single data frame
df7 <- c(repeat_date, repeat_lineage)
df7 <- data.frame(df7)


# Make a list of zero counts to be temporary placeholders.  
count_zero <- rep(c(0),each=length(df7$sample_date))
count_zero <- list(counts = count_zero)
df8 <- c(repeat_date, repeat_lineage, count_zero)
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

# Assign Data Frame
data <- df9

# Convert Variant Counts to Percentages
data2 <- data  %>%
  group_by(sample_date, lineage) %>%
  summarise(n = sum(counts)) %>%
  mutate(percentage =100 * n / sum(n),
         date=as.Date(sample_date, "%m/%d/%Y")) %>%
  ungroup() %>% as.data.frame()


# Plot Stacked Line Graph
ggplot(data2, aes(x=date, y=percentage, fill=lineage)) + geom_area(alpha=0.6 , size=.5, color="white") + theme_ipsum() +  ggtitle("Philadelphia Random Sampling SARS-CoV-2 Variants")

