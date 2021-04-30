# SarsCoV2-Graphs

These scripts are intended to produce graphs of SARS-CoV-2 lineages derived from John Everett's Pangolin-based pipeline.  

## Getting Started

To begin, direct the genomesMetaData.xlsx file from John Everett's pipeline, to the R working directory. Once this is done, Use R to run the program. 

### Prerequisites

To use the scripts associated with this project, the following R libraries must be installed:

* readxl
* plyr
* dplyr
* lubridate
* scales
* ggplot2
* viridis
* hrbrthemes

## Usage

Run excel2graph.r in R.

## Parameters

* title: The title to be printed above the graph.
* rationale_list: Samples that will be graphed by filtering the rationale type of the sample. 
* graph_type: Whether the lineages will be graphed as a bar graph ("bar") or stacked line graph ("line").
* date_start: The threshold for the first sample date to be graphed.
* date_stop: The threshold for the last sample date to be graphed.
* date_label: The x-axis format for displaying the date. Use the traditional R format for specifying dates.
* time_bin: The data should be binned, either "week" or "month".
* breaks_num: The number of labels that should be made on the X axis.
* bar_width: The bar width of the Genomes bar graph.
* scale: The scale for for the Genomes bar graph, either "linear" or "log".

## Example parameters


* title <- 'Philadelphia SARS-CoV-2 Variants'
* rationale_list <- list('All')
* graph_type <- 'bar'
* date_start <- '2020-01-01'
* date_stop <- '2021-03-27'
* date_label <- "%b %Y"
* time_bin <- 'week'
* breaks_num <- 6
* bar_width <- 5
* scale <- 'log'


## Authors

* **Andrew David Marques** - *Initial work* - [University of Pennsylvania](https://www.linkedin.com/in/andrew-marques-290a29164/)

* **John Everett** - *Collaborator* - [University of Pennsylvania](https://www.linkedin.com/in/everettjohn/)
