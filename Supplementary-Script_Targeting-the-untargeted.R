#General information####
#R-Script for the publication "Targeting the untargeted:  Uncovering the chemical complexity of root exudates"
#Katrin Möller1*, Annalena Ritter2*, Phillip J. Stobinsky2*, Kai Jensen1, Ina C. Meier2, Dirk Granse1, Harihar Jaishree Subrahmaniam1,2 

#1 Applied Plant Ecology, Institut für Pflanzenwissenschaften und Mikrobiologie, University of Hamburg, Ohnhorststraße 18, 22609 Hamburg 
#2 Functional Forest Ecology, Institut für Pflanzenwissenschaften und Mikrobiologie University of Hamburg, Ohnhorststraße 18, 22609 Hamburg 

#* K. Möller, A. Ritter, and P. J. Stobinsky contributed equally to this publication 

#Author for correspondence 
#Harihar Jaishree Subrahmaniam, Tel.: +49 40 42816-722, 
#Email: h.subrahmaniam@uni-hamburg.de, jaishree.subrahmaniam@gmail.com

#ORCID IDs 
#Katrin Möller: orcid.org/0009-0005-4457-9549 
#Annalena Ritter: orcid.org/0000-0003-1314-2705  
#Phillip J. Stobinsky: orcid.org/0009-0008-9804-9908 
#Kai Jensen: orcid.org /0000-0002-0543-070X  
#Ina C. Meier: orcid.org/0000-0001-6500-7519 
#Dirk Granse: orcid.org/0000-0002-8016-0703 
#Harihar Jaishree Subrahmaniam: orcid.org/0000-0002-1077-5024 

#Setup####
```{r Setup, Loading in all needed packages and the data file}
#Checking the working direktory
getwd()

#Loading all necessary packages
library(readxl)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(DataExplorer)
library(ggpubr)
library(multcomp)
library(multcompView)
library(broom)
library(tinytex)
library(bookdown)
library(GGally)
library(gridExtra)
library(ggforce)
library(scales)
library(patchwork)
library(writexl)
library(viridis)
library(RColorBrewer)
library(ggsci)
library(ggtree)
library(grid)
library(gtable)
library(egg)
library(ape)
library(pheatmap)
library(rotl)
library(cowplot)
library(ggdendro)
library(dendextend)
library(stringr)
library(phytools)

#Clearing the global environment
#rm(list=ls())

#Loading in the main data file
data <- read_excel("Supplementary-Table-1_Targeting-the-untargeted (1).xlsx")
```

#Figure creation####
```{r Figure 1, Plotting the publications in plant ecology (bar graph) and publications on untargeted root exudate metabolomics (line graph) over time, combined in one figure}
#Loading in a text file with data on the publications in plant ecology over time 
publications <- read_delim("C:/Root exudate chemistry commentary paper/Publications in plant ecology over time.txt", 
                           delim = "\t", escape_double = FALSE, 
                           trim_ws = TRUE)

#Calculating the number of publications on untargeted root exudate metabolomics per year
publication_data <- data %>%
  group_by(DOI, `Publication year`)%>% 
  summarise(dummy = 1) %>% 
  group_by(`Publication year`) %>%
  summarise(exudates_pub = sum(dummy)) 

#Creating a tibble with the needed variables for the figure, namely the cumulative sum of publications in plant ecology and untargeted root exudate metabolomics research
pub_data_full <- publications %>%
  full_join(publication_data, by = join_by(`Final Publication Year` == `Publication year`)) %>%
  arrange(`Final Publication Year`) %>%
  mutate_all(~replace(., is.na(.), 0)) %>% #replacing all empty cells with 0, else publication years
  mutate(exudate_cumsum = cumsum(exudates_pub),
         ecology_cumsum = cumsum(`Record Count`))

#Creating figure 1
plot_pub_per_year <- pub_data_full %>%
  ggplot(aes(x = `Final Publication Year`)) +
  geom_line(aes(y = ecology_cumsum / 2000), color = "red") +
  geom_col(aes(y = exudate_cumsum)) +
  theme_minimal() +
  labs(x = "Year") +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),   
    axis.line = element_line(), 
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), 
    axis.title.x = element_text(),
    axis.ticks = element_line(color = "black"), 
    axis.ticks.length = unit(0.15, "cm")) +
  scale_x_continuous(breaks = seq(2003, 2025, by = 1), expand = c(0,0)) +
  scale_y_continuous(name = "Cumulative sum of publications on\nuntargeted root exudate metabolomics (bars)", 
                     limits = c(0,60),
                     sec.axis = sec_axis(~ . * 2000, name = "Cumulative sum of publications on plant ecology (line)"), 
                     expand = c(0,0))

#Plotting figure 1
plot_pub_per_year
```

```{r Figure 2, Plotting plant utilisation and orders as stacked bar graphs}
#Counting distinct plant utilisation - plant order combinations
plant_type_order_count <- data %>%
  distinct(`DOI`, `Plant species`, `Plant type`, `Plant order`) %>%
  group_by(`Plant type`, `Plant order`) %>%
  summarise(Count = n()) %>%
  ungroup()

#Identifying and labelling plant orders that appear only once within each plant utilisation group
plant_type_order_count <- plant_type_order_count %>%
  group_by(`Plant type`) %>%
  mutate(Plant_order = ifelse(Count == 1, "Other", `Plant order`)) %>%
  ungroup()

#Reordering the plant utilisation bar graphs to show the group "Other" on the right most position on the x-axis
plant_type_order_count_reorder <- plant_type_order_count %>%
  group_by(`Plant type`, Plant_order) %>%
  summarise(Total_Count = sum(Count)) %>%
  ungroup() %>%
  mutate(`Plant type` = factor(`Plant type`, levels = c("Crop", "Model plant", "Other")))

#Reorder the plant orders by size in the largest bar graph (plant utilisation = Crop) and the other bar graphs according to this order
plant_type_order_count_reorder <- plant_type_order_count_reorder %>%
  arrange(`Plant type`, -Total_Count) %>%
  mutate(Plant_order = factor(Plant_order, levels = rev(unique(Plant_order))))

#Defining a colour-blind-friendly colour-palette used in the figures 2-4
custom_palette = brewer.pal(11, "Set3")

#Creating figure 2
plot_plant_type_order <- ggplot(plant_type_order_count_reorder, aes(x = `Plant type`, y = Total_Count, fill = Plant_order)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_minimal() +
  labs(x = "Plant type",
       y = "Number of experiments",
       fill = "Plant order") +
  theme(plot.title = element_blank(), 
        panel.grid.major = element_blank(),    
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(),
        axis.line = element_line(),    
        axis.ticks = element_line(color = "black")) +
  scale_fill_manual(values = custom_palette, 
                    labels = c("Other" = "Other order")) +
  scale_x_discrete(labels = c("Crop", "Model plant", "Other type"), expand = c(0.01, 0)) + 
  scale_y_continuous(breaks = seq(0, 70, by = 5), limits = c(0, 70), expand = c(0, 0)) 


#Plotting figure 2
plot_plant_type_order
```

```{r Figure 3, Plotting root exudate sampling approaches and solutions as stacked bar graphs}
#Choosing all relevant rows from the data file with distinct DOI, sampling approaches and sampling solutions
method_per_study <- data %>%
  distinct(`DOI`, `Sampling medium`, `Sampling approach`)

#Counting the number of studies using a distinct sampling media - sampling solution combination
n_per_method <- method_per_study %>%
  group_by(`Sampling medium`, `Sampling approach`) %>%
  summarize(`DOI` = n(), .groups = 'drop')

#Calculating the total number of studies using a distinct sampling method
total_per_method <- n_per_method %>%
  group_by(`Sampling approach`) %>%
  summarize(total_DOI = sum(`DOI`), .groups = 'drop')

#Reordering the sampling method based on the total number of studies using a distinct method
n_per_method <- n_per_method %>%
  mutate(`Sampling approach` = fct_reorder(`Sampling approach`, total_per_method$total_DOI[match(`Sampling approach`, total_per_method$`Sampling approach`)], .desc = TRUE))

#Creating figure 3
plot_n_per_method <- ggplot(n_per_method, aes(x = `Sampling approach`, y = `DOI`, fill = `Sampling medium`)) +
  geom_bar(stat = "identity") +
  labs(title = "Root exudate sampling approach and media",
       x = "Root exudate sampling approach",
       y = "Number of experiments",
       fill = "Root exudate sampling medium") +
  theme_minimal() +
  theme(plot.title = element_blank(),  
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.title.x = element_text(), 
        axis.title.y = element_text(vjust = 2),   
        axis.line = element_line(),
        axis.ticks = element_line(color = "black")) +
  scale_fill_manual(values = custom_palette) +  # Apply the custom 10-color palette
  scale_x_discrete(labels = c("Hydroponic", "Soil-hydroponic-hybrid", "Soil-based", "Agar", "Sorption", "Aeroponic"),  expand = c(0,0)) +
  scale_y_continuous(breaks = seq(0, 30, by = 5), limits = c(0, 30), expand = c(0, 0))

#Plotting figure 3
plot_n_per_method
```

```{r Figure 4, Plotting the root exudate sampling durations and sampling approach as stacked bar graphs}
#Choosing all relevant rows from the data file with distinct DOI, sampling approach and sampling duration
collectiontime_per_approach <- data %>%
  distinct(`DOI`, `Sampling duration`, `Sampling approach`)

#Defining an order for collection time displayed on the x-axis
collectiontime_per_approach <- collectiontime_per_approach %>%
  mutate(`Sampling duration` = factor(`Sampling duration`, 
                                      levels = c("<1 min", "2 min", "15 min", "1 h", "2 h", "3 h", "4 h", "6 h", "8 h", "10 h", "12 h", "22 h", "1 d", "1.25 d", "1.71 d", "2 d", "3 d", "6 d", "7 d", "13 d", "14 d", "21 d", "28 d", "32 d", "42 d", "45 d", "56 d", "undefinable", "N.A.")))

#Reordering the sampling approaches in the legend by size from top to bottom
collectiontime_per_approach <- collectiontime_per_approach %>%
  mutate(`Sampling approach` = fct_reorder(`Sampling approach`, `Sampling approach`, .fun = length, .desc = TRUE)) 

#Creating figure 4
plot_collectiontime_per_approach <- ggplot(collectiontime_per_approach, aes(x = `Sampling duration`, fill = `Sampling approach`)) +
  geom_bar(stat = "count") +
  theme_minimal() +
  labs(title = "Root exudate collection time per sampling method",
       x = "Root exudate sampling approach",
       y = "Number of experiments",
       fill = "Sampling approach") +
  theme(plot.title = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1),  
        axis.title.x = element_text(),
        axis.title.y = element_text(vjust = 2),   
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(),   
        axis.ticks = element_line(color = "black")) + 
  scale_fill_manual(values = custom_palette) +
  scale_x_discrete(labels = c("<1 min", "2 min", "15 min", "1 h", "2 h", "3 h", "4 h", "6 h", "8 h", "10 h", "12 h", "22 h", "1 d", "1.25 d", "1.71 d", "2 d", "3 d", "6 d", "7 d", "13 d", "14 d", "21 d", "28 d", "32 d", "42 d", "45 d", "56 d", "undefinable", "N.R."), expand = c(0,0)) +
  scale_y_continuous(breaks = seq(0, 12, by = 2), limits = c(0, 12), expand = c(0,0))

#Plotting figure 4
plot_collectiontime_per_approach
```

```{r Figure 5, Plotting the distribution of annotated metabolites to the total of reported metabolites in percent by analytical method as stacked bar graphs}
#Creating a custom colour-blind friendly colour-palette used in the figures 5+6
npg_colors <- pal_npg("nrc")(4)

#Assigning specific colors to each analytical method 
custom_colors <- c(
  "GC-MS" = npg_colors[1],
  "LC-MS" = npg_colors[2],
  "Mixed" = npg_colors[3],
  "Other" = npg_colors[4]
)

#Choosing all relevant rows from the data file and filtering out rows where the percentage of annotated compounds is not given (= NA)
percent_annotated_noNA <- data %>%
  distinct(`DOI`, `Sum of reported compounds`, `Sum of annotated compounds`, `Annotated of reported compounds [%]`, `Analytical method`, `Carbohydrates [n]`, `Fatty acids [n]`, `Amino acids and Peptides [n]`, `Terpenoids [n]`, `Shikimates and Phenylpropanoids [n]`, `Alkaloids [n]`, `Polyketides [n]`, `Other [n]`) %>%
  filter(!is.na(`Annotated of reported compounds [%]`)) 

#Creating figure 5
plot_percent_distribution <- ggplot(percent_annotated_noNA, aes(x = `Annotated of reported compounds [%]`, fill = `Analytical method`)) +
  geom_histogram(binwidth = 5, boundary = 0, position = "stack") +  
  theme_minimal() +
  labs(title = "Annotated metabolites [%]",
       x = "Annotated metabolites [%]",
       y = "Number of analyses",
       fill = "Analytical method") +
  theme(plot.title = element_blank(), 
        axis.title.x = element_text(size = 17),
        axis.title.y = element_text(size = 17),
        axis.text.x = element_text(size = 12, angle = 45, vjust = 1, hjust = 1),
        axis.text.y = element_text(size = 12),
        legend.text = element_text(size = 17), 
        legend.title = element_text(size = 17),
        panel.grid.major = element_blank(),    
        panel.grid.minor = element_blank(),   
        axis.line = element_line(),    
        axis.ticks = element_line(color = "black")) + 
  scale_fill_manual(values = custom_colors) +
  scale_x_continuous(breaks = seq(0, 100, by = 5), limits = c(0,100), expand = c(0,0)) + 
  scale_y_continuous(breaks = seq(0,22, by = 5), limits = c(0,22), expand = c(0,0)) +
  annotate("text", x = mean(range(percent_annotated_noNA$`Annotated of reported compounds [%]`)), y = Inf, vjust = 1, label = "(A)", size = 8)

#Plotting figure 5
plot_percent_distribution
```

```{r Figure 5.1, Plotting the distribution of annotated metabolites to the total of reported metabolites in percent by analytical method as stacked bar graphs}
#Choosing all relevant rows from the data file, filtering out rows where the percentage of annotated compounds is not given (= NA) and selecting all rows where the number of reported compounds is under or equally to 100 
percent_annotated_noNA_u100 <- data %>%
  distinct(`DOI`, `Sum of reported compounds`, `Sum of annotated compounds`, `Annotated of reported compounds [%]`, `Analytical method`, `Carbohydrates [n]`, `Fatty acids [n]`, `Amino acids and Peptides [n]`, `Terpenoids [n]`, `Shikimates and Phenylpropanoids [n]`, `Alkaloids [n]`, `Polyketides [n]`, `Other [n]`) %>%
  filter(!is.na(`Annotated of reported compounds [%]`)) %>%
  filter(`Sum of reported compounds` <= 100)

#Creating figure 5.1
plot_percent_distribution_u100 <- ggplot(percent_annotated_noNA_u100, aes(x = `Annotated of reported compounds [%]`, fill = `Analytical method`)) +
  geom_histogram(binwidth = 5, boundary = 0, position = "stack") +  # Stacked histogram
  theme_minimal() +
  labs(title = "Annotated metaboites [%] where -Total reported compounds- ≤100",
       x = "Annotated metabolites [%]",
       y = "Number of analyses \n(Total reported metabolites ≤100)",
       fill = "Analytical method") +
  theme(plot.title = element_blank(), 
        axis.title.x = element_text(size = 17),
        axis.title.y = element_text(size = 17),
        axis.text.x = element_text(size = 12, angle = 45, vjust = 1, hjust = 1),
        axis.text.y = element_text(size = 12),
        legend.text = element_text(size = 17), 
        legend.title = element_text(size = 17),
        panel.grid.major = element_blank(),    
        panel.grid.minor = element_blank(),   
        axis.line = element_line(),    
        axis.ticks = element_line(color = "black")) + 
  scale_fill_manual(values = custom_colors) +
  scale_x_continuous(breaks = seq(0, 100, by = 5), limits = c(0,100), expand = c(0,0)) + 
  scale_y_continuous(breaks = seq(0,22, by = 5), limits = c(0,22), expand = c(0,0)) +
  annotate("text", x = mean(range(percent_annotated_noNA$`Annotated of reported compounds [%]`)), y = Inf, vjust = 1, label = "(C)", size = 8)

#Plotting figure 5.1
plot_percent_distribution_u100
```

```{r Figure 5.2, Plotting the distribution of annotated metabolites to the total of reported metabolites in percent by analytical method as stacked bar graphs1}
#Choosing all relevant rows from the data file, filtering out rows where the percentage of annotated compounds is not given (= NA) and selecting all rows where the number of reported compounds is above to 100 
percent_annotated_noNA_a100 <- data %>%
  distinct(`DOI`, `Sum of reported compounds`, `Sum of annotated compounds`, `Annotated of reported compounds [%]`, `Analytical method`, `Carbohydrates [n]`, `Fatty acids [n]`, `Amino acids and Peptides [n]`, `Terpenoids [n]`, `Shikimates and Phenylpropanoids [n]`, `Alkaloids [n]`, `Polyketides [n]`, `Other [n]`) %>%
  filter(!is.na(`Annotated of reported compounds [%]`)) %>%
  filter(`Sum of reported compounds` > 100)

#Creating figure 5.2
plot_percent_distribution_a100 <- ggplot(percent_annotated_noNA_a100, aes(x = `Annotated of reported compounds [%]`, fill = `Analytical method`)) +
  geom_histogram(binwidth = 5, boundary = 0, position = "stack") +  # Stacked histogram
  theme_minimal() +
  labs(title = "Annotated metabolites [%] where -Total reported compounds- >100",
       x = "Annotated metabolites [%]",
       y = "Number of analyses \n(Total reported metabolites >100)",
       fill = "Analytical method") +  # Label for the fill legend
  theme(plot.title = element_blank(), 
        axis.title.x = element_text(size = 17),
        axis.title.y = element_text(size = 17),
        axis.text.x = element_text(size = 12, angle = 45, vjust = 1, hjust = 1),
        axis.text.y = element_text(size = 12),
        legend.text = element_text(size = 17), 
        legend.title = element_text(size = 17),
        panel.grid.major = element_blank(),    
        panel.grid.minor = element_blank(),   
        axis.line = element_line(),    
        axis.ticks = element_line(color = "black")) + 
  scale_fill_manual(values = custom_colors) +  # Apply the custom color mapping
  scale_x_continuous(breaks = seq(0, 100, by = 5), limits = c(0,100), expand = c(0,0)) + 
  scale_y_continuous(breaks = seq(0,22, by = 5), limits = c(0,22), expand = c(0,0)) +
  annotate("text", x = mean(range(percent_annotated_noNA$`Annotated of reported compounds [%]`)), y = Inf, vjust = 1, label = "(B)", size = 8)

#Plotting figure 5.2
plot_percent_distribution_a100
```

```{r Figure 6, Plotting the mean contribution of compound classes with standard deviation to the total number of annotated compounds; additionally calculating the median of these}
#Filtering out the analytical methods "mixed" and "other", choosing all relevant rows from the data file where the total contribution percentages are under 100 % and calculating the mean, median and standard deviation of all relevant metabolite classes
summary_table_percent_method <- data %>%
  filter(`Analytical method` != "Mixed" & `Analytical method` != "Other") %>%
  distinct(`DOI`,
           `Sum of reported compounds`, 
           `Sum of annotated compounds`, 
           `Carbohydrates of annotated [%]`, 
           `Fatty acids of annotated [%]`, 
           `Amino acids and Peptides of annotated [%]`, 
           `Terpenoids of annotated [%]`, 
           `Shikimates and Phenylpropanoids of annotated [%]`, 
           `Alkaloids of annotated [%]`, 
           `Polyketides of annotated [%]`, 
           `Other of annotated [%]`, 
           `Analytical method`) %>%
  filter(!if_any(
    c(`Carbohydrates of annotated [%]`, 
      `Fatty acids of annotated [%]`, 
      `Amino acids and Peptides of annotated [%]`, 
      `Terpenoids of annotated [%]`, 
      `Shikimates and Phenylpropanoids of annotated [%]`, 
      `Alkaloids of annotated [%]`, 
      `Polyketides of annotated [%]`, 
      `Other of annotated [%]`), 
    ~ . > 100)) %>%
  group_by(`Analytical method`) %>%
  summarise(
    mean_carbohydrate = mean(`Carbohydrates of annotated [%]`, na.rm = TRUE),
    median_carbohydrate = median(`Carbohydrates of annotated [%]`, na.rm = TRUE),
    sd_carbohydrate = sd(`Carbohydrates of annotated [%]`, na.rm = TRUE),
    mean_fatty_acid = mean(`Fatty acids of annotated [%]`, na.rm = TRUE),
    median_fatty_acid = median(`Fatty acids of annotated [%]`, na.rm = TRUE),
    sd_fatty_acid = sd(`Fatty acids of annotated [%]`, na.rm = TRUE),
    mean_amino_acids = mean(`Amino acids and Peptides of annotated [%]`, na.rm = TRUE),
    median_amino_acids = median(`Amino acids and Peptides of annotated [%]`, na.rm = TRUE),
    sd_amino_acids = sd(`Amino acids and Peptides of annotated [%]`, na.rm = TRUE),
    mean_terpenoid = mean(`Terpenoids of annotated [%]`, na.rm = TRUE),
    median_terpenoid = median(`Terpenoids of annotated [%]`, na.rm = TRUE),
    sd_terpenoid = sd(`Terpenoids of annotated [%]`, na.rm = TRUE),
    mean_shikimates = mean(`Shikimates and Phenylpropanoids of annotated [%]`, na.rm = TRUE),
    median_shikimates = median(`Shikimates and Phenylpropanoids of annotated [%]`, na.rm = TRUE),
    sd_shikimates = sd(`Shikimates and Phenylpropanoids of annotated [%]`, na.rm = TRUE),
    mean_alkaloid = mean(`Alkaloids of annotated [%]`, na.rm = TRUE),
    median_alkaloid = median(`Alkaloids of annotated [%]`, na.rm = TRUE),
    sd_alkaloid = sd(`Alkaloids of annotated [%]`, na.rm = TRUE),
    mean_polyketide = mean(`Polyketides of annotated [%]`, na.rm = TRUE),
    median_polyketide = median(`Polyketides of annotated [%]`, na.rm = TRUE),
    sd_polyketide = sd(`Polyketides of annotated [%]`, na.rm = TRUE),
    mean_other = mean(`Other of annotated [%]`, na.rm = TRUE),
    median_other = median(`Other of annotated [%]`, na.rm = TRUE),
    sd_other = sd(`Other of annotated [%]`, na.rm = TRUE)
  )


#Choosing all relevant rows from the data file where the total contribution percentages are under 100 % and calculating the mean, median and standard deviation of all relevant metabolite classes
summary_table_all_methods <- data %>%
  distinct(`DOI`, 
           `Sum of reported compounds`, 
           `Sum of annotated compounds`, 
           `Carbohydrates of annotated [%]`, 
           `Fatty acids of annotated [%]`, 
           `Amino acids and Peptides of annotated [%]`, 
           `Terpenoids of annotated [%]`, 
           `Shikimates and Phenylpropanoids of annotated [%]`, 
           `Alkaloids of annotated [%]`, 
           `Polyketides of annotated [%]`, 
           `Other of annotated [%]`) %>%
  subset(!is.na(`Carbohydrates of annotated [%]`)) %>%
  filter(!if_any(
    c(`Carbohydrates of annotated [%]`, 
      `Fatty acids of annotated [%]`, 
      `Amino acids and Peptides of annotated [%]`, 
      `Terpenoids of annotated [%]`, 
      `Shikimates and Phenylpropanoids of annotated [%]`, 
      `Alkaloids of annotated [%]`, 
      `Polyketides of annotated [%]`, 
      `Other of annotated [%]`), 
    ~ . > 100)) %>%
  summarise(
    mean_carbohydrate = mean(`Carbohydrates of annotated [%]`, na.rm = TRUE),
    median_carbohydrate = median(`Carbohydrates of annotated [%]`, na.rm = TRUE),
    sd_carbohydrate = sd(`Carbohydrates of annotated [%]`, na.rm = TRUE),
    mean_fatty_acid = mean(`Fatty acids of annotated [%]`, na.rm = TRUE),
    median_fatty_acid = median(`Fatty acids of annotated [%]`, na.rm = TRUE),
    sd_fatty_acid = sd(`Fatty acids of annotated [%]`, na.rm = TRUE),
    mean_amino_acids = mean(`Amino acids and Peptides of annotated [%]`, na.rm = TRUE),
    median_amino_acids = median(`Amino acids and Peptides of annotated [%]`, na.rm = TRUE),
    sd_amino_acids = sd(`Amino acids and Peptides of annotated [%]`, na.rm = TRUE),
    mean_terpenoid = mean(`Terpenoids of annotated [%]`, na.rm = TRUE),
    median_terpenoid = median(`Terpenoids of annotated [%]`, na.rm = TRUE),
    sd_terpenoid = sd(`Terpenoids of annotated [%]`, na.rm = TRUE),
    mean_shikimates = mean(`Shikimates and Phenylpropanoids of annotated [%]`, na.rm = TRUE),
    median_shikimates = median(`Shikimates and Phenylpropanoids of annotated [%]`, na.rm = TRUE),
    sd_shikimates = sd(`Shikimates and Phenylpropanoids of annotated [%]`, na.rm = TRUE),
    mean_alkaloid = mean(`Alkaloids of annotated [%]`, na.rm = TRUE),
    median_alkaloid = median(`Alkaloids of annotated [%]`, na.rm = TRUE),
    sd_alkaloid = sd(`Alkaloids of annotated [%]`, na.rm = TRUE),
    mean_polyketide = mean(`Polyketides of annotated [%]`, na.rm = TRUE),
    median_polyketide = median(`Polyketides of annotated [%]`, na.rm = TRUE),
    sd_polyketide = sd(`Polyketides of annotated [%]`, na.rm = TRUE),
    mean_other = mean(`Other of annotated [%]`, na.rm = TRUE),
    median_other = median(`Other of annotated [%]`, na.rm = TRUE),
    sd_other = sd(`Other of annotated [%]`, na.rm = TRUE)
  ) %>%
  pivot_longer(cols = everything(), names_to = c(".value", "Component"), names_pattern = "(mean|sd|median)_(.*)") %>%
  mutate(`Analytical method` = "All methods combined")


#Ordering the compound classes by the mean from calculations including all analytical methods
ordered_components <- summary_table_all_methods %>%
  arrange(desc(mean)) %>%
  pull(Component)

#Changing the table with excluded analytical methods into long format
long_format_method <- summary_table_percent_method %>%
  pivot_longer(
    cols = starts_with("mean_") | starts_with("sd_"),
    names_to = c(".value", "Component"),
    names_pattern = "(mean|sd)_(.*)"
  )

#Combining both tables with excluded analytical methods and all methods included
combined_summary <- long_format_method %>%
  bind_rows(summary_table_all_methods)

#Ensuring the component is ordered by the mean value from the table with all methods included
combined_summary <- combined_summary %>%
  mutate(Component = factor(Component, levels = ordered_components))


#Assigning specific colors to each analytical method
custom_colors_1 <- c(
  "GC-MS" = npg_colors[1],
  "LC-MS" = npg_colors[2],
  "All methods combined" = npg_colors[5]
)

#Creating figure 6
plot_contr_com_method <- ggplot(combined_summary, aes(x = Component, y = mean, fill = `Analytical method`)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = mean, ymax = mean + sd), 
                position = position_dodge(width = 0.9), width = 0) +
  geom_errorbar(aes(ymin = mean + sd, ymax = mean + sd), 
                position = position_dodge(width = 0.9), width = 0.25) +
  labs(title = "Mean contribution of compound classes grouped by method with mean and SD",
       x = "Compound class",
       y = "Mean contribution to annotated metabolites [%]",
       fill = "Analytical method") +
  theme_minimal() +
  theme(plot.title = element_blank(), 
        axis.title.x = element_blank(),
        axis.title.y = element_text(vjust = 2),   
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "right",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), 
        axis.line = element_line(),  
        axis.ticks = element_line(color = "black")) +
  scale_fill_manual(values = custom_colors_1) +
  scale_x_discrete(labels = c("Other", "Shikimates and  \n phenylpropanoids", "Carbohydrates", "Amino acids and Peptides", "Fatty acids", "Polyketides", "Terpenoids", "Alkaloids"), expand = c(0,0)) +
  scale_y_continuous(breaks = seq(0,80, by = 10), limits = c(0,80), expand = c(0,0))

#Plotting figure 6
plot_contr_com_method
```

```{r Figure 7, Plotting the distribution of annotated metabolites to the total of reported metabolites in percent by analytical method as stacked bar graphs}
initDB <- data.frame(Type = 'LC-MS', Label = '(A) LC-MS', 
                     includeOthers = F,
                     excludeFirstAuthor = c('Herz'),
                     excludeSpecies = c('Dactylis glomerata'),
                     FigureHeight = 185, FigureWidth = 160,
                     C_TextDispWidthFactor = 4.4, #3.6,
                     C_TextBackgroundWidthFactor = 3.5, #3.3,
                     C_DO_DISPLAY_LEGEND_IN_ALL_PLOTS = F,
                     C_TopLabelPanelHeight = 0.06)
initDB <- rbind(initDB, data.frame(Type = 'GC-MS', Label = '(B) GC-MS', 
                                   includeOthers = F,
                                   excludeFirstAuthor = c('Herz'),
                                   excludeSpecies = c(''),
                                   FigureHeight = 150, FigureWidth = 160,
                                   C_TextDispWidthFactor = 4,
                                   C_TextBackgroundWidthFactor = 3.5, 
                                   C_DO_DISPLAY_LEGEND_IN_ALL_PLOTS = T,
                                   C_TopLabelPanelHeight = 0.07))

# Define your species list
species_list <- c(
  'Pennisetum glaucum', 
  'Brassica napus', 
  'Arabidopsis thaliana', 
  'Medicago sativa', 
  'Pinus sylvestris', 
  'Zea mays', 
  'Cunninghamia lanceolata ', 
  # 'Sorghum spp.', 
  'Achillea millefolium', 
  'Alopecurus pratensis', 
  'Arrhenatherum elatius', 
  'Dactylis glomerata', 
  'Galium mollugo', 
  'Galium verum', 
  'Lolium perenne', 
  'Plantago lanceolata', 
  'Poa pratensis', 
  'Ranunculus acris', 
  'Solanum lycopersicum', 
  'Cyperus alternifolius', 
  'Pistia stratiotes', 
  'Cucumis sativus', 
  'Pisum sativum', 
  'Vicia faba', 
  'Lupinus albus', 
  'Lolium arundinaceum', 
  'Avena strigosa', 
  'Phacelia tanacetifolia', 
  'Sinapis alba', 
  'Trifolium alexandrinum', 
  'Chrysopogon zizanioides', 
  'Solanum tuberosum', 
  'Musa acuminata', 
  'Robinia pseudoacacia', 
  'Triticum aestivum', 
  'Gossypium hirsutum', 
  #'Vitis riparia × Vitis labrusca', 
  'Sedum alfredii', 
  #'Vitis spp.', 
  'Hordeum vulgare', 
  'Brachypodium distachyon', 
  'Medicago truncatula', 
  'Sorghum bicolor', 
  'Glycine max', 
  'Nicotiana tabacum', 
  'Panax notoginseng', 
  'Perilla frutescens', 
  'Bouteloua gracilis', 
  'Panicum virgatum', 
  'Vigna radiata', 
  'Vigna unguiculata', 
  'Cedrus deodara', 
  'Cinnamomum camphora', 
  #  'Liriodendron chinense × tulipifera', 
  'Metasequoia glyptostroboides', 
  'Populus euramericana', 
  'Pterocarya stenoptera', 
  #'Sabina chinensis', 
  'Avena barbata', 
  'Solanum melongena',
  'Juniperus chinensis'
)

if (1 == 2) {
  # Match species names with Open Tree of Life taxonomy
  resolved_names <- tnrs_match_names(species_list)  # Warning message: Some names were duplicated: ‘sorghum bicolor’
  
  # Filter out problematic OTT IDs
  #problematic_ott_ids <- c("3915043", "329913")
  #resolved_names[resolved_names$ott_id %in% problematic_ott_ids, ]
  #valid_resolved_names <- resolved_names[!resolved_names$ott_id %in% problematic_ott_ids, ]
  valid_resolved_names <- valid_resolved_names[!is.na(valid_resolved_names$ott_id), ]
  
  # Get the OTT IDs for the matched species
  ott_ids <- valid_resolved_names$ott_id
  
  # Retrieve the phylogenetic tree
  if (length(ott_ids) > 0) {
    tree <- tol_induced_subtree(ott_ids = ott_ids)
    plot(tree, main = "Phylogenetic Tree from Open Tree of Life", cex = 0.7)
    
    # Save the tree in Newick format
    write.tree(tree, file = "phylogenetic_tree_from_OTL.newick")
  } else {
    message("No valid OTT IDs to generate a tree.")
  }
}

for (mT in 1:nrow(initDB)) {
  
  #mT <- 1
  methodType <- initDB[mT, 'Type']
  
  # Define the compound contribution columns
  
  compound_columns <- c(
    "Carbohydrates",
    "Fatty acids",
    "Amino acids and Peptides",
    "Terpenoids",
    "Shikimates and Phenylpropanoids",
    "Alkaloids",
    "Polyketides"
  )
  
  if (initDB[mT, 'includeOthers']) { displayColnames <- c(displayColnames, "Others") }
  compound_columns <- paste0(compound_columns, ' of annotated [%]')
  
  # Shorten the compound class names for better display in the heatmap
  displayColnames <- c(
    "Carbohydrates",
    "Fatty acids",
    "Amino acids and peptides",
    "Terpenoids",
    "Shikimates & Phenylprop.",
    "Alkaloids",
    "Polyketides"
  )
  if (initDB[mT, 'includeOthers']) { displayColnames <- c(displayColnames, "Others") }
  
  hmDataOrder <- displayColnames 
  
  
  tree <- read.tree(, file = "phylogenetic_tree_from_OTL.newick")
  
  # Load your heatmap data
  #Data_file_17_09_2024 <- read_excel("./Data_file_17.09.2024.xlsx")
  #Data_file_17_09_2024 <- read_excel("./Data_file_exudates_untargeted.xlsx")
  Data_file_17_09_2024 <- read_excel("Supplementary-Table-1_Targeting-the-untargeted (1).xlsx")
  
  
  #str(Data_file_17_09_2024)
  #Data_file_17_09_2024[, compound_columns]
  
  ####### Data_file_17_09_2024[Data_file_17_09_2024$`Plant species` %in% c('Arabidopsis thaliana'), ]
  ##### unique(Data_file_17_09_2024[!(Data_file_17_09_2024$`Plant species` %in% c(species_list)), 'Plant species'])
  ##### testDBframme <- data.frame(species = species_list)
  ##### testL <- as.list(unique(Data_file_17_09_2024[, 'Plant species']))
  ##### testDBframme[!testDBframme[, 'species'] %in%  testL$`Plant species`,]
  
  #Data_file_17_09_2024$`New exudation analytical method condensed`
  Data_file_17_09_2024 <- Data_file_17_09_2024 |> 
    filter(`Analytical method` %in% c(methodType)) |> 
    filter(!`First author` %in% initDB[mT, 'excludeFirstAuthor']) |> 
    filter(!`Plant species` %in% initDB[mT, 'excludeSpecies']) 
  
  
  
  # filter(!`Plant species` %in% 'Dactylis glomerata')
  
  colnames(Data_file_17_09_2024)
  # 
  # compound_columns <- c(
  #   "Carbohydrate contribution to total annotated [%]",
  #   "Fatty acid contribution to total annotated [%]",
  #   "Animo acid and peptide contribution to total annotated [%]",
  #   "Terpenoid contribution to total annotated [%]",
  #   "Shikimates and Phenylpropanoids contribution to total annotated [%]",
  #   "Alkaloid contribution to total annotated compounds [%]",
  #   "Polyketide contribution to total annotated compounds [%]",
  #   "Other contribution to total annotated compounds [%]"
  # )
  
  
  # Aggregate data by plant species and calculate the mean contribution
  aggregated_data <- Data_file_17_09_2024 %>%
    group_by(`Plant species`) %>%
    summarise_at(vars(all_of(compound_columns)), mean, na.rm = TRUE)
  
  # Calculate the number of experiments per plant species
  study_counts <- Data_file_17_09_2024 %>%
    group_by(`Plant species`) %>%
    summarise(study_count = n(), total_annotated = mean(`Sum of annotated compounds`))
  #  summarise(study_count = n(), total_annotated = sum(`Total number of annotated compounds [n]`))
  
  #aggregated_data_mat$log_total_annotated <- 1
  # Total number of annotated compounds (your actual data)
  ####total_annotated <-Data_file_17_09_2024$`Total number of annotated compounds [n]`
  #total_annotated <- c(51, 37, 632, 42, 101, 29, 75, 97, 21, 34, 14, 34, 318, 275, 21, 36, 42, 38, 90, 99, 82, 34, 132, 38, 162, 308, 22, 34, 62, 1179, 1179, 56, 157, 1179, 39, 775, 15, 60, 74, 38, 34, 34, 39, 376, 34, 52, 51, 89339, 62, 46, 7, 58, 48, 73, 72, 11, 11, 436, 222, 1636)
  
  # Merge the study counts with the aggregated data
  aggregated_data_with_study_counts <- merge(aggregated_data, study_counts, by = "Plant species")
  
  # Modify the plant species names to include the number of experiments in parentheses
  aggregated_data_with_study_counts$`Plant species` <- paste0(
    aggregated_data_with_study_counts$`Plant species`, " (", aggregated_data_with_study_counts$study_count, ")"
  )
  
  # Set row names to plant species and remove the column from the data frame
  aggregated_data_mat <- as.data.frame(aggregated_data_with_study_counts)
  rownames(aggregated_data_mat) <- aggregated_data_mat$`Plant species`
  aggregated_data_mat$`Plant species` <- NULL
  #aggregated_data_mat$study_count <- NULL  # Remove the study count column as it's now in the row names
  
  
  colnames(aggregated_data_mat) <- c(displayColnames, 'study_count', 'total_annotated')
  
  # Log transform the total annotated compounds to handle large ranges
  #### log_total_annotated <- log10(total_annotated + 1)  # Adding 1 to avoid log(0)
  #### aggregated_data_mat$log_total_annotated <- log10(total_annotated + 1)  # Adding 1 to avoid log(0)
  #aggregated_data_mat$log_total_annotated <- 1
  
  aggregated_data_mat$log_total_annotated <- log10(aggregated_data_mat$total_annotated + 1) 
  
  # Create a data frame for the annotation (total annotated compounds)
  annotation_row <- data.frame(Total_Annotated_Compounds = aggregated_data_mat$log_total_annotated)
  rownames(annotation_row) <- rownames(aggregated_data_mat)
  
  # Define appropriate color breaks for the total annotated compounds
  total_compound_breaks <- c(1, 2, 3, 4, 5, 6, 7)  # Log scale values corresponding to actual values
  
  # Create a color palette for the total annotated compounds
  annotation_colors <- list(
    Total_Annotated_Compounds = colorRampPalette(c("lightyellow", "orange", "red", "darkred"))(length(total_compound_breaks) - 1)
  )
  
  # Plot the heatmap with custom color breaks for total annotated compounds
  pheatmap(
    aggregated_data_mat[,displayColnames],
    cluster_rows = TRUE,             # Cluster plant species
    cluster_cols = TRUE,             # Cluster compound classes
    scale = "column",                # 'none', # Normalize by columns (compound classes)
    color = colorRampPalette(c("blue", "white", "red"))(50),  # Blue-white-red scale for contributions
    fontsize_row = 8,                # Font size for plant species
    fontsize_col = 8,                # Font size for compound classes
    main = "Heatmap of Compound Class Contributions by Plant Species",
    display_numbers = FALSE,         # Optionally display numbers inside the heatmap cells
    angle_col = 45,                  # Slant the column names at a 45-degree angle
    labels_row = rownames(aggregated_data_mat),  # Display plant species names with study counts
    legend_labels = "Contribution (%)",          # Adjust legend label
    annotation_row = annotation_row,  # Add the annotation for total annotated compounds
    annotation_colors = annotation_colors,  # Use the defined colors for the annotation
    annotation_legend = TRUE  # Display the annotation legend
  )
  
  
  
  
  # We will remove the OTT IDs and other unnecessary parts from the tree species names
  tree$tip.label <- tree$tip.label %>%
    str_replace_all("_ott\\d+$", "") %>%  # Remove OTT IDs at the end
    str_replace_all("_", " ") %>%         # Replace underscores with spaces
    tolower() %>%                         # Convert to lowercase to match heatmap names
    str_squish()                          # Remove any extra spaces
  
  # Step 2: Clean the species names in the heatmap data
  # Assuming the heatmap data is already in lowercase and formatted as needed
  rownames(aggregated_data_mat) <- rownames(aggregated_data_mat) %>%
    str_replace_all("_", " ") %>%         # Replace underscores with spaces
    tolower() %>%                         # Convert to lowercase
    str_squish()                                  # Remove extra spaces
  
  
  aggRownames <- str_split(rownames(aggregated_data_mat), pattern = ' ', simplify = T)
  aggRownamesA <- paste(aggRownames[,1], aggRownames[,2], sep = ' ')
  rownames(aggregated_data_mat) <- aggRownamesA
  # Ensure species names match between heatmap and tree
  common_species <- intersect(tree$tip.label, aggRownamesA)
  
  # Check if common species exist
  print(paste("Number of common species:", length(common_species)))
  print(common_species)
  
  
  # Check species names in the tree after cleaning
  print("Species names in tree:")
  print(tree$tip.label)
  
  # Check species names in the heatmap data
  print("Species names in heatmap:")
  print(rownames(aggregated_data_mat))
  
  
  # Step 3: Subset the data to include only common species
  if (length(common_species) > 0) {
    # Subset the heatmap data to include only common species
    aggregated_data_mat_subset <- aggregated_data_mat[rownames(aggregated_data_mat) %in% common_species, ]
    
    # Reorder the heatmap rows to match the order of species in the tree
    # reordered_species <- match(cleaned_tree_tips, rownames(aggregated_data_mat_subset))
    #  aggregated_data_mat_subset <- aggregated_data_mat_subset[reordered_species, ]
    
    # Step 4: Log transform the total annotated compounds (retaining your previous customization)
    #dg log_total_annotated <- log10(total_annotated + 1)  # Adding 1 to avoid log(0)
    #dg   annotation_row <- data.frame(Total_Annotated_Compounds = aggregated_data_mat_subset$log_total_annotated)
    #dg   rownames(annotation_row) <- rownames(aggregated_data_mat_subset)
    
    # Define color palette for annotations
    annotation_colors <- list(
      Total_Annotated_Compounds = colorRampPalette(c("lightyellow", "orange", "red", "darkred"))(100)
    )
    
    # Step 5: Plot the heatmap aligned with the phylogenetic tree using phylo.heatmap
    cleaned_tree_tips <- tree$tip.label
    tree_subset <- drop.tip(tree, tree$tip.label[!cleaned_tree_tips %in% common_species])  # Subset the tree
    # Ensure the order of species in the heatmap matches the tree
    ordered_heatmap <- aggregated_data_mat_subset[match(tree_subset$tip.label, rownames(aggregated_data_mat_subset)), ]
    
    tree_subset$tip.label <- paste0(str_to_sentence(tree_subset$tip.label), ' (', ordered_heatmap$study_count,   ')')
    rownames(ordered_heatmap) <- tree_subset$tip.label
    
    # dg added
    annotation_row <- data.frame(Total_Annotated_Compounds = ordered_heatmap$log_total_annotated)
    rownames(annotation_row) <- rownames(ordered_heatmap)
    
    ordered_heatmap$Zero <- NA
    # Plot the phylogenetic tree and heatmap side by side
    ##  layout(matrix(c(1,2), 1, 2, byrow = TRUE))
    ## show the regions that have been allocated to each plot
    #layout.show(2)
    ##  phylo.heatmap(tree_subset, ordered_heatmap[,  c('log_total_annotated', 'Zero')], fsize = 0.7, fsize.tip = 0.6, 
    ##                colors = colorRampPalette(c("lightyellow", "orange", "red", "darkred"))(100), # Blue-white-red scale for contributions
    ##                legend=T
    ##                )
    
    phylo.heatmap(tree_subset, ordered_heatmap[, displayColnames], fsize = 0.7, fsize.tip = 0.6,
                  standardize = F,
                  scale = "column",                # Normalize by columns (compound classes)
                  colors = colorRampPalette(c("blue", "white", "red"))(50),  # Blue-white-red scale for contributions
                  main = "Heatmap of Compound Class Contributions by Plant Species",
                  legend=T)
    
    #tips = paste(str_to_sentence(rownames(ordered_heatmap)), 'X'),
    #display_numbers = T,         # Optionally display numbers inside the heatmap cells
    #node.numbers = T)
    
    #              labels_row = str_to_sentence(rownames(ordered_heatmap)),  # Display plant species names with study counts
    #              legend_labels = "Contribution (%)",          # Adjust legend label
    #              annotation_row = annotation_row,  # Add the annotation for total annotated compounds
    #              annotation_colors = annotation_colors,  # Use the defined colors for the annotation
    #              annotation_legend = TRUE  # Display the annotation legend
    
    
  } else {
    stop("No common species found between tree and heatmap.")
  }
  
  
  
  gcPhylo <- tree_subset
  
  if (is.null(gcPhylo$edge.length)) {
    gcPhylo$edge.length <- rep(1, nrow(gcPhylo$edge)) # Assign length 1 to all branches
  }
  ultrametric_tree <- chronos(gcPhylo, lambda = 1)
  
  dhc <- as.dendrogram(ultrametric_tree)
  
  myTreeData <- dendro_data(dhc, type = "rectangle")
  
  #myTreeData$labels
  myTreeDataA <- label(myTreeData)
  myTreeDataA$dislabel <- paste0("italic('", 
                                 paste0(str_split(myTreeDataA$label, pattern = ' ', simplify = T)[,1], ' ', 
                                        str_split(myTreeDataA$label, pattern = ' ', simplify = T)[,2]), "')~'",
                                 str_split(myTreeDataA$label, pattern = ' ', simplify = T)[,3], "'")
  
  ordered_heatmap$Species <- rownames(ordered_heatmap)
  ordered_heatmap$ASp <- paste0('A', str_pad(1:nrow(ordered_heatmap), 4, pad = 0), ' ', rownames(ordered_heatmap))
  ordered_heatmap$Blank <- NA
  ordered_heatmap$Annotated <- NA
  
  ordered_heatmap$total_annotated_label <- as.integer(round(ordered_heatmap$total_annotated, digits = 0))
  
  hmData <- pivot_longer(ordered_heatmap[,c('ASp', 'Species', displayColnames, 'Blank', 'Annotated', 'total_annotated_label', 'log_total_annotated')], 
                         cols = c(displayColnames, 'Blank', 'Annotated', 'total_annotated_label', 'log_total_annotated'), names_to = "X", values_to = 'Component class contribution')
  
  hmData$X <- factor(hmData$X, 
                     levels = c(hmDataOrder, 'Blank', 'Annotated', 'total_annotated_label', 'log_total_annotated'), 
                     labels =  c(hmDataOrder, ' ', 'Mean annotated', 'total_annotated_label', 'log_total_annotated'))
  
  lsitASp <- str_pad(paste0(str_split(unique(hmData$ASp), pattern = ' ', simplify = T)[,2], ' ', 
                            str_split(unique(hmData$ASp), pattern = ' ', simplify = T)[,3], ' ', 
                            str_split(unique(hmData$ASp), pattern = ' ', simplify = T)[,4]), 
                     width = 35, pad ='-', side = 'left', use_width = T)
  hmData$ASp <- factor(hmData$ASp, levels = unique(hmData$ASp), labels = lsitASp)
  
  myTreeDataA$yw <- 0
  
  for (i in 1:nrow(myTreeDataA)) {
    myTreeDataA$yw[i] <-strwidth(myTreeDataA$dislabel[i])
  }
  myTextDispWidth <- (max(myTreeDataA$yw) * initDB[mT, 'C_TextDispWidthFactor']) + 0.5
  
  hmDataCount <- hmData |> filter(X %in% c('total_annotated_label'))
  hmDataCount$X <- 'Mean annotated'
  hmDataCol <- hmData |> filter(X %in% c('log_total_annotated'))
  hmDataCol$X <- 'Mean annotated'
  
  # Funktion für Farbverlauf erstellen
  create_gradient <- function(colors, n) {
    gradient_function <- colorRampPalette(colors)
    colors_list <- gradient_function(n)
    return(colors_list)
  }
  # Beispiel: Farbverlauf mit drei Ausgangsfarben
  colors_three <- c("lightyellow", "orange", "red")  # Die gewählten Ausgangsfarben
  maxGradColors <- 100
  gradient_three <- create_gradient(colors_three, maxGradColors+1)  # Erzeugt 100 interpolierte Farben
  hmDataCol$FillColors <- gradient_three[as.integer((hmDataCol$`Component class contribution` - 
                                                       min(hmDataCol$`Component class contribution`)) *maxGradColors
                                                    / (max(hmDataCol$`Component class contribution`) - 
                                                         min(hmDataCol$`Component class contribution`)))+1] 
  
  C_DO_DISPLAY_LEGEND_POS <- 'none'
  if (initDB[mT, 'C_DO_DISPLAY_LEGEND_IN_ALL_PLOTS']) { C_DO_DISPLAY_LEGEND_POS <- 'bottom' } 
  
  
  
  pPhyloHeat <- hmData |> filter(!(X %in% c('total_annotated_label', 'log_total_annotated'))) |> 
    ggplot(aes(X, ASp, fill= `Component class contribution`)) + 
    geom_tile(width=1) + xlab('') + ylab('') +  
    coord_cartesian(clip = "off") +  # Clipping ausschalten
    #geom_rect(aes(xmin = 4, xmax = 6, ymin = -1, ymax = 1), fill = "white", color = "black") +
    
    geom_rect(data =  data.frame(xmin = -(myTextDispWidth), ymin = 0,
                                 xmax = -(myTextDispWidth+max(segment(myTreeData)$yend)+0.5), 
                                 ymax = max(segment(myTreeData)$xend)+1), 
              mapping = aes(xmin = xmin, ymin = ymin, xmax = xmax, ymax = ymax), 
              inherit.aes = F, fill='white')+
    geom_segment(data =  segment(myTreeData), 
                 mapping = aes(x = -(y+myTextDispWidth), y = x, xend = -(yend+myTextDispWidth), yend = xend), inherit.aes = F)+
    geom_point(data =  myTreeDataA, 
               mapping = aes(x = -(myTextDispWidth-0.1), y = x), inherit.aes = F, size = 0.8, color = 'black')+
    geom_rect(data =  myTreeDataA, 
              mapping = aes(xmin = 0.75, xmax = -(yw* initDB[mT, 'C_TextBackgroundWidthFactor']), ymin = x-0.5, ymax = x+0.5), inherit.aes = F,fill = 'white')+
    geom_text(data = myTreeDataA, 
              aes(x = 0.5, y = x, label = dislabel), 
              hjust = 1, vjust = 0.5, parse = T,  inherit.aes = F) +
    
    geom_tile(data = hmDataCol,  
              mapping = aes(x = X, y = ASp), width=1.3, show.legend = F, 
              inherit.aes = 'none', fill = hmDataCol$FillColors) + 
    geom_text(data = hmDataCount, 
              mapping = aes(x = X, y = ASp, label = `Component class contribution`), size = 3, show.legend = F,
              inherit.aes = 'none') +
    
    scale_fill_gradient2(name = "Component class contribution [%]", low="white", mid="blue", high = 'darkblue', midpoint = 50, na.value = 'white') + 
    
    
    theme_pubclean () +
    theme(
      legend.position = C_DO_DISPLAY_LEGEND_POS,            # Legende unten
      legend.title = element_text(margin = margin(b = 16)),
      axis.text.y = element_blank(),
      #   axis.title.y = element_text(margin = margin(l = 20)),
      #axis.text.y = element_text(face = "italic", hjust = 1), # family = "Courier"),
      axis.ticks.length = unit(0, "pt"),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1.07),  # x-Achsen-Beschriftung drehen
      panel.grid.major = element_blank(),     # Haupt-Hilfslinien ausblenden
      panel.grid.minor = element_blank(),      # Neben-Hilfslinien ausblenden
      panel.border = element_blank(),         # Rand des Panels ausblenden
      axis.line = element_blank(),             # Achsenlinien ausblenden
      plot.margin = unit(c(0, 0.2, 0, -0.5), "cm")
    ) 
  
  #  print(pPhyloHeat)
  
  combined_plotA <- ggarrange(NULL, pPhyloHeat, heights = c(initDB[mT, 'C_TopLabelPanelHeight'], 1), 
                              ncol = 1,  labels = c(initDB[mT, 'Label'], ''), hjust = 0.001)
  
  # ggsave(paste0('./phyloheat', methodType, '.pdf'), combined_plot, device = 'pdf', width = 240, height = initDB[mT, 'FigureHeight'], units = 'mm', bg = 'white')#  print(combined_plot)
  
  ggsave(paste0('./phyloheat', methodType, '.pdf'), combined_plotA, device = 'pdf', width = initDB[mT, 'FigureWidth'], height = initDB[mT, 'FigureHeight'], units = 'mm', bg = 'white')
  ggsave(paste0('./phyloheat', methodType, '.png'), combined_plotA, device = 'png', width = initDB[mT, 'FigureWidth'], height = initDB[mT, 'FigureHeight'], units = 'mm', bg = 'white')
  #  ggsave(paste0('./phyloheat', methodType, '.eps'), combined_plotA, device = 'eps', width = initDB[mT, 'FigureWidth'], height = initDB[mT, 'FigureHeight'], units = 'mm', bg = 'white')
}
```

#Data table exports####
```{r Table exports}

#Table export for figure 1
write_xlsx(pub_data_full, "Graph1.xlsx")

#Table export for figure 2
write_xlsx(planttype_order_count_reorder, "Graph2.xlsx")

#Table export for figure 3
write_xlsx(n_per_method, "Graph3.xlsx")

#Table export for figure 4
write_xlsx(collectiontime_per_approach, "Graph4.xlsx")

#Table export for figure 5
write_xlsx(percent_annotated_noNA, "Graph5.xlsx")

#Table export for figure 5.1
write_xlsx(percent_annotated_noNA_u100, "Graph5.1.xlsx")

#Table export for figure 5.2
write_xlsx(percent_annotated_noNA_a100, "Graph5.2.xlsx")

#Table export for figure 6
write_xlsx(combined_summary, "Graph6.xlsx")
```