# Intro to Spatial Autocorrelation in R
Riley Watts
October 20th, 2024

## Introduction

In 1970, Waldo Tobler stated that the First Law of Geography is “everything is related to everything else, but near things are more related than distant things” (Tobler, 1970). We can explore how this argument applies to different variables through the concept of spatial autocorrelation (SAC), which explores how geographically proximate locations exhibit similar (or dissimilar) values for a given variable. Positive SAC will indicate that similar values are clustered together, while negative values indicate a dispersed pattern. Both patterns usually point to underlying factors that are shaping how variables are distributed across space.

This tutorial is an introduction to exploring SAC using R to visualize spatial distribution and run two key statistics that quantify variation between values at locations in proximity to each other: Moran’s I and Local Moran’s I. In recent years, spatial autocorrelation has been applied in numerous fields to explore the spatial dynamics behind socio-economic and demographic phenomena, such as income inequality and linguistic diversity (Rey & Janikas, 2006; Khan & Siddique, 2021). Our example uses 2016 Canadian census from the county of Antigonish, Nova Scotia that describes income levels and French knowledge language. Census data provides a wealth of demographic and socio-economic variables at small spatial scales, allowing us to explore spatial patterns that might not be visible at larger scales.

## Setting up your project

### Libraries

One of your first steps in most R projects will be loading libraries. These are collections of functions and datasets that broaden the range of what we can do with R. Here's an overview of the libraries we'll be using:

-   **sf** is short for “simple features” and allows us to list these features as records in a data.frame which makes handling and analyzing spatial data much easier **([an sf cheat sheet can be found here](https://github.com/rstudio/cheatsheets/blob/main/sf.pdf))**

-   **sp** offers allows us to use traditional spatial data classes such as point, lines, and polygons and provides methods to analyze them ([**learn more about sp and it’s functions here**](https://www.rdocumentation.org/packages/sp/versions/2.1-4)**;** [original paper](https://cran.r-project.org/web/packages/sp/vignettes/intro_sp.pdf)).

-   **e1071** provides functions for statistical analysis or machine learning tasks and in this case provides the **skewness()** function as part of descriptive statistics (Meyer et al., 2022).

-   **spdep** provides tools for spatial dependency and autocorrelation analysis and includes tools we'll use such as **nb2listw(), moran.test(), and localmoran()** (Bivand et al., 2013).

-   **tmap** is a mapping tool that includes functions for thematic mapping, and will allow us to visualize the data and analysis. Any functions beginning with **tm\_** are sourced from this package (Tennekes, 2018).

-   **shinyjs** allows us to use interactive tools called Shiny applications that have been built using the Shiny package in R and JavaScript. In this tutorial we'll use a tool called pallete_explorer to select palettes for our map (Dean Attail, 2021).

-   **raster** is usually used for manipulating rasters in R, but also provides the **crs()** function for matching coordinate systems that we'll be using (Hijmans, 2022).

To use these libraries, there are two steps:

1\. Install the relevant libraries using the command install.packages(“name of library”). In the code below, remove the '\#' before each package you need to install to run the code. If you've already installed these libraries, you can delete this part of the code or leave the hashtags, you only need to do this once for your computer.

2.Load the libraries using the command library(“selected library”). You’ll have to do this each time you start a new session in R.

```{r Libraries, echo=TRUE, eval=TRUE, message=FALSE, warning=FALSE}
#Install packages if not already installed:
#install.packages("knitr")
#install.packages("sf")
#install.packages("tmap")
#install.packages("spdep")
#install.packages("shinyjs")
#install.packages("e1071")
#install.packages("raster")

#Load in libraries:
library("knitr")
library("sf")
library("tmap")
library("spdep")
library("e1071")
library("raster")
library("shinyjs")
```

### Working Directory

The next key part of R set up will be making sure you’ve correctly set your working directory. This is the folder on your computer where R will look for files or save output. Think of it as the "home base" for your project. To set it, you can use this command. Replace "path/to/your/folder" with the actual location where your data is stored. You can find your current working directory with:

```{r Read in data, echo=TRUE, eval=TRUE, warning=FALSE}
#Set working directory. Make sure to replace the path. 
setwd("C:/path to your folder")

```

Before moving on it's a good idea to double-check that your working directory is where you want it. Using this command will tell you:

``` r
#Find out where you're working directory is currently
getwd()
```

    ## [1] "C:/Users/watts/Desktop/Fall 2024 Classes/Geog 418/Assignment 3 Riley Watts"


### Reading and Cleaning Data

Now we’re ready to load our data into R. For this project we’re working with two files:

-   **Shapefile:** This contains the geographic boundaries for the census areas, and we'll load it as an sf (simple features) object.

-   **CSV**: This contains the attributes for each area. In this case they are census statistics (like income and French knowledge) for each area, and we'll load it as a dataframe.

The "lda_000b16a_e.shp" shapefile delineates the census dissemination area boundaries that we'll be using in this analysis as well as GEO UIDs for each. The dissemination areas (DA) used in this analysis represent small geographic regions, making them appropriate for detecting localized patterns of clustering or dispersion in socio-economic attributes like income and language knowledge as analyzed here.

The "census_data.csv" file has 11 columns, with attribute data including: the GEO UID, Province code, Province name, CD code, CD name, DA name, Population, Land area, Median total income, Income Sample Size, French Knowledge, and Language Sample Size.

Here's how to do it:

``` r
#From the working dir read in the csv
csv <- read.csv("census_data.csv") 

#Data source is the working dir (where the layer is), layer is the name of the file 
shp <- st_read("lda_000b16a_e.shp")
```

    ## Reading layer `lda_000b16a_e' from data source 
    ##   `C:\Users\watts\Desktop\Fall 2024 Classes\Geog 418\Assignment 3 Riley Watts\lda_000b16a_e.shp' 
    ##   using driver `ESRI Shapefile'
    ## Simple feature collection with 56589 features and 22 fields
    ## Geometry type: MULTIPOLYGON
    ## Dimension:     XY
    ## Bounding box:  xmin: 3689439 ymin: 659338.9 xmax: 9015737 ymax: 5242179
    ## Projected CRS: PCS_Lambert_Conformal_Conic
    
Before we analyze the data, cleaning it will allow us to identify the information we need, remove anything we don’t, and make sure the format is appropriate for our analysis.

To make sure information is identifiable, we’ll create a vector of the columns in our csv file with appropriate names and apply these new names to our a dataframe. Next, we’ll remove data we don’t want in the analysis. The code below creates a new field by adding a column assigned the number of ID characters in the GEO UID and removes rows with a GEO UID that isn't 8 characters by creating a new dataframe (csv_clean) with only rows that have the appropriate value.

To map our census data we’ll have to merge it with the geographic features in our shapefile using the **merge()** command. We can do this using the common “GEO UID” attribute. From the merged dataframe, we’ll narrow our analysis down to a single city by creating a **subset()** of the data.

For raw count datasets such as the census file used, standardizing by converting them to rates will allow us to compare across regions with different populations and scales. Below, we convert our “French Knowledge” and “Language Sample Size” attributes into a new Percent French Speakers attribute.

```{r Clean data, echo=TRUE, eval=TRUE, warning=FALSE}
#New column names
cols <- c("GEO UID", "Province code", "Province name", "CD code",
        "CD name", "DA name", "Population", "Land area", 
        "Median total income", "Income Sample Size", "French Knowledge", 
        "Language Sample Size")

#Apply those names to dataframe
colnames(csv) <- cols

#Add column to count number of ID charactors
csv$len <- nchar(csv$`GEO UID`)

#Remove IDs with less than 8 numbers
csv_clean <- subset(csv, csv$len == 8)

#Merge spatial and aspatial data
census_DAs <- merge(shp, csv_clean, 
                    by.x = "DAUID", 
                    by.y = "GEO UID", 
                    all.x = TRUE)

#Subset for Antigonish 
Municp <- subset(census_DAs, census_DAs$CDNAME == "Antigonish")

#Convert to rate
Municp$PercFrench <- (Municp$`French Knowledge`/Municp$`Language Sample Size`)*100
```

The final step to make sure our analysis is accurate will be removing any missing (NA) values that will be irrelevant and change our results. To make sure we’re only working with areas that have complete data for both of our variables, we’ll select only rows with values from each column using the **!is.na()** function.

```{r NA Remove, echo=TRUE, eval=TRUE, warning=FALSE}
#Remove Income NA
Income_noNA <- Municp[which(!is.na(Municp$`Median total income`)),]

#Remove French NA
French_noNA <- Municp[which(!is.na(Municp$`PercFrench`)),]
```

### Descriptive Statistics

After cleaning the data it's a good idea to understand their characteristics using some basic descriptive stats. What's the average value? How varied is the data? Is it skewed?

To do this we'll run four functions for each variable, **mean()**, **sd()**,**skewness()**, and **kurtosis()**. Assigning the output of each it's own name and creating **data.frame()** will allow us to display it as a table using the **kable()** function. Note that we've also rounded each value to 2 decimals using the **round()** function.

``` r
#Calculate descriptive stats for Income
meanIncome <- mean(Income_noNA$`Median total income`)
stdevIncome <- sd(Income_noNA$`Median total income`)
skewIncome <- skewness(Income_noNA$`Median total income`)
kurtIncome <- kurtosis(Income_noNA$`Median total income`)

#Calculate descriptive stats for French
meanFrench <- mean(French_noNA$`PercFrench`)
stdevFrench <- sd(French_noNA$`PercFrench`)
skewFrench <- skewness(French_noNA$`PercFrench`)
kurtFrench <- kurtosis(French_noNA$`PercFrench`)

#Create dataframe for display in table
data <- data.frame(Variable = c("Income", "French Language"),
                   Mean = c(round(meanIncome,2), round(meanFrench,2)),
                   StandardDeviation = c(round(stdevIncome,2), round(stdevFrench,2)),
                   Skewness = c(round(skewIncome,2), round(skewFrench,2)), 
                   Kurtosis = c(round(kurtIncome,2), round(kurtFrench,2)))

#Produce table
kable(data, caption = paste0("Descriptive statistics for selected ", 2016 , "census variables"))
```

| Variable        |     Mean | StandardDeviation | Skewness | Kurtosis |
|:----------------|---------:|------------------:|---------:|---------:|
| Income          | 33124.00 |           4923.24 |    -0.57 |     0.08 |
| French Language |    10.24 |              8.47 |     2.40 |     6.72 |

Descriptive statistics for selected 2016census variables


Measures of central tendency and dispersion like the ones in the table above are crucial for understanding how well the models run by statistics such as Moran's I will accurately represent the data. For example, Moran's I assumes normally distributed underlying data. In that case, the large skewness (2.40) and kurtosis (6.72) of our french language dataset may lead to exaggerated autocorrelation values (Griffith, 2003). In contrast, the low skewness and kurtosis of the median income mean Moran's I will likely be more accurate. We'll expand on this further in the analysis section.

## Analysis

Now we have relevant data in an appropriate format and we understand the general tendencies of our variables, so we're ready to start our analysis.

### Making your map

The first part of the code below opens a tool called **palette_explorer()** which will allow you to select the appropriate colour pattern from ColorBrewer or Virdis according to the number of classes and the type of data pattern you're working with (sequential, categrical, or diverging). At the bottom of the window select "tmap layer function code" under code generator to see an example of the appropriate code for each set of colour options and adjust the next part of the code accordingly.

There are two tmap modes you can choose from. To create an interactive html map that will allow you to zoom in and view the data on a finer scale and compare it against a basemap, remove the hashtag in front of **tmap_mode("view")**. In this case, the **tm_layout()** parameters below aren't necessary. If you want to create a static map that can be printed as a png, leave the hashtag or delete that section of code.

There are three sections to mapping each variable using **'tmap'**: - **tm_shape()**: This function loads the spatial object that we'll be mapping. - **tm_polygons()**: This function has a set of parameters for filling our polygons. We can select: the column with the variable that we'd like to map,the title of the map, the style of classification, color palette (using the code we got from palette_explorer), the transparency of the borders (denoted by .alpha), and the color of any NA values. - **tm_layout()**: This function allows us to modify the layout of our map. In this case we've made sure the legend will appear in the top left of the map to best fit the negative space left by the Anitgonish area.

Printing the maps out next to each other will allow us to compare them more easily, so we can use **tmap_arrange()** to place them next to each other, formatting it similar to how we'd format a table.

```{r StudyArea, echo=TRUE, eval=TRUE, warning=FALSE, fig.cap="Antigonish census dissemination areas showing median total income (left) and percentage of respondants with knowledge of french (right)."}
#Choose a pallete
tmaptools::palette_explorer() #Tool for selecting pallettes

#Create an interactive map
#tmap_mode("view")

#Map median Income
map_Income <- tm_shape(Income_noNA) + 
  tm_polygons(col = "Median total income", 
              title = "Median total income", 
              style = "jenks", 
              palette = "PuRd", n = 6,
              border.alpha = 0,
              colorNA = "grey") +
  tm_layout(legend.position = c("LEFT", "BOTTOM"))

#Map French Knowledge
map_French <- tm_shape(French_noNA) + 
  tm_polygons(col = "PercFrench", 
              title = "Percentage with \n French Knowledge", 
              style = "jenks", 
              palette = "Blues", n = 6,
              border.alpha = 0,
              colorNA = "grey") +
  tm_layout(legend.position = c("LEFT", "BOTTOM"))

#Print maps side by side
tmap_arrange(map_Income, map_French, ncol = 2, nrow = 1)

```
![variablemapsstatic](https://github.com/user-attachments/assets/555aa9e1-7a6b-419c-9ccc-7d364f1e2f4c)

## Neighbourhood matrix

You can imagine a weighted neighbourhood matrix as a table that will tell us how close or connected different areas are to each other. Each cell in the table is assigned a number value that describes how close one area is to another (Cliff & Ord, 1981). Areas that are next to each other get a higher weight while areas that are farther apart get lower weights or zeros. Defining and weighing a spatial connection between a location/area and its neighbours is the first step to testing out Tobler's Law and beginning our analysis of spatial autocorrelation.

A key step in setting up this matrix is deciding what will be counted as a neighbour when we assign the values. In the code chunk below we explore two options, each named after the potential movements of a chess piece.

-   **Queen's Weight**: In a queen weighting matrix, polygons will be considered neighbours as long as they share a single boundary point. This is easiest to imagine on a square chess board, where the queen can move diagonally across squares that only meet at the corner as well as squares that share a side.

-   **Rook's Weight**: In a rook weighting matrix, polygons must share more than one boundary point, sharing a side. Imagining this on a square chess board, rooks can only move in a straight horizontal or vertical line where the squares share sides.

Since we've loaded the **'spdep'** package, we can use the **poly2nb()** function to define neighbours (Bivand et al., 2013). The default for this function is queen weighting, so we don't need to add anything but our desired spatial dataset as the variable. To change to rooks weighting, we can add 'queen = FALSE'. Our next step will be to visualize these relationships, so we can use **nb2lines()** from the same package to convert the neighbours list it created into a LINESTRING object, as long as we assign it a list ('.nb') and coordinates.

To assign it coordinates we can use the **st_coordinates()** function from the **'sf'** package to extract the coordinates of our object, in this case the '*noNA'* polygons*.* Using **st_centroid()**, we can calculate the geometric center (centroid) of those polygons to ensure the lines are mapped clearly.

Finally, we should assign our output a Coordinate Reference System (CRS) to ensure it aligns correctly on our map. The **crs()** function is part of '**sp**' and we can use it to ensure that our '.net' line object has the same CRS as our polygon features.

```{r Neighbours, echo=TRUE, eval=TRUE, warning=FALSE}

#Income Neighbours - Queens weight
Income.nb <- poly2nb(Income_noNA)
# Use st_coordinates to get the coordinates
Income.net <- nb2lines(Income.nb, coords=st_coordinates(st_centroid(Income_noNA)))
crs(Income.net) <- crs(Income_noNA)

#Income Neighbours - Rooks weight
Income.nb2 <- poly2nb(Income_noNA, queen = FALSE)
Income.net2 <- nb2lines(Income.nb2, coords=st_coordinates(st_centroid(Income_noNA)))
crs(Income.net2) <- crs(Income_noNA)

#French Neighbours - Queens weight
French.nb <- poly2nb(French_noNA)
French.net <- nb2lines(French.nb, coords=st_coordinates(st_centroid(French_noNA)))
crs(French.net) <- crs(French_noNA)

#French Neighbours - Rooks weight
French.nb2 <- poly2nb(French_noNA, queen = FALSE)
French.net2 <- nb2lines(French.nb2, coords=st_coordinates(French_noNA))
crs(French.net2) <- crs(French_noNA)

```

Returning to the **'tmap'** package, we can visualize which polygons have been included or excluded as neighbours based on each method. This time, **tm_shape()** is used to overlay multiple spatial objects, the census boundaries from the 'Income_noNA' data and the '.net' connections drawn between polygons.

There's no polygons function this time since our focus is the **tm_borders()** of the census areas and the **tm_lines** that connect the areas if they are counted as neighbours. The other addition for these maps is adjusting the line weight using '**lwd ='**. The code below creates one map for queen weight, one weight for rook, and one that combines the lines from both, using different line weights to differentiate. We don't need to repeat the process with the French data since they use the same boundaries.

Finally, we can arrange our maps next to each other like we did before, this time creating three columns for the three maps.

```{r Neighboursmap, echo=TRUE, eval=TRUE, warning=FALSE, fig.cap="Kamloops census dissemination areas showing median total income neighbours queens weight (left)  rooks weight (middle) and the combination of the two (right)."}

#Make queens map
IncomeQueen <- tm_shape(Income_noNA) + tm_borders(col='lightgrey') + 
              tm_shape(Income.net) + tm_lines(col='red', lwd = 0.1, alpha = 0.75)

#Make rooks map
IncomeRook <- tm_shape(Income_noNA) + tm_borders(col='lightgrey') + 
              tm_shape(Income.net2) + tm_lines(col='blue', lwd = 0.1, alpha = 0.7)

#Make combined map
IncomeBoth <- tm_shape(Income_noNA) + tm_borders(col='lightgrey') + 
               tm_shape(Income.net) + tm_lines(col='red', lwd = 0.25, alpha = 0.7) +
               tm_shape(Income.net2) + tm_lines(col='blue', lwd = 0.1, alpha = 0.7)

#Print maps in a three pane figure
tmap_arrange(IncomeQueen, IncomeRook, IncomeBoth, ncol = 3, nrow = 1)

```
![neighboursmapsstatic](https://github.com/user-attachments/assets/1bd606b1-bc59-44a6-9f0c-2e3e4e58764b)

Once neighbours have been defined, there are there are three different approaches we can take to assigning the weights, styles: "B", "W", and "C". The simplest approach is "B" or "Binary", where there's only two possible values, 1 or 0: 1 if it's a neighbour, 0 if it's not. However, in cases where features have an unequal number of neighbours, you may want to standardize the total weight of their neighbours by making the weight of each neighbour for a feature add up to 1 by using the "W" style. For example, if a neighbour is one of four for a feature, it will be assigned a weight of 0.25. If it is one of two, it will weigh 0.5. the final option is global standardization, "C", where we would assign the same value to all neighbours in the dataset based off the total number of neighbours. Since census boundaries are irregular in size and shape and their distribution is outside of our control, we use row standardization in the code below (Bivand et al., 2013).

To create our weights matrix, we'll us the **nb2listw()** function from the **'spdep'** library. This function takes the neighbourhood links we created earlier (.nb) and assigns the weights based off our chosen style. **For our analysis we've chosen Queen weighting because...** In the case that an area has no neighbours, we'll set the **zero.policy = TRUE** to ensure that the process runs smoothly.

To view a sample of the weight distributions, we'll use **head(Income.lw[["weights"]])** to extract entries from the matrix, and **[c(1:3)]** to specify that we just want to see the first three.

``` r
#Create Income weights matrix
Income.lw <- nb2listw(Income.nb, zero.policy = TRUE, style = "W")

#Create French weights matrix
French.lw <- nb2listw(French.nb, zero.policy = TRUE, style = "W")

#
head(Income.lw[["weights"]])[c(1:3)]
```

    ## [[1]]
    ## [1] 0.3333333 0.3333333 0.3333333
    ## 
    ## [[2]]
    ## [1] 0.5 0.5
    ## 
    ## [[3]]
    ## [1] 0.25 0.25 0.25 0.25

## Global Moran’s I

Our next step is to calculate a Global Moran's I statistic using the weights and neighbours that we defined. This is a global statistic that will indicate whether similar values in our variables are clustered together or spread out randomly (Anselin, 1995). It works by multiplying the spatial weight by both the differences between the mean value of that variable and a single feature and the difference between the mean and the neighbours of that feature. In the equation below, the mean value of our variable is denoted by $x$, the value of the feature of interest by $x_i$ and its neighbours values by $x_j$. To standardize this value so we can compare and infer across datasets, the equation uses a denominator. As a result, we can infer from looking at a high $I$ value that the data is positively spatially correlated and from a low value that it is negatively spatially correlated.

$$
I = \frac{\sum_{i=1}^n\sum_{j=1}^nW_{i,j}(x_i - \bar{x})(x_j - \bar{x})}{(\sum_{i=1}^n\sum_{j=1}^nW_{i,j})\sum_{i=1}^n(x_i - \bar{x})^2}
$$

To calculate the $I$ using R, we'll use the function **moran.test()**, selecting our variables from the chosen dataset, putting in the weights matrix we just selected as a parameter, and setting the 'zero.policy = TRUE' again to ensure it runs smoothly in zones where there are no neighbours. To make sure we have statistically significant results, we can't just calculate the observed test statistic. We have to perform a Z-test that will compare it to a theoretical value we would expect from a random distribution and take into account the variance in our values. The expected random value is calculated using the following equation, where $n$ is the number of values:

$$
E(I) = \frac{-1}{(n-1)}
$$

To extract these calculations from our test, the below code assigns the $I$ value to 'mI', our expected random value to 'eI', and the variance of $I$ to 'var'.

```{r Global Morans I, echo=TRUE, eval=TRUE, warning=FALSE}
#Calculate Global Moran's I for Income
miIncome <- moran.test(Income_noNA$`Median total income`, Income.lw, zero.policy = TRUE)

#Extract Global Moran's I results for Income
mIIncome <- miIncome$estimate[[1]]
eIIncome <- miIncome$estimate[[2]]
varIncome <- miIncome$estimate[[3]]

#Calculate Global Moran's I for French
miFrench <- moran.test(French_noNA$PercFrench, French.lw, zero.policy = TRUE)

#Extract Global Moran's I results for French
mIFrench <- miFrench$estimate[[1]]
eIFrench <- miFrench$estimate[[2]]
varFrench <- miFrench$estimate[[3]]
```

With these variables, our results are $E(I)$ of -0.0256 for French, an $I$ of 0.176, and a variance of 0.000764. For income, an $E(I)$ of -0.00278, an $I$of 0.0.381, and a variance of 0.011. At first glance, we can see that for both datasets the observed *I* is higher than that expected random *I*, meaning the distribution is leaning towards clustered.

To understand what a relative low or high value is for this particular dataset, it is helpful to calculate the range of possible $I$ values we could get. Unlike the theoretical range of Moran's $I$ from -1 to 1 as with some simpler correlation coefficient, instead it depends on the values in the weight matrix, and data doesn't often exhibit perfect spatial patterns(Anselin, 1995).

We can define our own function to calculate that range by using **function()** and defining its parameter as **.lw**, which is expected to be a spatial weights list object. Then we'll use the **listw2mat()** function from **'spdep'** to convert our spatial weights list into a matrix format **wmat**, which will make it easier to perform further calculations. We also need to make sure that the matrix is balanced, i.e., the weight being assigned to feature A as a neighbour of feature B is the same when it's the other way around, smoothing out irregularities that come from imperfectly balanced real-world data. This is done by adding the matrix we just created to its transpose and dividing it by 2 with **t()/2**. Next we can use the **eigen()** function, extract its values, and calculate the **range()** of these values. These values describe the overall spatial pattern and structure of the matrix and their range will become our Moran's $I$ range. To ensure the function outputs this range, we'll use the **return()** function.

```{r Global Morans Range, echo=TRUE, eval=TRUE, warning=FALSE}
#Function to calculate the range of global Moran's I
moran.range <- function(lw) {
  wmat <- listw2mat(lw)
  return(range(eigen((wmat + t(wmat))/2)$values))
}

#Calculate the range for the Income variable
range <- moran.range(Income.lw)
minRange <- range[1]
maxRange <- range[2]
```

The minimum for Income is -0.581 and the maximum is 1.02. Rerunning the code for French would give us an almost identical range because the number of features that defines the range is the same (Anselin, 1995).

Revisiting the values we extracted from our initial test, we can now calculate a Z-score to determine if our results are statistically significant. We'll perform a Z-test using the following formula:

Our null hypothesis will be that the data is randomly distributed, and our alternative hypothesis will be that the data is spatially autocorrelated. If we want to reject our null hypothesis, we will have to surpass a confidence threshold. Here we'll use an $\alpha$ value of 0.05, and if our Z-score falls above 1.96 or below - 1.96, we can say it is significantly spatially autocorrelated. A value greater than +1.96 would imply positive autocorrelation, and a value less than -1.96 would imply negative autocorrelation.

The equation for this Z-test is:

$$
Z = \frac{I- E(I)}{\sqrt{V(I)}}
$$

Translated into R, that becomes:

```{r Global Morans ZScore, echo=TRUE, eval=TRUE, warning=FALSE}
#Calculate z-test for Income
zIncome <- (mIIncome - eIIncome) / (sqrt(varIncome))

#Calculate z-test for French
zFrench <- (mIFrench - eIFrench) / (sqrt(varFrench))
```

The resulting Z-scores of 3.97 for income and 2.31 for french language knowledge indicate that both demographics are significantlly postively spatially autocorrelated, or clustered. Each returning a p-value of \< 0.01. However, revisiting our descriptive statistics, we know that the skewness and outliers (indicated by kurtosis) of the french dataset are likely to inflate autocorrelation values for Moran's *I*. In particular, as a global statistic that is providing a singular value to describe the whole area, Moran's *I* will be heavily influenced by outliers.

## Local spatial autocorrelation

The global Moran's $I$ statistic gave us a single value to describe the overall tendency in the distribution of each variable, but what if we want to identify neighbourhoods that are particularly isolated in their language knowledge or median income? Or maybe understand where the most homogenous areas of the data are? Especially in the case of the french language set, we know that there may be many areas of no or negative spatial autocorrelation that are masked by the high outliers when a global statistic is run. In this case we can use what's called a local autocorrelation to identify hotspots with high values, coldspots with low values, and outliers where the values are not in-line with their neighbours. This is the Local Moran's $I$ test, also called Local Indicators of Spatial Association (LISA), developed by Anselin (1995).

Because we are still asking whether features and their neighbours are alike, you'll notice many of the same features in the Local Moran's $I$ that were in the global equation. It's the arrangement of these features that makes the difference.

$$
I_i = \frac{x_i - \bar{x}}{S_i^2}\sum{_{j=1}^n}W_{i,j}(x_j - \bar{x})\space \space where \space \space S_i^2 = \frac{\sum_{i=1}^n (x_i - \bar{x})^2}{n-1} 
$$

As with our last test, **'spdep'** provides a convenient **localmoran()** fuction that only requires us to define our variable and the weight matrix as parameters. We should also still extract the resulting values, which this time include Z- and P-values as well as the $I$, Expected $I$, and variance.

```{r Local Morans I, echo=TRUE, eval=TRUE, warning=FALSE}
#Calculate LISA test for Income
lisa.testIncome <- localmoran(Income_noNA$`Median total income`, Income.lw)

#Extract LISA test results for Income
Income_noNA$Ii <- lisa.testIncome[,1]
Income_noNA$E.Ii<- lisa.testIncome[,2]
Income_noNA$Var.Ii<- lisa.testIncome[,3]
Income_noNA$Z.Ii<- lisa.testIncome[,4]
Income_noNA$P<- lisa.testIncome[,5]

#Calculate LISA test for Income
lisa.testFrench <- localmoran(French_noNA$PercFrench, French.lw)

#Extract LISA test results for Income
French_noNA$Ii <- lisa.testFrench [,1]
French_noNA$E.Ii<- lisa.testFrench [,2]
French_noNA$Var.Ii<- lisa.testFrench [,3]
French_noNA$Z.Ii<- lisa.testFrench [,4]
French_noNA$P<- lisa.testFrench [,5]
```

We can revisit our intro to **'tmaps'** to visualize and better understand our test results, this time we will add a title and manual breaks in the polygon section, and add a compass and scale bar.

```{r MappingLocalMoransI, echo=TRUE, eval=TRUE, warning=FALSE, fig.cap="Laval census dissemination areas showing LISA z-scores for median total income (left) and percentage of respondants with knowledge of french (right)."}

#Map LISA z-scores for Income
map_LISA_Income <- tm_shape(Income_noNA) +
  tm_polygons(col = "Z.Ii",
              title = "Local Moran's I Z-Scores",
              style = "fixed",
              border.alpha = 0.1,
              midpoint = NA,
              colorNA = NULL,
              breaks = c(min(Income_noNA$Z.Ii),-1.96,1.96,max(Income_noNA$Z.Ii)),
              palette = "-RdBu", n = 3)+
  tm_compass(position=c("right", "top"))+
  tm_scale_bar(position=c("right", "bottom"))+
  tm_legend(position = c("left", "top"))

#Map LISA z-scores for French
map_LISA_French <- tm_shape(French_noNA) +
  tm_polygons(col = "Z.Ii",
              title = "Local Moran's I Z-Scores",
              style = "fixed",
              border.alpha = 0.1,
              midpoint = NA,
              colorNA = NULL,
              breaks = c(min(French_noNA$Z.Ii),-1.96,1.96,max(French_noNA$Z.Ii)),
              palette = "-RdBu", n = 3)+
  tm_compass(position=c("right", "top"))+
  tm_scale_bar(position=c("right", "bottom"))+
  tm_legend(position = c("left", "top"))

#Plot maps in a 2 pane figure
tmap_arrange(map_LISA_Income, map_LISA_French, ncol = 2, nrow = 1)
```
![LISAmapsstatic](https://github.com/user-attachments/assets/8eb60a38-f9b0-4a13-a23e-b024c6917953)

The results show that there are a few area of particular clustering where areas with similar incomes are grouped together as well as a couple of census areas that have negative spatial autocorrelation with their neighbours, being surrounded by much different values. There are fewer areas of significant clustering for French knowledge.

### Moran Scatter Plots

While maps can show us the spatial distribution of the values, scatter plots can visualize the quantitative relationship between a variable and its spatial lag, also known as the average of the neighbouring values (Bivend et al., 2013). They also give us a way to visualize trends in the data on multiple scales, highlighting both global regressions and the distribution of individual points.

Let's start with the code and then we'll walk through the result.

Returning to the **'spdep'** package, we can use the **moran.plot()** function to generate a scatter plot of Moran's *I* values with just the spatial variable and a spatial weights list. As with before we'll set the zero policy to true to handle areas with no neighbours.

The spatial check ('spChk =') option can help identify errors such as data points with mismatched neighbours, improper connections, or no corresponding entries. If we were completely confident in our data structure or wanted to save computation time on a huge dataset, we could set it to FALSE and disable all checks. Setting it to NULL will perform basic checks to ensure that spatial relationships are reasonably aligned without compromising performance. We'll use that option here since we have already explored the data through other analyses.

The other parameters are for presentation. We'll set labels to NULL to avoid crowding the scatter plot with labels for each data point, set our axis labels to represent the measured variables, and run the function without any additional messages by setting quiet to NULL.

```{r MoransIScatter, echo=TRUE, eval=TRUE, warning=FALSE, fig.cap= "Moran's I scatter plot for median total income."}
#Create Moran's I scatter plot for Income
moran.plot(Income_noNA$`Median total income`, Income.lw, zero.policy=TRUE, spChk=NULL, labels=NULL, xlab="Median Total Income ($)", 
           ylab="Spatially Lagged Median Total Income ($)", quiet=NULL)
```

Here's a breakdown of how to interpret the resulting plot:

-   **Axes:** the x-axis represents the values of the variable for each census tract while the y-axis represents the spatial lag of this variable

-   **Dashed lines:** These lines indicate the mean value, with the distance from a point on the plot to the vertical line representing the distance between a census tracts value and the mean and its distance to the horizontal line representing the difference between the spatial lag and the mean.

-   **Quadrants:** The mean lines divides the plot into quadrants that allow us to quickly understand the surroundings of each data point. Points in the bottom left and top right represent relatively clustered areas, with both the tract and its neighbours having values above (top right) or below (bottom left) the mean. In contrast, relative outliers are shown in the top left and bottom right, where the tract and its neighbours move away from the mean in opposite directions. We can see in the plot above that most data points appear in the low-low or high-high quadrants, with only five outliers in the high-low or low-high areas.

-   **Regression line:** This line represents a linear regression analysis between the variable and its spatially lagged values (its neighbours). A positive slope indicates a positive spatial autocorrelation/clustering, moving through the bottom left and top right quadrants that have high-high and low-low values respectively. A negative slope will indicate a negative spatial autocorrelation/dispersion, moving through the top left and bottom right quadrants with high-low and low-high values. We can see that the plot above has a positive slope, indicating positive auto-correlation.

-   **Diamonds:** The diamonds represent statistically significant results.

We can conclude that median income in Antigonish is positively spatially autocorrelated, with a few outlier areas that are different than their neighbours. This analysis identifies one statistically significant outlier with a high median income tract surrounded by lower incomes, and one significant instance of clustering, with a group of low-income census tracts.

Let's repeat this analysis for our French knowledge variable.

```{r MoransIScatter2, echo=TRUE, eval=TRUE, warning=FALSE, fig.cap= "Moran's I scatter plot for percentage of respondants with knowledge of french."}
#Create Moran's I scatter plot for French
moran.plot(French_noNA$PercFrench, French.lw, zero.policy=TRUE, spChk=NULL, labels=NULL, xlab="Respondants with knowledge of French (%)", 
           ylab="Spatially Lagged knowledge of French (%)", quiet=NULL)
```

Again, the impact of skewness and kurtosis in the French dataset on our analysis is much clearer to see in this plot than in the maps we produced. The regression identifies a positive spatial autocorrelation, with two statistically significant instances of clustering where both the census tract and its neighbours have a high percentage of people with french knowledge. There is also one outlier, with a low percentage of people with french knowledge surrounded by census tracts with higher percentages. We can see that while there is a positive regression, there is also a quite even spread of points in the positive and negative correlation quadrants.

## Summary

Our analysis of spatial autocorrelation (SAC) for both income and French language knowledge in Antigonish revealed significantly positive SAC for both variables.

For income, this pattern likely reflects broader socio-economic trends due to factors such as housing, employment opportunities, and access to resources (Khan & Siddique, 2021).

For french, this pattern may be linked to historical migration patterns, cultural factors, or social networks that reinforce linguistic clustering (Kiernan, 2014). However, it is important to note that the skewness and kurtosis of the French language dataset impacted these results. The dataset's skewness, indicating a large concentration of low percentages of French speakers, and its kurtosis, highlighting the presence of extreme values, suggest that a few areas have disproportionately high percentages of French speakers compared to the rest. These characteristics can inflate the degree of observed autocorrelation, making it crucial to consider data transformations or alternative methods to address non-normality.

### Next Steps in Spatial Analysis

While Moran’s I and Local Moran’s I (LISA) provided valuable insights, there's a few more advanced techniques that could be a next step for us to explore these patterns:

1.  **Geary’s C** is a global spatial autocorrelation measure that focuses on local contrasts rather than global averages like Moran's *I* does. As a result, it is more sensitive to local variations than Moran's I and could complement this analysis by offering additional insights into local patterns (Geary, 1954).

2.  **Getis-Ord Gi**\*: This hotspot analysis technique could be applied to identify regions with statistically significant clusters of high or low values of income and French language knowledge, helping to pinpoint specific areas of interest (Getis & Ord, 1992).

3.  **Multiscale Geographically Weighted Regression (MGWR)**: A more sophisticated approach that allows for modeling the relationships between variables at different spatial scales. MGWR could reveal how the relationship between income and French language knowledge changes across different regions (Fotheringham et al., 2017).

4.  **Spatial Filtering Techniques**: If skewness continues to impact the results, applying spatial filtering methods, such as Eigenvector Spatial Filtering (ESF), could help remove spatial autocorrelation from the dataset, allowing for more accurate modeling of non-spatial processes (Griffith, 2003).

## References

Anselin, L. (1995). Local Indicators of Spatial Association—LISA. Geographical Analysis, 27(2), 93-115. 

Bivand, R., Pebesma, E., & Gómez-Rubio, V. (2013). Applied spatial data analysis with R, Second edition. Springer, NY. https://asdar-book.org/

Cliff, A., & Ord, J. (1981). Spatial Processes: Models and Applications. Pion. 

Dean Attali (2021). shinyjs: Easily Improve the User Experience of Your Shiny Apps in Seconds. R package version 2.1.0. https://CRAN.R-project.org/package=shinyjs

Fotheringham, A. S., Yang, W., & Kang, W. (2017). Multiscale Geographically Weighted Regression (MGWR). Annals of the American Association of Geographers, 107(6), 1247-1265. 

Geary, R. C. (1954). The Contiguity Ratio and Statistical Mapping. The Incorporated Statistician, 5(3), 115-145. 

Getis, A., & Ord, J. K. (1992). The analysis of spatial association by use of distance statistics. Geographical Analysis, 24(3), 189-206. 

Griffith, D. A. (2003). Spatial Autocorrelation and Spatial Filtering: Gaining Understanding Through Theory and Scientific Visualization. Springer Science & Business Media. 

Hijmans, R.J. (2022). raster: Geographic Data Analysis and Modeling. R package version 3.5-15. https://CRAN.R-project.org/package=raster

Khan, M.S. & Siddique, A.B. (2021). Spatial Analysis of Regional and Income Inequality in the United States. Economies 9(4), 159. 

Meyer, D., Dimitriadou, E., Hornik, K., Weingessel, A., & Leisch, F. (2022). e1071: Misc Functions of the Department of Statistics, Probability Theory Group (Formerly: E1071), TU Wien. R package version 1.7-9. https://CRAN.R-project.org/package=e1071

Tennekes, M. (2018). "tmap: Thematic Maps in R." Journal of Statistical Software, 84(6), 1-39. doi:10.18637/jss.v084.i061

Tobler, W. R. (1970). A computer movie simulating urban growth in the Detroit region. Economic Geography, 46(sup1), 234-240.
