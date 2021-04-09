library(animation)
library(gtools)
library(fBasics) # for Heaviside function
library(doParallel)
library(gplots)

seasonCols <- c(
	"Spring" = "#9FCE62", 
	"Calving" = "#FFE548", 
	"postCalving" = "#FFA251",
	"Summer" = "#DC4A45",
	"Fall" = "#8E4D1D",
	"Winter" = "#618CCF")


###############################################################################################
# Parameter definition
###############################################################################################

# Define annual cycle
DOY <- c(1:365)

# What habitat (range) do caribou occupy throughout the year?
# See Fig. 11 of https://www.enr.gov.nt.ca/sites/enr/files/resources/draft_-_caribou_range_assessment_and_technical_information.pdf

# Used for temperature data from MERRA for the different ranges
bouRange <- c(
	rep("Winter", 109), # January - April 19
	rep("Spring", 43), # 20 April - June 1
	rep("Calving", 27), # 2 June - 28 Jun (includes post-calving migration)
	rep("Summer", 70), # 29 Jun - 6 September
	rep("Fall", 85), # 7 September - 30 November
	rep("Winter", 31) #  1 December - 31 December
)

if(length(bouRange) != length(DOY)) warning("\n\n\n\n\n***********\n\n\n\nbouRange not the right length. Seasons don't add up.\n\n\n\n\n***********\n\n\n\n")

# Breeding date, when animals move up a class, is June 7 (DOY = 158)
breedDOY <- as.numeric(strftime(as.Date("1985-06-07"), "%j"))

# L4 resume development at the start of spring migration
L4startDOY <- as.numeric(strftime(as.Date("1985-04-20"), "%j"))

#------------------------------------------------------------------------------
# Movement speed
#------------------------------------------------------------------------------
# What is the migration speed of caribou throughout the year

# 1) Based on data
# bouMove <- read.csv("parameterization/dailyMovementRates.csv")
# # How long is the migration based on the Gunn daily movement rates?
# sum(round(bouMove$movementRate)) # 3089 km

# Interpolated for seasons.
migSpeed <- 14

bouMove <- c(
	rep(0, (109 - 11)),# winter
	rep(migSpeed, 11 + 43), #spring
	rep(0, 15), # calving
	rep(migSpeed, 27-15), # post-calving migration
	rep(0, 70), # summer
	rep(migSpeed, 96),# Fall migration
	rep(0, 20)) # winter

###############################################################################################
# Initial population size
# based on Boulanger et al. 2011 J Wildlife Man
###############################################################################################

# Distribution of herd among stages/sex 
bouDist1985 <- c(yearling = 88000, bull = 138000, calf = 176000, cow = 240000)

# Proportion of adults that are female for breeding purposes
propFemale <- bouDist1985['cow']/(bouDist1985['bull'] + bouDist1985['cow'])
