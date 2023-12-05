DATA

"FAs_concentration.csv" is yellow perch and zooplankton fatty acid data  as concentrations

"FAs_percent.csv" is yellow perch and zooplankton fatty acid data as percent composition

"fish_pop.csv" is yellow perch population data from the mesocosms
	- YP.start is the number of yellow perch estimated to be present at the start of the experiment
	- YP.end is the number of yellow perch that were present at the end of the experiment
	- SS is the number of spot-tail shiners that found their way into the mesocosms and were recovered at the end of the experiment
	- BB is the number of burbot that found their way into the mesocosms and were recovered at the end of the experiment

"perch_biometrics.csv" has the yellow perch length and weight data
	- TL is total length and FL is fork length

"perch_diet.csv" has the results from the yellow perch stomach contents analysis

"zoop_biomass_2021.csv" is summary biomass data from the zooplankton community data from Desiree Langenfeld

"Perch Stocking Data.xlsx" is an Excel spreadsheet containing the mesocosm stocking and initial mortality data

CODE

Prefered order of operation for running the scripts, due to some inter-reliancies:

1. "Perch Survival and Growth.R"
2. "Perch Diet Analysis. R"
3. "Fatty Acids by Relative Composition.R"
4. "Fatty Acids by Concentration.R"