Disentangling influence over group speed and direction reveals multiple patterns of influence in moving meerkat groups
Readme file 

Raw data
 There is one RData file per meerkat group, located in data/movement/level0. Each file contains:
- allX: N * T matrix containing the x coordinates in UTM S34 of the individual trajectories for that group. N is the number of individuals in the group and T is the number of time steps at which this group was recorded (each column therefore represents one second of data)
- allY: N * T matrix containing the y coordinates in UTM S34 of the individual trajectories for that group.
- indInfo: data frame with N rows containing information about each individual, such as date of birth, social status, etc…
- movSummary: N * D data frame indicating whether each individual was recorded o r not on each day. D is the number of days for which data were collected for that group.
- timeLine: vector of length T indicating the date and time in UTC of each data point.
- dates:  vector of length D indicating all dates at which data were recorded for that group, in YYYYMMDD format.
- dayIdx: vector of length D+1 giving indexes for the location of the first data point of each day in the timeLine vector, plus the index of the very last data point.
Scripts
Here is a short description of each of the scripts. They should be executed in order, since each one relies on the output from the previous one. The working directory should be changed to the location of the project folder on the user’s machine.
• 1-Pre-processing: this script cleans the data. It is removing moments that would not reflect typical meerkat group movement (low number of individuals recorded, alarms, roving individuals, group encounters, inactivity…), as detailed in the main text.  Takes as input the RData coordinate files in data/movement/level0 as well as the scan data in data/scans, and outputs the cleaned coordinate files in data/movement/level1.

• 2-Computing_spatial_metrics: this script computes group and individual spatial metrics (mainly velocity vectors) from which the influence scores are going to be calculated, by spatially discretizing the trajectories. Takes as input the cleaned RData coordinate files in data/movement/level1 and outputs a table in output/ with each row being one second for one individual and columns being the different spatial metrics. 

• 3-Logistic_modelling: this script models individual turn and speed influence using a modified logistic regression function, then uses these model fits to get individual turn and speed influence scores. Takes as input the spatial metrics table, and outputs 4 tables (one for each influence type) with the individuals as rows and the fitted model parameters and influence scores as columns. 

• 4-Bootstraps: this script gets confidence intervals on the influence scores by bootstrapping the data. It first looks at the distribution of autocorrelation lag for the turn and speed influence across dates and individuals, then uses the mean of these distributions as the size of data chunks that are drawn during the bootstrap. It takes as input the spatial metric table and the logistic fits tables, and overwrites the latter with 2 additional columns for each table (lower and upper confidence interval of the influence score).

• 5-Figures: this script outputs the figures and statistical tests results found in the main text. It takes as input the spatial metric table, the logistic fits tables, and the cleaned RData coordinate files in data/movement/level.

• 6-Varying_discretization_step_length: this script checks for robustness of the results by re-running the analysis with different values for the spatial discretization threshold (5, 15, and 20 meters instead of 10 meters in the main text). It takes as input the cleaned RData coordinate files in data/movement/level1 and outputs spatial metrics tables, logistics fits tables and figures for each of the discretization thresholds. 

• 7-Supplemental_figures: this script outputs the figures found in the supplementary material. It takes as input the spatial metric table, the logistic fits tables and the cleaned RData coordinate files in data/movement/level1.
