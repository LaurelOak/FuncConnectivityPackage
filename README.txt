README.TXT

This software and data package contains Matlab code routines needed to develop functional connectivity networks from water quality sample data. 

Code was developed by Laurel Larsen (laurel@berkeley.edu) for Larsen, Newman, Saunders, and Harvey, in review, "Complex networks of functional connectivity in an isolated wetland reconnected to its floodplain."

Users intending to use this code to develop similar functional connectivity networks with other data should first generate a spreadsheet in which each row corresponds to a different sample. Data columns should include sample date and all sample analytical information (i.e., analyte concentrations). The following pre-processing steps should then be undertaken:

1. Generate new columns representing sample values at previous time steps for each site and each analyte.
2. Use these columns to generate new columns representing the change in concentration at each site.
3. Use statistical software of your choice (we used JMP Statistical Software, SAS, Inc., Cary, NC) to develop an ANOVA model of the data that includes a site effect and a date effect (equivalent to a repeated-measures ANOVA). Calculate the residuals from the ANOVA model and add the residuals to the spreadsheet as a separate column. 
4. Repeat the ANOVA modeling process on the change in concentration, and generate a separate column with with residuals from that analysis.
5. For each sample date, copy the values of the analyses from all other sites as separate columns. (In other words, you are generating columns for every variable that might form a link to the of interest.)

The spreadsheet that you have generated above is the input to the set of Matlab codes that we include in this package. (Our exact input is the series of 'InputSS___.csv.' spreadsheets. There are four inputs, corresponding to the four different subsets of data for which functional connectivity networks were generated.)

The included series of codes conducts the set of four stepwise multiple linear regressions listed in Fig. 1 of the paper. Different codes were developed for different analytes, since each analyte's data is located in different columns of the input spreadsheet. Users of this code will need to adapt the introductory section of the code to ensure that the data are pulled from the correct set of columns. The FindCorrelationsTP.m code is fully commented to aid users; the codes for the other analytes are structured identically but are not fully commented. Note that the code requires Matlab's Statistics Toolbox to run. All code is supplied under the 2-clause BSD license, an approved open-source license. 

Results from the automated stepwise multiple linear regressions are saved as .xlsx files, in which the selected model is highlighted. We include these spreadsheets as an example of how other researchers may structure their own notes recording the outcomes of the program. Adjacency matrices in which the weights represent the correlation between nodes in the final functional connectivity networks are also included in this package, as .mat files


Package contents:
1. INPUT DATA
InputSSAllData.csv: Input spreadsheet for the analysis, which includes original data and the results of the ANOVA modeling (i.e., residuals), as described above. This file is the complete dataset. Analyte concentrations are given in mg/L.

InputSSNoFire.csv: Input spreadsheet for the analysis, not including the year of data (2011) impacted by the fire. Analyte concentrations are given in mg/L.

InputSSPostGap.csv: Input spreadsheet for the analysis, restricted to just the period of time after the levee gap was opened up. Analyte concentrations are given in mg/L.

InputSSNoFireNoFlood.csv: Input spreadsheet for the analysis for the pre-gap period of time, excluding the effects of the fire (2011). Analyte concentrations are given in mg/L.

2. ANALYSIS CODE

NetworkCalculations_TP.m: The only fully commented code supplied in this package. The other codes are included because they draw input data from different columns of the input spreadsheet, corresponding to the other analytes, but otherwise, all other codes are identical. 

NetworkCalculations_DOC.m

NetworkCalculations_logSO4.m

NetworkCalculations_Cl.m

NetworkCalculations_Ca.m

NetworkCalculations_logNH3.m

3. RESULTS
3.1. MODEL SELECTION
All of the spreadsheets contained within this folder house the results of the stepwise multiple linear regression models, with the selected models highlighted. There are four tabs on the spreadsheet, corresponding to the four periods of time analyzed. Each water quality analyte has its own spreadsheet, listed in the title. Results of the forwards and backwards regressions on the "diff" and "nondiff" data for each site are included within each tab, together with the significant variables included in each model, their coefficients, and their p-values. These spreadsheets are in Microsoft Excel format.

3.2. ADJACENCY MATRICES
Adjacency matrices for the final selected model, for each analyte, for each time period are included in this directory. File names correspond to the the analyte and time period. The matrices are symmetric, with rows/columns corresponding to sites. Sites are listed in the following (alphabetical) order:
C1
C2
DB1
DB2
DB3
L67C Canal
RS1
RS2
S1
UB1
UB2
UB3

