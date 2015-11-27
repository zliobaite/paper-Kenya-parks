# Computational experiments for the Kenya national parks manuscript

This repository provides source code for reproducing the computational experiments reported in the manuscript **Modeling local climate from dental traits of mammals in 13 Kenyan national parks**, currently under review. 

Scripts are numbered in the order of appearance in the manuscript text. 

	run1_prepare_data.R
	
Extracts dental features from occurence matrix and trait data. Produces a dataset describing sites. Produces Table 1 (located in the `results` directory).

	run2_correlation_analysis.R
	
Computes correlations and produces Figure 1 and Figure A2. 

	run3_regression_analysis.R	
	
Produces and tests LARS and PLS regression models. Plots Figure 2 and Figure A1 (R2 accuracies). Produces decision formulas for Section 4.3.
	
	run4_residual_analysis.R
	
Produces residual analysis, plots Figure 3. 

	run5_residual_table.R
	
Produces Figure 5 and residuals for Tables 2 and 3.