# 2026-08-14
# CORRE Guillaume PhD

This folder contains material to reproduce the benchmark of GENETHOFF with GUIDE-seq (v1.0.2 / V2) and iGUIDE-seq (v1.2.0).

'benchmarking_steps.sh'                     	:   Explains how to collect the libraries from SRA and how to start each pipeline.
'PRJNA506241_experiments.csv'    				:   Contains library metadata
'Benchmarking_analysis_script.R' 				:   Script to aggregate cleavage site fro mthe 4 pipelines, make plots and calculate Kendall's W score 
'complete_OT_table_all_methods_10bp.csv'		:	Table with all cleavage sites from the 4 pipelines, clustered in a 10bp window.
'complete_OT_table_all_methods_wide_10bp.csv'	:	Same as above but in a 'wide' format.

Folders with pipelines name contain : 
	- the manifest files required to start each pipeline and/or sample informations.
	- the output from the pipeline containing positions of leavage sites (for iGUIDEseq, those were extract from the report html file).