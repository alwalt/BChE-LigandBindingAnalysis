## Analysis of Small Ligands Binding to Butyrylcholinesterase

### Introduction:
This repository contains a collection of scripts and tools developed over the course of studying and analyzing the binding of small ligands to the protein Butyrylcholinesterase (BChE) via molecular dynamics simulations (MD).  These tools aid in various tasks such as data preparation, analysis, clustering, docking, visualization, and more.

### File Descriptions:

- **Kmeans_clustering_BChE.pl**: Perl script for K-means clustering specific to BChE protein system.
  
- **Multiple_ICM_Dockings.pl**: Automates multiple dockings using the ICM software.
  
- **Multiple_kmeans.py**: Python script to perform K-means clustering on multiple datasets.
  
- **am1_bcc.py**: Python script related to AM1-BCC charge calculation.
  
- **am1_bcc_img_prep.py**: Prepares images for visualization post AM1-BCC calculations.
  
- **atom2img.py**: Converts atomic data to image representation for visualization.
  
- **auto-plot-script.py**: Automated plotting tool for visualizing MD data.
  
- **backbone_vs_chains.pl**: Analyzes backbone vs. chains in the protein structure.
  
- **batch_prep_vectors.sh**: Batch script to prepare vector datasets.
  
- **batch_sum_kmeans.sh**: Batch script to summarize K-means results.
  
- **batch_summarize_kmeans.sh**: Advanced batch script to summarize K-means clustering results.
  
- **batch_trjconv.sh**: Script to automate trajectory conversions in batch.
  
- **contacts_analysis.pl**: Perl script to analyze contact points in protein structures.
  
- **dock_queue_1.pl**, **dock_queue_2.pl**, **dock_queue_3.pl**: Queue management scripts for docking simulations.
  
- **final_log_loop.sh**: Shell script to manage final logging operations.
  
- **final_str_analysis_with_key.pl**: Comprehensive structure analysis with reference to a key.
  
- **image.sh**: Shell script for visualization and image handling.
  
- **just_amino.pl**: Perl script to isolate amino acids in protein structure.
  
- **kmeans_analysis.1.pl** & **kmeans_analysis.pl**: Scripts for analyzing K-means clustering results.
  
- **kmeans_var.pl**: Analyzes variance in K-means clustering.
  
- **make_rms_rg_xvgs_index_mod.1.py** & **make_rms_rg_xvgs_index_mod.py**: Scripts for generating RMS, RG, and XVGS indices.
  
- **plot_runs.sh**: Shell script for plotting simulation runs.
  
- **prep.sh**: General preparation script.
  
- **prep_vectors_for_Kmeans.pl**: Prepares vector data for K-means clustering.
  
- **str_analysis_with_key.pl**: Protein structure analysis with reference to a key.
  
- **summarize_kmeans_results.pl**: Summarizes K-means clustering results.
  
- **vector_file_checkup.pl**: Checks vector data files for inconsistencies or errors.
  
- **xtc_to_pdb.pl**: Converts XTC trajectory files to PDB format.

### Usage:
Each script may have its own usage syntax. Refer to the comments or documentation within individual scripts for specific usage instructions.

### Requirements:
Ensure you have the necessary dependencies, software, and libraries installed. Some scripts might depend on external software or libraries for MD simulations, docking, or data visualization.

### License:
Please see the accompanying `LICENSE` file for licensing information.

### Contributions:
Contributions and feedback are welcome. Please raise issues on GitHub or submit pull requests if you have improvements or bug fixes.

### Note:
This project was last updated 6 years ago. Some tools might be outdated or require adjustments to work with newer versions of dependencies or software.
