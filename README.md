This repository contains the scripts and data for the following paper:

Jesser et al., Hemispheric lateralization determines poststroke cortical reorganization in the human brain, bioRxiv. DOI: (*will be added after uploading*)

Folders:
- data
  - controls
  - patients
  *.mat files with arrays of size 116 ROIs x 100 volumes per subject
- preprocessing
  - SPM12 preprocessing and job scripts for the healthy controls (Script\_control\*.m) and the stroke patients (Script\_mpi\*.m for scans acquired at the Max Planck Institute, and Script_ukt*.m for scans acquired at the Universitaetsklinikum Tuebingen)
 - analysis
   - Custom MATLAB scripts for the calculation of the graph metrics
