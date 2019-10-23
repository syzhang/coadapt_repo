# Cognitive control of brain-machine interfaces for pain.

## Introduction

This is the code used for the analyses in the paper 'Cognitive control of brain-machine interfaces for pain' (submitted, currently in review). An older preprint version of the paper can be found on biorxiv: [Endogenous controllability of closed-loop brain-machine interfaces for pain](https://www.biorxiv.org/content/10.1101/369736v1).

The pipeline takes the fMRI dicom files generated in the experiment described in the paper above as input. The data needs to be separated according to experimental days:
* day 1: decoder construction sessions,
* day 2: adaptive control sessions.

The code can only work when the data is organised following the structure described below. The pipeline will then be able to preprocess the dicom files, extract features from ROI masks, conduct cross-validation, and train SLR decoders used in the real-time experiment.

## Prerequisites
To run the code in this repo, you will need the following software:
* MATLAB
* SPM12
  
## Usage
Run the following functions in MATLAB in the following order:
```matlab
run_preprocessing_00(); % returns preprocessed .nii files and normalised ROI masks
run_make_feature_01(); % returns ROI features extracted
run_cv_multiple_02(); % returns cross-validation results
run_train_classifier_03(); % returns template and weight files required for real-time experiment
```
Please change data paths and parameter settings within the scripts.

## Folder structure
```
├───dat                             # data and output folders
│   └───01                          # subject 01
│       ├───1stEPI                      # 1st EPI (.dcm) reference scan
│       │   ├───d1                          # day 1 reference
│       │   └───d2                          # day 2 reference
│       ├───1stNii                      # 1st EPI (.nii) reference scan
│       │   ├───d1
│       │   └───d2
│       ├───d1                          # original EPI data (.dcm) from day 1
│       ├───d2                          # original EPI data (.dcm) from day 2
│       ├───dummy                       # storing dummy scans 
│       │   └───d1
│       │   └───d2
│       ├───EPI                         # EPI data and output in .nii 
│       │   ├───d1                          # original EPI data converted to .nii (day 1)
│       │   ├───d2                          # original EPI data converted to .nii (day 1)
│       │   ├───tmpl_d1_d1                  # realigned day 1 EPI with day 1 reference scan
│       │   ├───tmpl_d1_d2                  # realigned day 2 EPI with day 1 reference scan
│       │   ├───tmpl_d2_d1                  # realigned day 1 EPI with day 2 reference scan
│       │   └───tmpl_d2_d2                  # realigned day 2 EPI with day 2 reference scan
│       ├───params                      # parameters generated
│       │   ├───roi_signal_scantrend        # raw signals extracted from ROIs
│       │   │   └───tmpl_d1_d1
│       │   │       └───001                     # ROI number (can run multiple at once)
│       │   │   └───...                         # remaining ROI signals
│       │   ├───scan_idx_d1                 # day 1 labels
│       │   └───scan_idx_d2                 # day 2 labels
│       ├───ROI                         # processed ROI files
│       │   ├───d1                          # ROIs normalised to day 1 reference
│       │   └───d2                          # ROIs normalised to day 2 reference
│       └───t1                          # Structural t1 scans + normalisation output
│           └───1                           # original t1 data (.dcm)
│
│   └───02                          # subject 02
│
│   └───...                         # remaining subjects
│
│   └───MNI_ROI_to_use              # ROI masks in MNI space
├───res                             # result folder containing trained classifiers
│   ├───cv_results                  # cross-validation results
│   └───weight_txt                  # trained SLR classifier weights (.txt) + reference scan copy
│       └───01                          # subject 01
│       └───...                         # remaining subjects
├───script                          # scripts for project
├───SLR                             # trained SLR weights (.mat) as backup
│   └───01                               # subject 01
│   └───...                              # remaining subjects
└───SLR1.51                         # SLR toolbox used
    └───TEST
```
## Citation
Zhang S, Yoshida W, Mano H, Yanagisawa T, Shibata K, Kawato M, Seymour B. 2018. Endogenous controllability of closed-loop brain-machine interfaces for pain. bioRxiv 369736. doi:10.1101/369736
