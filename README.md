# Precise and Quantitative Chlorosis Severity Assessment Framework (PQCSAF) using Evolutionary Superpixels
# Developed by Sourav Samanta, Sanjoy Pratihar, and Sanjay Chatterji
## Description

This repository contains MATLAB code for implementing the core modules of the **PQCSAF (Precise and Quantitative Chlorosis Severity Assessment Framework)**, a system designed to assess chlorosis severity in plant leaves. The approach utilizes an evolutionary-optimized superpixel segmentation method based on **Cuckoo Search-optimized SLIC** and performs classification of superpixels to generate severity scores.

> ⚠️ **Note**:  This repository includes code for segmentation, feature extraction, classification, and severity estimation. A separate function is provided for feature selection, along with the accompanying dataset.



## Files Included

- `main_script.m` – Entry point for running the core framework.
- `cs_slic.m` – Custom function implementing Cuckoo Search-optimized SLIC for superpixel segmentation.
- `leaf_superpixel.m` – Function to extract individual leaf regions from the segmented image.
- `feature_extraction.m` – Function to compute Color GLCM texture features from each superpixel.
- `classification_superpixel.m` – Function to classify each superpixel based on extracted features.
- `severity_estimation.m` – Function to compute the chlorosis severity score from classification results.
- `display_score.m` – Function to visualize the severity scores of the leaf image.

- `The multics_feature_selection() function is designed for feature selection 
- `using a multi-swarm cuckoo search optimization approach.

---

## Requirements

- MATLAB R2016a or later
- MATLAB implementation of:
  - Simple Linear Iterative Clustering (SLIC)
  - Cuckoo Search (CS) algorithm
  - Color GLCM feature extraction

---

## Disclaimer

This code is intended for academic and research use. For inquiries regarding the complete pipeline or potential collaborations, please contact the authors.

