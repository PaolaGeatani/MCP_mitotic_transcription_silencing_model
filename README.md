# Modelling mitotic transcription silencing in early _Drosophila_ embryos
This page contains the code used to estimate the parameters governing mitotic transcriptional silencing from MCP transcription curves in blastoderm _Drosophila_ embryos in the paper "Distinct regulation across the transcription cycle drives mitotic transcriptional silencing".

## Overview
The code performs single-nucleus fitting of the MCP fluorescence signal to a mathematical equation composed of a sigmoidal decay (parameters kf, tm, tau) and an exponential decay term (ts and kappa).

## Contents
- `scripts/MCP_tracks_fitting_transcriptional_silencing_model.R`: per-nucleus model fitting, goodness-of-fit, correlation analysis, diagnostic plots

## Input format
The MCP tracks are stored in an excel file, one tab per nucleus. Each tab contains columns relative to time, subtracted fluorescence intensity, normalized fluorescence intensity (additional decay in the code). Each sheet is named as date_embryonumber_nucleusletter (e.g. Aug08_1_a).

## Usage
Set the input file with the data and the output directory where you want to store the results. Run the script: the fitted parameters will be saved in an excel file, one raw per nucleus.

## Requirements
R (version 4.3.2). Required packages are reported at the beginning of the code.

## Citation
If you want to use this code, please contact us and cite us (citation details will be added after publication).
