# xe-calibration-consortium

## Set up

In order to read MRD files you will need to download/clone the ismrmrd repository and add it to your MATLAB path: https://github.com/ismrmrd/ismrmrd

## Description

- calibration_production_v2302.m: script to extract 129Xe calibration parameters and visualize calibration scan data in real time.
- bonusSpectraCalibration_v2302.m: script for analyzing bonus calibration spectra, currently only accepts Siemens twix files
- MRD_data_inspector.m: script for inspecting data and header variables from MRD file
- demo_siem_file_sort_v2108.m: script for inspecting data and header variables from Siemens twix files
- NOTE: Remaining files are required functions and classes for calibration_production_v2302.m and bonusSpectraCalibration_v2302.m

## Usage

Run scripts in MATLAB, a pop-up will prompt you to select an input file.
