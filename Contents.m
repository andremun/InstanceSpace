% Instance Space Analysis Toolkit
% Version 2.0.0 (R2025a) 24-Jul-2026
%
% Standalone pipeline functions
%   PRELIM      - Data preprocessing (outlier bounding, normalisation)
%   SIFTED      - Feature selection (correlation + GA clustering)
%   PILOT       - 2D/3D projection (PBLDR algorithm)
%   CLOISTER    - Empirical bound estimation
%   PYTHIA      - Algorithm selection (classifier training/prediction)
%   TRACE       - Footprint construction (TRACE3 algorithm)
%   FILTER      - Density-based instance subsetting for small-scale experiments
%
% Class interface
%   InstanceSpace - Full pipeline wrapper class (build/explore/save/load)
%
% Backward-compatibility wrappers
%   buildIS     - Thin wrapper around InstanceSpace.build()
%   exploreIS   - Thin wrapper around InstanceSpace.explore()
%
% Utilities
%   ISAdefaults          - Fill in missing opts fields with default values
%   ISAvalidateOpts      - Validate user-supplied opts fields
%   ISAsubsetData         - Subset rows (and optionally feature columns) of a data struct
%   ISAgetClassifierFcn   - Resolve a classifier name from the registry
%   ISAmigrateModel       - Migrate legacy model.mat to current format
%
% Output
%   scriptcsv   - Write a model or explore() result to CSV files
%   scriptpng   - Write a model or explore() result to PNG figures
%   scriptweb   - Write colour-scaled CSV data for MATILDA's web tools
%   scriptfcn   - Shared plotting/CSV-writing helpers
%   ISArecallView - Snap a reopened 3D footprint figure to its stored viewpoint
%
% Example
%   example           - Getting-started example using the reference dataset
%   test_integration  - Exhaustive option-coverage regression suite
