% Instance Space Analysis Toolkit
% Version 0.9.1 (R2025a) 20-Aug-2026
%
% Standalone pipeline functions
%   PRELIM      - Data preprocessing (outlier bounding, normalisation)
%   SIFTED      - Feature selection (correlation + GA clustering)
%   PILOT       - 2D/3D projection (PBLDR algorithm)
%   CLOISTER    - Empirical bound estimation
%   PYTHIA      - Algorithm selection (classifier training/prediction)
%   TRACE       - Footprint construction (TRACE3 algorithm)
%   TRACE_legacy - Pre-refactor DBSCAN + alpha-shape footprint algorithm
%                  (opts.trace.method = 'legacy', 2D only)
%   FILTER      - Density-based instance subsetting for small-scale experiments
%   PILOTviewpoint - Optimal 2D camera viewpoint(s) of a 3D PILOT projection
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
%   liveDemoIS        - Interactive, stage-by-stage walkthrough (open in
%                        MATLAB's Live Editor)
