% -------------------------------------------------------------------------
% Instance Space Analysis (ISA) Toolkit
% Copyright (c) 2026 Mario Andres Munoz Acosta and contributors
% School of Computing and Information Systems
% The University of Melbourne, Australia
%
% SPDX-License-Identifier: LicenseRef-PolyForm-Noncommercial-1.0.0
% License: https://polyformproject.org/licenses/noncommercial/1.0.0/
%
% You may use, modify, and redistribute this software for non-commercial
% research and educational purposes only. Commercial use requires prior
% written permission. See the LICENSE file for full terms.
%
% Reference:
%   Simpson, C., Munoz, M.A., Kandanaarachchi, S. & Campello, R.J.G.B.
%   (2025). ISA3: A 3-dimensional expansion of Instance Space Analysis.
%   Machine Learning, 114, 240. https://doi.org/10.1007/s10994-025-06871-5
%
%   Smith-Miles, K. & Munoz, M.A. (2023). Instance Space Analysis for
%   Algorithm Testing. ACM Computing Surveys, 55(12), Article 255.
%   https://doi.org/10.1145/3572895
% -------------------------------------------------------------------------
function [fitFcn, p1label, p2label] = ISAgetClassifierFcn(name)
% ISAgetClassifierFcn  Resolve a classifier name from the registry.
%
%   [fitFcn, p1label, p2label] = ISAgetClassifierFcn(name)
%
%   name     : string matching an entry in the classifier registry
%   fitFcn   : MATLAB fitc* function handle
%   p1label  : human-readable label for the first Sobol-tuned hyperparameter
%   p2label  : human-readable label for the second Sobol-tuned hyperparameter
%              (may be 'N/A' for classifiers that have only one tunable parameter)
%
%   Registry:
%   +-----------+---------------+---------------------+------------------+
%   | name      | MATLAB fn     | param 1             | param 2          |
%   +-----------+---------------+---------------------+------------------+
%   | 'knn'     | fitcknn       | NumNeighbors [1,25] | Distance (cat.)  |
%   | 'svm'     | fitcsvm       | BoxConstraint log2  | KernelScale log2 |
%   | 'tree'    | fitctree      | MinLeafSize [1,100] | N/A              |
%   | 'nb'      | fitcnb        | Bandwidth  log10    | N/A              |
%   | 'linear'  | fitclinear    | Lambda     log10    | N/A              |
%   | 'ensemble'| fitcensemble  | NumLearningCycles [10,200] | MinLeafSize [1,20]|
%   +-----------+---------------+---------------------+------------------+
%
%   fitcecoc is excluded: PYTHIA trains one binary classifier per algorithm;
%   multi-class ECOC machinery is never required.

switch lower(name)
    case 'knn'
        fitFcn  = @fitcknn;
        p1label = 'NumNeighbors';
        p2label = 'Distance';
    case 'svm'
        fitFcn  = @fitcsvm;
        p1label = 'BoxConstraint';
        p2label = 'KernelScale';
    case 'tree'
        fitFcn  = @fitctree;
        p1label = 'MinLeafSize';
        p2label = 'N/A';
    case 'nb'
        fitFcn  = @fitcnb;
        p1label = 'Bandwidth';
        p2label = 'N/A';
    case 'linear'
        fitFcn  = @fitclinear;
        p1label = 'Lambda';
        p2label = 'N/A';
    case 'ensemble'
        fitFcn  = @fitcensemble;
        p1label = 'NumLearningCycles';
        p2label = 'MinLeafSize';
    otherwise
        error('ISA:ISAgetClassifierFcn:unknownClassifier', ...
            ['Unknown classifier ''%s''. ' ...
             'Valid options: knn, svm, tree, nb, linear, ensemble.'], name);
end
end
