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
%   | 'ensemble'| fitcensemble  | NumCycles [10,200]  | MinLeafSize [1,20]|
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
