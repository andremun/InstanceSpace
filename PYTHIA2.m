function out = PYTHIA2(Z, Y, Ybin, Ybest, algolabels, opts)
% PYTHIA2  Deprecated. Use PYTHIA with opts.classifier = 'knn'.
%
% This wrapper exists for backward compatibility. All new code should call
% PYTHIA directly with the desired classifier type set in opts.classifier.
warning('ISA:PYTHIA2:deprecated', ...
    ['PYTHIA2 is deprecated and will be removed in a future release. ' ...
     'Call PYTHIA with opts.classifier = ''knn'' instead.']);
opts.classifier = 'knn';
out = PYTHIA(Z, Y, Ybin, Ybest, algolabels, opts);
end
