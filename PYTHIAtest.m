function out = PYTHIAtest(model, Z, Y, Ybin, Ybest, algolabels)
% PYTHIAtest  Deprecated. Use PYTHIA with a trained model (7-arg eval mode).
%
%   This function exists only for backward compatibility. All new code that
%   calls exploreIS will use PYTHIA eval mode directly.
warning('ISA:PYTHIAtest:deprecated', ...
    ['PYTHIAtest is deprecated and will be removed in a future release. ' ...
     'Use PYTHIA(Z, Y, Ybin, Ybest, algolabels, opts, trainedPythia) instead.']);
opts = struct('verbose', true);
out = PYTHIA(Z, Y, Ybin, Ybest, algolabels, opts, model);
end
