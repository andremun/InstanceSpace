function [X, out] = SIFTED2(X, Y, Ybin, featlabels, opts)
% SIFTED2  Deprecated. Use SIFTED instead.
%
% This wrapper exists for backward compatibility. All new code should call
% SIFTED directly; SIFTED2 was promoted and renamed as part of the v1.7
% refactor (spec Phase 6).
warning('ISA:SIFTED2:deprecated', ...
    'SIFTED2 is deprecated and will be removed in a future release. Call SIFTED instead.');
[X, out] = SIFTED(X, Y, Ybin, featlabels, opts);
end
