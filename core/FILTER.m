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
function [subsetIndex,isDissimilar,isVISA,unif] = FILTER(X,Y,Ybin,opts)
% FILTER  Density-based instance subsetting for small-scale experiments.
%
%   [subsetIndex,isDissimilar,isVISA,unif] = FILTER(X,Y,Ybin,opts)
%
%   For every pair of instances closer than opts.mindistance in feature
%   space, the second is marked redundant according to opts.type (see
%   below); the caller (buildIS/InstanceSpace) typically keeps the
%   complement of subsetIndex, i.e. drops the redundant instances.
%
%   Inputs
%     X, Y, Ybin - (ninst x nfeats)/(ninst x nalgos)/(ninst x nalgos)
%                  feature, performance, and good-performance matrices
%     opts       - struct with fields:
%                    mindistance double  feature-space distance threshold
%                                        below which two instances are
%                                        considered too close
%                    type        char    extra condition (on top of feature
%                                        closeness) required before an
%                                        instance is marked redundant:
%                                        'Ftr' (none), 'Ftr&AP' (similar
%                                        algorithm performance too),
%                                        'Ftr&Good' (both instances good
%                                        on every algorithm), or
%                                        'Ftr&AP&Good' (both)
%
%   Outputs
%     subsetIndex  - (ninst x 1) logical; true where the instance was
%                    found redundant against an earlier, kept instance
%     isDissimilar - (ninst x 1) logical; false where an instance triggered
%                    a feature-space-closeness check against another
%     isVISA       - (ninst x 1) logical; true where two instances were
%                    feature-close but did not meet opts.type's extra
%                    condition, so neither was marked redundant despite
%                    the proximity ("visually important, but not
%                    subsetted away")
%     unif         - feature-space uniformity of the retained (non-
%                    redundant) subset: 1 minus the coefficient of
%                    variation of nearest-neighbour distances (closer to 1
%                    = more evenly spread, closer to 0 = clustered).
%                    Previously computed but assigned to an undefined
%                    workspace variable (model.data.unif) and discarded;
%                    now a real output.

[ninst,nalgos] = size(Y);
nfeats = size(X,2);

subsetIndex = false(ninst,1);
isDissimilar = true(ninst,1);
isVISA = false(ninst,1);
gamma = sqrt(nalgos/nfeats)*opts.mindistance;

% Precompute every pairwise feature/performance distance once instead of
% calling pdist2 per pair inside the O(ninst^2) loop below (that call
% overhead, repeated up to ninst^2 times, was FILTER's dominant cost).
% O(ninst^2) memory; fine up to a few thousand instances -- a KD-tree
% (knnsearch) would be needed to scale further.
Dx = squareform(pdist(X));
Dy = squareform(pdist(Y));
% Db(i,j) = all(Ybin(i,:) & Ybin(j,:)), which -- since all(A&B) is true
% iff all(A) and all(B) are both true -- reduces to "both i and j are
% good on every algorithm", a per-instance property rather than a real
% pairwise one. isGood(ii) && isGood(jj) below is that, computed once.
isGood = all(Ybin, 2);

% The elimination itself stays a sequential double loop, not just the
% distance lookups: which instances end up marked redundant depends on
% the running state of subsetIndex (an instance already marked redundant
% is skipped as both a future ii and jj), so this greedy process isn't
% safe to vectorise away without changing which instances get kept.
for ii=1:ninst
    if ~subsetIndex(ii)
        for jj=ii+1:ninst
            if ~subsetIndex(jj) && Dx(ii,jj) <= opts.mindistance
                isDissimilar(jj) = false;
                Db = isGood(ii) && isGood(jj);
                switch opts.type
                    case 'Ftr'
                        subsetIndex(jj) = true;
                    case 'Ftr&AP'
                        if Dy(ii,jj) <= gamma
                            subsetIndex(jj) = true;
                            isVISA(jj) = false;
                        else
                            isVISA(jj) = true;
                        end
                    case 'Ftr&Good'
                        if Db
                            subsetIndex(jj) = true;
                            isVISA(jj) = false;
                        else
                            isVISA(jj) = true;
                        end
                    case 'Ftr&AP&Good'
                        if Db
                            if Dy(ii,jj) <= gamma
                                subsetIndex(jj) = true;
                                isVISA(jj) = false;
                            else
                                isVISA(jj) = true;
                            end
                        else
                            isVISA(jj) = true;
                        end
                    otherwise
                        disp('Invalid flag!')
                end
            end
        end
    end
end

% Assess the uniformity of the retained (non-redundant) subset, reusing
% the already-computed Dx submatrix instead of a second pdist call.
Dkept = Dx(~subsetIndex, ~subsetIndex);
Dkept(1:size(Dkept,1)+1:end) = NaN; % diagonal -> NaN, excludes self-distance from min
nearest = min(Dkept,[],2,'omitnan');
unif = 1-(std(nearest)./mean(nearest));
fprintf('[FILTER] Uniformity of the instance subset: %.4g\n', unif);

end