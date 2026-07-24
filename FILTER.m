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
function [subsetIndex,isDissimilar,isVISA] = FILTER(X,Y,Ybin,opts)
% FILTER  Density-based instance subsetting for small-scale experiments.
%
%   [subsetIndex,isDissimilar,isVISA] = FILTER(X,Y,Ybin,opts)
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
%
%   Note: this also prints the feature-space uniformity of the retained
%   subset (of the instances NOT marked redundant) to the console.

[ninst,nalgos] = size(Y);
nfeats = size(X,2);

subsetIndex = false(ninst,1);
isDissimilar = true(ninst,1);
isVISA = false(ninst,1);
gamma = sqrt(nalgos/nfeats)*opts.mindistance;

for ii=1:ninst
    if ~subsetIndex(ii)
        for jj=ii+1:ninst
            if ~subsetIndex(jj)
                Dx = pdist2(X(ii,:),X(jj,:));
                Dy = pdist2(Y(ii,:),Y(jj,:));
                Db = all(Ybin(ii,:) & Ybin(jj,:));
                if  Dx <= opts.mindistance
                    isDissimilar(jj) = false;
                    switch opts.type
                        case 'Ftr'
                            subsetIndex(jj) = true;
                        case 'Ftr&AP'
                            if Dy <= gamma
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
                                if Dy <= gamma
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
end

% Assess the uniformity of the data
D = squareform(pdist(X(~subsetIndex,:)));
ninst = size(D,1);
D(eye(ninst,'logical')) = NaN;
nearest = min(D,[],2,'omitnan');
model.data.unif = 1-(std(nearest)./mean(nearest));
disp(['Uniformity of the instance subset: ' num2str(model.data.unif,4)]);

end