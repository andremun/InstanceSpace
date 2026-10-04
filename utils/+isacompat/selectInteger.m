function [ind, info] = selectInteger(fcn, upper, settings, workers)
% selectInteger  Integer feature-search backend retaining MATLAB native GA.
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
%   Smith-Miles, K. & Munoz, M.A. (2023). Instance Space Analysis for
%   Algorithm Testing. ACM Computing Surveys, 55(12), Article 255.
%   https://doi.org/10.1145/3572895
% -------------------------------------------------------------------------

if isacompat.isOctave()
    [ind,info]=isacompat.integerSearch(fcn,upper,settings);
else
    options=optimoptions('ga','FitnessLimit',settings.FitnessLimit, ...
        'FunctionTolerance',settings.FunctionTolerance,'MaxGenerations',settings.MaxGenerations, ...
        'MaxStallGenerations',settings.MaxStallGenerations,'PopulationSize',settings.PopulationSize, ...
        'UseParallel',workers>0);
    d=numel(upper);
    ind=ga(fcn,d,[],[],[],[],ones(1,d),upper,[],1:d,options);
    info=struct('backend','matlab-ga');
end
end
