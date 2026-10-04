function [best, info] = integerSearch(fitness, upper, options)
% integerSearch  Bounded integer GA, or exact enumeration for a small space.
% Uses elitism, tournament selection, uniform crossover and uniform mutation.
% The objective/folds are shared with MATLAB; stochastic trajectories differ.
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

upper=upper(:)'; d=numel(upper); count=prod(upper);
if count<=options.PopulationSize
    candidates=ones(count,d); ids=(0:count-1)'; stride=1;
    for j=1:d
        candidates(:,j)=1+mod(floor(ids/stride),upper(j)); stride=stride*upper(j);
    end
    scores=zeros(count,1);
    for i=1:count, scores(i)=fitness(candidates(i,:)); end
    [value,index]=min(scores); best=candidates(index,:);
    info=struct('backend','exact-integer-enumeration','fitness',value,'generations',1);
    return
end
population=1+floor(rand(options.PopulationSize,d).*upper);
stall=0; previous=Inf;
for generation=1:options.MaxGenerations
    scores=zeros(size(population,1),1);
    for i=1:numel(scores), scores(i)=fitness(population(i,:)); end
    [scores,order]=sort(scores); population=population(order,:);
    best=population(1,:); value=scores(1);
    if value<=options.FitnessLimit, break; end
    if previous-value<options.FunctionTolerance, stall=stall+1; else, stall=0; end
    if stall>=options.MaxStallGenerations, break; end
    previous=value;
    next=population;
    for i=2:size(population,1)
        ids=randi(size(population,1),2,2);
        [~,a]=min(scores(ids(:,1))); [~,b]=min(scores(ids(:,2)));
        child=population(ids(a,1),:); other=population(ids(b,2),:);
        cross=rand(1,d)<0.5; child(cross)=other(cross);
        mutate=rand(1,d)<1/d; replacement=1+floor(rand(1,d).*upper);
        child(mutate)=replacement(mutate); next(i,:)=child;
    end
    population=next;
end
info=struct('backend','serial-integer-ga','fitness',value,'generations',generation);
end
