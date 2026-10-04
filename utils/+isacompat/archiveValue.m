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
function out=archiveValue(value,decode)
% Versioned plain-data encoding. Only explicitly supported types reconstruct.
if decode && isstruct(value) && isscalar(value) && isfield(value,'isaArchiveType')
    switch value.isaArchiveType
        case 'knn'
            filename=[tempname '.mat']; cleanup=onCleanup(@() delete(filename));
            fid=fopen(filename,'wb'); fwrite(fid,value.payload,'uint8'); fclose(fid);
            out=loadmodel(filename);
        case 'alpha'
            out=isacompat.AlphaShape(value.points); out.Alpha=value.alpha;
            out.RegionThreshold=value.threshold;
        case 'partition', out=isacompat.Partition(value.masks);
        case 'categorical', out=categorical(value.values,value.categories,value.categories);
        case 'string', out=string(value.values);
        otherwise, error('ISA:compat:archiveType','Unknown archive type.');
    end
elseif ~decode && isa(value,'ClassificationKNN')
    filename=[tempname '.mat']; cleanup=onCleanup(@() delete(filename));
    savemodel(value,filename); fid=fopen(filename,'rb'); bytes=fread(fid,Inf,'*uint8'); fclose(fid);
    out=struct('isaArchiveType','knn','payload',bytes);
elseif ~decode && isa(value,'isacompat.AlphaShape')
    out=struct('isaArchiveType','alpha','points',value.Points,'alpha',value.Alpha,'threshold',value.RegionThreshold);
elseif ~decode && (isa(value,'cvpartition') || isa(value,'isacompat.Partition'))
    masks=[]; for k=1:value.NumTestSets, masks(:,k)=test(value,k); end
    out=struct('isaArchiveType','partition','masks',logical(masks));
elseif ~decode && isa(value,'categorical')
    out=struct('isaArchiveType','categorical','values',{cellstr(value)},'categories',{categories(value)});
elseif ~decode && isa(value,'string')
    out=struct('isaArchiveType','string','values',{cellstr(value)});
elseif isstruct(value)
    out=value; fields=fieldnames(value);
    for k=1:numel(value)
        for j=1:numel(fields), out(k).(fields{j})=isacompat.archiveValue(value(k).(fields{j}),decode); end
    end
elseif iscell(value)
    out=cell(size(value));
    for k=1:numel(value), out{k}=isacompat.archiveValue(value{k},decode); end
elseif isnumeric(value) || islogical(value) || ischar(value)
    out=value;
else
    error('ISA:compat:archiveType','Unsupported archive value of class %s.',class(value));
end
end
