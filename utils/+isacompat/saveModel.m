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
function saveModel(filename,model)
if isacompat.isOctave()
    archiveVersion=2; payload=isacompat.archiveValue(model,false,archiveVersion);
    temporary=[tempname(fileparts(filename)) '.mat'];
    cleanup=onCleanup(@() removeTemporary(temporary));
    save(temporary,'archiveVersion','payload','-mat7-binary');
    [ok,msg]=movefile(temporary,filename,'f');
    if ~ok, error('ISA:compat:archiveWrite','%s',msg); end
else
    save(filename,'-struct','model','-v7.3');
end
end
function removeTemporary(filename)
if isfile(filename), delete(filename); end
end
