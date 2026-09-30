classdef ExportReviewTest < matlab.unittest.TestCase
% ExportReviewTest  Regression cases for the September 2026 review.
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


    methods (Test)
        function testGeometryExports(tc)
            fixture = tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            root = [fixture.Folder filesep]; rng(2); Z = rand(30,3); Y = rand(30,2);
            opts = ISAdefaults(struct()); opts.pythia.skip = true;
            m.data = struct('X',Z,'Xraw',Z,'Y',Y,'Yraw',Y,'Ybin',true(30,2), ...
                'Ybest',min(Y,[],2),'P',ones(30,1),'beta',true(30,1),'numGoodAlgos',2*ones(30,1), ...
                'featlabels',{{'x','y','z'}},'algolabels',{{'a','b'}}, ...
                'instlabels',{cellstr("i"+(1:30)')});
            m.featsel.idx = 1:3; m.pilot.Z = Z;
            m.pythia = PYTHIA(Z,Y,m.data.Ybin,m.data.Ybest,m.data.algolabels,opts.pythia);
            m.trace = TRACE(Z,m.data.Ybin,[],m.data.P,m.data.beta,m.data.algolabels,opts.trace);
            m.cloist = CLOISTER(Z,eye(3),opts.cloister);
            scriptcsv(m,root);
            V = readmatrix([root 'footprint_a_good.csv']); V = V(:,2:end);
            F = readmatrix([root 'footprint_a_good_faces.csv']);
            tc.verifyTrue(all(F(:)>=1 & F(:)<=size(V,1)));
            meshVolume = abs(sum(dot(V(F(:,1),:),cross(V(F(:,2),:),V(F(:,3),:),2),2)))/6;
            tc.verifyEqual(meshVolume,volume(m.trace.good{1}.polygon),'AbsTol',1e-9);
            writematrix(42,[root 'user_results.csv']);
            m.data.algolabels{2} = 'c'; m.trace.good{1}.polygon = []; m.trace.best{1}.polygon = [];
            scriptcsv(m,root);
            tc.verifyFalse(isfile([root 'footprint_b_good.csv']));
            tc.verifyEqual(height(readtable([root 'footprint_a_good.csv'])),0);
            tc.verifyTrue(isfile([root 'user_results.csv']));
            m.pilot.Z = Z(:,1:2);
            m.trace = TRACE(m.pilot.Z,m.data.Ybin,[],m.data.P,m.data.beta,m.data.algolabels,opts.trace);
            m.cloist = CLOISTER(Z,[1 0 0;0 1 0],opts.cloister);
            scriptcsv(m,root);
            tc.verifyFalse(isfile([root 'footprint_a_good_faces.csv']));
            manifest = jsondecode(fileread([root 'geometry_manifest.json']));
            tc.verifyEqual(manifest.dimensions,2);
        end
    end
end
