classdef ReviewFixTest < matlab.unittest.TestCase
% ReviewFixTest  Regression cases for the September 2026 review.
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
        function testCVSelectionSummary(tc)
            rng(17); Z = rand(30,2); Y = rand(30,2); labels = Y<0.5;
            opts = ISAdefaults(struct()); opts.pythia.classifier = 'tree';
            opts.pythia.params = ones(2,1); opts.pythia.kFold = 3;
            out = PYTHIA(Z,Y,labels,min(Y,[],2),{'a','b'},opts.pythia);
            scores = out.Ysub .* max(out.precision',0);
            scores(isnan(scores)) = 0;
            [best, selection] = max(scores,[],2); selection(best<=0) = 0;
            tc.verifyEqual(out.selection0CV,selection);
            tc.verifyEqual(out.summary{1,8},'CV_model_precision');
            eval = PYTHIA(Z,Y,labels,min(Y,[],2),{'a','b'},opts.pythia,out);
            tc.verifyEqual(eval.summary{1,8},'Test_model_precision');
        end
        function testConstantScale(tc)
            opts = ISAdefaults(struct());
            p = opts.perf; p.AbsPerf = true; p.auto = true; p.bound = true; p.norm = true;
            X = [ones(20,1), (1:20)']; Y = [2*ones(20,1), (1:20)'];
            [a,b,t] = PRELIM(X,Y,p);
            [c,d] = PRELIM(X,Y,p,t);
            tc.verifyEqual(a,c,'AbsTol',1e-8); tc.verifyEqual(b,d,'AbsTol',1e-8);
            tc.verifyEqual(t.sigmaX(1),1); tc.verifyEqual(t.sigmaY(1),1);
            opts.pythia.skip = true;
            q = PYTHIA(X,Y,Y<5,min(Y,[],2),{'a','b'},opts.pythia);
            tc.verifyEqual(q.sigma(1),1);
        end
        function testPLSProjectionMean(tc)
            rng(12);
            X = rand(30,4)+[4 8 2 6]; Y = rand(30,3);
            opts = ISAdefaults(struct()); opts.pilot.method = 'pls';
            out = PILOT(X,Y,{'a','b','c','d'},opts.pilot);
            tc.verifyEqual((X-out.Xmean)*out.A', out.Z, 'AbsTol', 1e-10);
            b = CLOISTER(X,out.A,opts.cloister,out.Xmean);
            a = CLOISTER(X,out.A,opts.cloister);
            tc.verifyEqual(b.Zedge, a.Zedge-out.Xmean*out.A', 'AbsTol', 1e-10);
        end
        function testSelectionIgnoresTestLabels(tc)
            opts = ISAdefaults(struct());
            trained.classifiers = {struct('constant',true,'value',false), ...
                                   struct('constant',true,'value',false)};
            trained.precision = [1;1];
            trained.defaultAlgorithm = 2;
            Z = [0 1; 1 0];
            a = PYTHIA(Z, [1 9 1;1 9 1], logical([1 0 1;1 0 1]), [1;1], {'a','b','new'}, opts.pythia, trained);
            b = PYTHIA(Z, [9 1 1;9 1 1], logical([0 1 1;0 1 1]), [1;1], {'a','b','new'}, opts.pythia, trained);
            tc.verifyEqual(a.selection1, [2;2]);
            tc.verifyEqual(a.selection1, b.selection1);
            trained = rmfield(trained, {'precision','defaultAlgorithm'});
            c = PYTHIA(Z, [9 1;9 1], logical([0 1;0 1]), [1;1], {'a','b'}, opts.pythia, trained);
            tc.verifyEqual(c.selection1, [1;1]);
        end
    end
end
