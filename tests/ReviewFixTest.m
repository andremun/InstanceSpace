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
        function testUnobservedEvaluationRow(tc)
            p = struct('MaxPerf',false,'AbsPerf',true,'epsilon',5,'betaThreshold',.5,'auto',false);
            [~,~,out] = PRELIM([1;2],[NaN NaN;1 2],p);
            tc.verifyTrue(isnan(out.Ybest(1))); tc.verifyEqual(out.P(1),0);
            opts = ISAdefaults(struct());
            trained.classifiers = {struct('constant',true,'value',true),struct('constant',true,'value',true)};
            trained.precision = [1;1]; trained.defaultAlgorithm = 1;
            evaluated = PYTHIA([1;2],[NaN NaN;1 2],out.Ybin,out.Ybest,{'a','b'},opts.pythia,trained);
            tc.verifyEqual(evaluated.summary{end-1,4},1);
            tc.verifyEqual(evaluated.summary{end,4},1);
        end
        function testInvalidCVRejected(tc)
            opts = ISAdefaults(struct()); opts.pythia.classifier = 'svm';
            opts.pythia.params = [-1 1]; opts.pythia.kFold = 2; opts.pythia.verbose = false;
            Z = [(1:8)',[1;3;2;4;6;5;8;7]]; Y = (1:8)';
            tc.verifyError(@() PYTHIA(Z,Y,Y<5,Y,{'a'},opts.pythia),'ISA:PYTHIA:invalidCV');
        end
        function testScoreMetadata(tc)
            rng(14); Z = rand(20,2); Y = rand(20,1);
            opts = ISAdefaults(struct()); opts.pythia.classifier = 'svm';
            opts.pythia.params = [1 1]; opts.pythia.kFold = 2;
            out = PYTHIA(Z,Y,Y<0.5,Y,{'a'},opts.pythia);
            tc.verifyEqual(out.scoreTypeCV,{'probability'});
            tc.verifyTrue(all(out.Pr0subIsProbability));
            tc.verifyEqual(out.scoreType,{'probability'});
            tc.verifyTrue(all(out.Pr0hat>=0 & out.Pr0hat<=1));
        end
        function testTraceAcceptanceStatus(tc)
            Z = [0 0;1 0;0 1;1 1;.5 .5]; Z = [Z;Z];
            labels = [true(5,1);false(5,1)]; opts = ISAdefaults(struct());
            opts.trace.PI = .9; opts.trace.minAreaFrac = 0;
            out = TRACE(Z,labels,labels,ones(10,1),true(10,1),{'a'},opts.trace);
            tc.verifyFalse(out.good{1}.accepted);
            tc.verifyEqual(out.good{1}.terminationReason,'spectrumExhausted');
            tc.verifyLessThan(out.good{1}.purity,opts.trace.PI);
        end
        function testRegretWeights(tc)
            opts = ISAdefaults(struct()); opts.pythia.useweights = true;
            Y = [1 11;3 8;2 4;4 10]; Ybest = min(Y,[],2);
            out = PYTHIA([1 0;0 1;1 1;2 1],Y,true(4,2),Ybest,{'a','b'},opts.pythia);
            tc.verifyEqual(out.W,[2 10;2 5;2 2;2 6]);
        end
        function testSelectorRecall(tc)
            opts = ISAdefaults(struct());
            trained.classifiers = {struct('constant',true,'value',true),struct('constant',true,'value',true)};
            trained.precision = [1;1]; trained.defaultAlgorithm = 1;
            out = PYTHIA([0 1;1 0],ones(2),true(2),ones(2,1),{'a','b'},opts.pythia,trained);
            tc.verifyEqual(out.summary{end,9},100);
            trained.classifiers{1}.value = false; trained.classifiers{2}.value = false;
            out = PYTHIA([0 1;1 0],ones(2),true(2),ones(2,1),{'a','b'},opts.pythia,trained);
            tc.verifyEqual(out.summary{end,9},0);
        end
        function testTraceSmallEvaluation(tc)
            opts = ISAdefaults(struct()); rng(8);
            for dims = [2 3]
                Z = rand(40,dims); y = true(40,1);
                trained = TRACE(Z,y,y,ones(40,1),y,{'a'},opts.trace);
                for n = [1 2 4]
                    q = repmat(Z(1,:),n,1); labels = true(n,1);
                    evaluated = TRACE(q,labels,labels,ones(n,1),labels,{'a'},opts.trace,trained);
                    tc.verifyEqual(evaluated.space.measure,trained.space.measure);
                end
            end
        end
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
