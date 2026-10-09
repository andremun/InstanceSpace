classdef BoundaryPolicyTest < matlab.unittest.TestCase
    % Shared analytically labelled MATLAB/Python boundary-distance contract.
    methods (Test)
        function sharedContract(tc)
            source=jsondecode(fileread(fullfile(fileparts(mfilename('fullpath')),...
                'contracts','trace_boundary_cases.json')));
            transforms=[1 0;16 1000;0.125 -16];
            for c=1:numel(source.cases)
                item=source.cases{c};
                for k=1:size(transforms,1)
                    scale=transforms(k,1); offset=transforms(k,2);
                    if strcmp(item.kind,'polygon')
                        shell=item.shell*scale+offset;
                        geometry=polyshape(shell);
                        hole=squeeze(item.holes(1,:,:))*scale+offset;
                        geometry=subtract(geometry,polyshape(hole));
                    else
                        geometry=alphaShape(item.vertices*scale+offset,Inf);
                    end
                    q=item.queries*scale+offset; tolerance=item.tolerance*scale;
                    tc.verifyEqual(ISAfootprintContains(geometry,q),logical(item.exact));
                    tc.verifyEqual(ISAfootprintContains(geometry,q,tolerance),logical(item.tolerant));
                    enlarged=ISAfootprintContains(geometry,[q;ones(1,size(q,2))*1e12],tolerance);
                    tc.verifyEqual(enlarged(1:end-1),logical(item.tolerant));
                    for row=1:size(q,1)
                        tc.verifyEqual(ISAfootprintContains(geometry,q(row,:),tolerance),logical(item.tolerant(row)));
                    end
                end
            end
        end
        function fittedToleranceSurvivesPersistence(tc)
            opts=ISAdefaults(struct()); tc.verifyEqual(opts.trace.boundaryTolerance,0);
            fp=struct('polygon',polyshape([0 0;1 0;1 1;0 1]),'measure',1,...
                'measureLabel','Area','elements',1,'goodElements',1,'density',1,'purity',1);
            trained=struct('good',{{fp}},'best',{{fp}},'hard',fp,'space',fp,...
                'summary',{cell(2,11)},'boundaryTolerance',0.001);
            filename=[tempname '.mat']; cleanup=onCleanup(@() delete(filename));
            save(filename,'trained'); restored=load(filename);
            q=[0.5 -0.0005;0.5 -0.002];
            evaluated=TRACE(q,true(2,1),[],ones(2,1),true(2,1),{'a'},opts.trace,restored.trained);
            tc.verifyEqual(evaluated.good{1}.elements,1);
            tc.verifyEqual(evaluated.good{1}.goodElements,1);
            tc.verifyEqual(evaluated.boundaryTolerance,0.001);
        end
        function invalidTolerance(tc)
            for value=[-1 Inf NaN]
                tc.verifyError(@() ISAfootprintContains(polyshape([0 0;1 0;0 1]),[0 0],value),...
                    'ISA:footprint:tolerance');
                tc.verifyError(@() ISAvalidateOpts(struct('trace',struct('boundaryTolerance',value))),...
                    'ISA:ISAvalidateOpts:notPositive');
            end
        end
    end
end
