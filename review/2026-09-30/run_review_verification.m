function run_review_verification(group)
% Run all repository tests in four disjoint groups in clean MATLAB sessions.
% Group 1 runs non-option tests. Groups 2-4 partition independent option cases.
% Reports go to tempdir, one JUnit, coverage, and JSON result per group.
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(root,fullfile(root,'tests'),fullfile(root,'core'),fullfile(root,'utils'),fullfile(root,'output'));
set(groot,'defaultFigureVisible','off');
import matlab.unittest.TestSuite
import matlab.unittest.TestRunner
import matlab.unittest.plugins.XMLPlugin
import matlab.unittest.plugins.CodeCoveragePlugin
import matlab.unittest.plugins.codecoverage.CoberturaFormat
suite = TestSuite.fromFolder(fullfile(root,'tests'),'IncludingSubfolders',true);
option = startsWith({suite.Name},'PipelineOptionsTest/');
if group == 1
    suite = suite(~option);
else
    suite = suite(option);
    suite = suite(group-1:3:end);
end
runner = TestRunner.withTextOutput();
prefix = fullfile(tempdir,sprintf('isa-review-group-%d',group));
runner.addPlugin(XMLPlugin.producingJUnitFormat([prefix '.xml']));
files = {};
for folder = {'core','output','utils'}
    listing = dir(fullfile(root,folder{1},'*.m'));
    files = [files, fullfile({listing.folder},{listing.name})]; %#ok<AGROW>
end
files = [files,fullfile(root,{'InstanceSpace.m','buildIS.m','exploreIS.m'})];
runner.addPlugin(CodeCoveragePlugin.forFile(files,'Producing',CoberturaFormat([prefix '-coverage.xml'])));
results = runner.run(suite);
report = struct('names',{{results.Name}},'passed',[results.Passed], ...
    'failed',[results.Failed],'incomplete',[results.Incomplete],'duration',[results.Duration]);
fid = fopen([prefix '.json'],'w');
cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s',jsonencode(report));
fprintf('REVIEW_GROUP_%d: %d/%d passed\n',group,sum([results.Passed]),numel(results));
assertSuccess(results);
end
