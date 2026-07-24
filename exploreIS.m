function out = exploreIS(rootdir)
% exploreIS  Thin backward-compatibility wrapper around InstanceSpace (spec §7.5).
%
% Requires model.mat to already exist in rootdir (written by buildIS).
% Preserves the pre-Phase-7 exploreIS(rootdir) calling convention for
% callers -- notably the MATILDA web platform -- that invoke this entry
% point directly. New code should use InstanceSpace directly:
%
%   obj = InstanceSpace.load(rootdir);
%   obj = obj.explore(rootdir);
%   out  = obj.getResults(1);
%
% -------------------------------------------------------------------------
%
% By: Mario Andres Munoz Acosta
%     School of Mathematics and Statistics
%     The University of Melbourne
%     Australia
%     2019
%
% -------------------------------------------------------------------------

if ~(endsWith(rootdir, '/') || endsWith(rootdir, '\'))
    rootdir = [rootdir '/'];
end
if ~isfile([rootdir 'model.mat']) || ~isfile([rootdir 'metadata_test.csv'])
    error(['Please place the datafiles in the directory ''' rootdir '''']);
end

obj = InstanceSpace.load(rootdir);
obj = obj.explore(rootdir);
out = obj.testResults{end};

fprintf('EOF:SUCCESS\n');
end
