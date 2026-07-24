function model = buildIS(rootdir)
% buildIS  Thin backward-compatibility wrapper around InstanceSpace (spec §7.5).
%
% Preserves the pre-Phase-7 buildIS(rootdir) calling convention (a plain
% function taking a directory and returning the in-memory model struct)
% for callers -- notably the MATILDA web platform -- that invoke this
% entry point directly. New code should use InstanceSpace directly:
%
%   obj = InstanceSpace(rootdir);
%   obj = obj.build();
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
if ~isfile([rootdir 'metadata.csv']) || ~isfile([rootdir 'options.json'])
    error(['Please place the datafiles in the directory ''' rootdir '''']);
end

obj = InstanceSpace(rootdir);
obj = obj.build();
model = obj.model;

fprintf('EOF:SUCCESS\n');
end
