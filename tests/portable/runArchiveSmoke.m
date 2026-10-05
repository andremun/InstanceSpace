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
function runArchiveSmoke()
% Literal marker-like data must remain data, alongside real archived objects.
filename=[tempname '.mat']; cleanup=onCleanup(@() removeFile(filename));
% Native flattened MAT files can legitimately have archive-like custom fields.
for versionValue={1,999,'custom'}
    native=struct('archiveVersion',versionValue{1},'data',[1 NaN 3]);
    save(filename,'-struct','native','-v7');
    assert(isequaln(isacompat.loadModel(filename),native));
    native.payload=struct('userData',true);
    save(filename,'-struct','native','-v7');
    assert(isequaln(isacompat.loadModel(filename),native));
end
native=struct('archiveVersion','ordinary custom field');
save(filename,'-struct','native','-v7');
assert(isequaln(isacompat.loadModel(filename),native));
tags={'knn','alpha','partition','categorical','string','struct','unknown',[],42,{'nested'}};
for j=1:numel(tags)
    literal=struct(); literal.isaArchiveType=tags{j};
    literal.value=struct('isaArchiveType','alpha','points',[0 0;1 0;0 1]);
    literal.payload={struct('isaArchiveType','knn'),NaN};
    % Check top-level collisions as well as nested cells/struct arrays.
    roundTrip(filename,literal);
    model=struct(); model.literal=literal;
    model.cells={literal,{literal}};
    model.array=repmat(literal,[2 1 2]);
    model.empty=repmat(literal,[0 3 2]);
    model.fieldless=repmat(struct(),[2 0]);
    model.archiveVersion=999; model.payload=literal;
    roundTrip(filename,model);
end
% These are valid object-marker shapes too: the old decoder silently
% changed their types instead of merely rejecting an unknown tag.
roundTrip(filename,struct('isaArchiveType','partition','masks',logical([1 0;0 1])));
roundTrip(filename,struct('isaArchiveType','string','values',{{'plain data'}}));
roundTrip(filename,struct('isaArchiveType','alpha','points',[0 0;1 0;0 1], ...
    'alpha',1,'threshold',0));
% Portable object encoding is a shared contract even when MATLAB uses native files.
partition=isacompat.Partition(logical([1 0;0 1;1 0]));
restored=isacompat.archiveValue(isacompat.archiveValue(partition,false),true);
for k=1:2
    assert(isequal(test(restored,k),test(partition,k)));
    assert(isequal(training(restored,k),~test(partition,k)));
end
geometry=isacompat.AlphaShape([0 0;1 0;0 1]);
restored=isacompat.archiveValue(isacompat.archiveValue(geometry,false),true);
assert(area(restored)==.5 && isequal(inShape(restored,[.1 .1;2 2]),[true;false]));
for value={categorical({'a';'b';'a'}),string({'one','two'})}
    restored=isacompat.archiveValue(isacompat.archiveValue(value{1},false),true);
    assert(isequal(cellstr(restored),cellstr(value{1})));
end
expectArchiveError(@() isacompat.archiveValue(@sin,false),'ISA:compat:archiveType');
expectArchiveError(@() isacompat.archiveValue(struct('x',1),true),'ISA:compat:archiveType');
expectArchiveError(@() isacompat.archiveValue(struct('isaArchiveType','unknown'),true),'ISA:compat:archiveType');
expectArchiveError(@() isacompat.archiveValue(struct('isaArchiveType','struct','value',42),true),'ISA:compat:archiveType');
expectArchiveError(@() isacompat.archiveValue(1,false,1),'ISA:compat:archiveVersion');
expectArchiveError(@() isacompat.archiveValue(1,true,3),'ISA:compat:archiveVersion');
if isacompat.isOctave()
    X=[0 0;1 0;0 1;2 2;3 2;2 3]; Y=logical([0;0;0;1;1;1]);
    classifier=fitcknn(X,Y,'NumNeighbors',3,'Weights',[1;2;3;4;5;6]);
    geometry=isacompat.AlphaShape([0 0;1 0;0 1]);
    partition=isacompat.Partition(logical([1 0;0 1;1 0]));
    cats=categorical({'a';'b';'a'}); strings=string({'one','two'});
    original=struct('classifier',classifier,'geometry',geometry, ...
        'partition',partition,'cats',cats,'strings',strings,'literal',literal);
    isacompat.saveModel(filename,original);
    restored=isacompat.loadModel(filename);
    checkObjects(original,restored,X);
    assert(isequaln(restored.literal,literal));
    % Fixed schema-1 layout: do not generate the legacy struct envelope with
    % the schema-2 encoder, which would make the compatibility test circular.
    archiveVersion=1;
    payload=struct('metadata',struct('note','legacy data','values',[1 NaN 3]));
    fields={'classifier','geometry','partition','cats','strings'};
    for j=1:numel(fields)
        payload.(fields{j})=isacompat.archiveValue(original.(fields{j}),false);
    end
    save(filename,'archiveVersion','payload','-mat7-binary');
    legacy=isacompat.loadModel(filename);
    assert(isequaln(legacy.metadata,payload.metadata));
    checkObjects(original,legacy,X);
    % Re-saving a legacy archive upgrades it to the collision-safe schema.
    isacompat.saveModel(filename,legacy); disk=load(filename); assert(disk.archiveVersion==2);
    checkObjects(original,isacompat.loadModel(filename),X);
    archiveVersion=3; save(filename,'archiveVersion','payload','-mat7-binary');
    try
        isacompat.loadModel(filename);
        error('ISA:portable:missingError','Expected an unsupported schema error.');
    catch err
        assert(strcmp(err.identifier,'ISA:compat:archiveVersion'));
    end
end
% Exact portable envelopes still require Octave and a supported schema.
if ~isacompat.isOctave()
    for archiveVersion=[1 2 3]
        payload=struct('data',42);
        save(filename,'archiveVersion','payload','-v7');
        try
            isacompat.loadModel(filename);
            error('ISA:portable:missingError','Expected a runtime/schema error.');
        catch err
            assert(strcmp(err.identifier,'ISA:compat:archiveVersion'));
        end
    end
end
fprintf('[PORTABLE] PASS: native custom fields, literal markers, objects and schema-1 migration.\n');
end
function roundTrip(filename,original)
encoded=isacompat.archiveValue(original,false);
assert(isequaln(isacompat.archiveValue(encoded,true),original));
if isacompat.isOctave()
    isacompat.saveModel(filename,original); disk=load(filename);
    assert(disk.archiveVersion==2);
    assert(isequaln(isacompat.loadModel(filename),original));
end
end
function checkObjects(original,restored,X)
assert(isa(restored.classifier,'ClassificationKNN'));
[a,b]=predict(original.classifier,X); [c,d]=predict(restored.classifier,X);
assert(isequal(a,c) && max(abs(b(:)-d(:)))<1e-12);
assert(isa(restored.geometry,'isacompat.AlphaShape'));
assert(area(original.geometry)==area(restored.geometry));
assert(isequal(inShape(original.geometry,X),inShape(restored.geometry,X)));
for k=1:original.partition.NumTestSets
    assert(isequal(test(original.partition,k),test(restored.partition,k)));
    assert(isequal(training(original.partition,k),training(restored.partition,k)));
end
assert(isequal(cellstr(original.cats),cellstr(restored.cats)));
assert(isequal(cellstr(original.strings),cellstr(restored.strings)));
end
function removeFile(filename)
if isfile(filename), delete(filename); end
end

function expectArchiveError(f,id)
try, f(); catch err, assert(strcmp(err.identifier,id)); return; end
error('ISA:portable:missingError','Expected %s.',id);
end
