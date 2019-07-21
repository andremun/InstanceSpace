function [out] = testIS(rootdir)
% -------------------------------------------------------------------------
% testIS.m
% -------------------------------------------------------------------------
%
% By: Mario Andres Munoz Acosta
%     School of Mathematics and Statistics
%     The University of Melbourne
%     Australia
%     2019
%
% -------------------------------------------------------------------------

startProcess = tic;
% -------------------------------------------------------------------------
% Collect all the data from the files
modelfile = [rootdir 'results.mat'];
datafile = [rootdir 'metadata.csv'];
optsfile = [rootdir 'options.json'];
if ~isfile(modelfile) || ~isfile(optsfile)
    error('Please place the datafiles in the directory specified.');
end
opts = jsondecode(fileread(optsfile));
disp('-------------------------------------------------------------------------');
disp('Listing options used:');
optfields = fieldnames(opts);
for i = 1:length(optfields)
    disp(optfields{i});
    disp(opts.(optfields{i}));
end

disp('-------------------------------------------------------------------------');
model = load(modelfile);
Xbar = readtable(datafile);
varlabels = Xbar.Properties.VariableNames;
isname = strcmpi(varlabels,'instances');
isfeat = strncmpi(varlabels,'feature_',8);
isalgo = strncmpi(varlabels,'algo_',5);
issource = strcmpi(varlabels,'source');
instlabels = Xbar{:,isname};
if isnumeric(instlabels)
    instlabels = num2cell(instlabels);
    instlabels = cellfun(@(x) num2str(x),instlabels,'UniformOutput',false);
end
if any(issource)
    S = categorical(Xbar{:,issource});
    sourcelabels = cellstr(unique(S));
end
X = Xbar{:,isfeat};
Y = Xbar{:,isalgo};

if model.opts.bound.flag
    himask = bsxfun(@gt,X,model.bound.hibound);
    lomask = bsxfun(@lt,X,model.bound.lobound);
    X = X.*~(himask | lomask) + bsxfun(@times,himask,model.bound.hibound) + ...
                                bsxfun(@times,lomask,model.bound.lobound);
end

if model.opts.norm.flag
    X = bsxfun(@minus,X,minX)+1;
    X = bsxfun(@rdivide,bsxfun(@power,X,model.norm.lambdaX)-1,model.norm.lambdaX);
    X = bsxfun(@rdivide,bsxfun(@minus,X,model.norm.muX),model.norm.sigmaX);
    
    Y(Y==0) = eps; % Assumes that Y is always positive and higher than 1e-16
    Y = bsxfun(@rdivide,bsxfun(@power,Y,model.norm.lambdaY)-1,model.norm.lambdaY);
    Y = bsxfun(@rdivide,bsxfun(@minus,Y,model.norm.muY),model.norm.sigmaY);
end

featlabels = 
X = X(:,model.featsel.idx);
Z = X*model.pbldr.A';

% -------------------------------------------------------------------------
% Removing the template data such that it can be used in the labels of
% graphs and figures.
featlabels = strrep(featlabels,'feature_','');
algolabels = strrep(algolabels,'algo_','');

end