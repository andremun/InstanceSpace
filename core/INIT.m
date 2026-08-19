function [data, extra] = INIT(rootdir, opts, trainedModel)
% INIT  Load and filter instance-space metadata (build or explore time).
%
%   Training mode  (2 args): [data, extra] = INIT(rootdir, opts)
%     Reads rootdir/metadata.csv, applies opts.selvars.feats/.algos
%     column filtering, and drops instances/features with too many
%     missing values (opts.prelim.nanThreshold).
%
%   Evaluation mode (3 args): [data, extra] = INIT(rootdir, opts, trainedModel)
%     Reads rootdir/metadata_test.csv, applies the same opts.selvars
%     filtering, validates the feature set/order against trainedModel
%     (a previously-built InstanceSpace object's obj.model), and
%     reconciles algorithm columns against trainedModel.data.algolabels:
%     known algorithms line up by name at their trained column position,
%     unseen ones are appended as new columns (NaN for training-only
%     algorithms).
%
%   Before this function existed, this logic was two independent,
%   drifted implementations -- the version embedded in
%   InstanceSpace.runPrelim (build time) and a separate inline
%   reimplementation in InstanceSpace.evaluateTestSet (explore time).
%   INIT is a pure extraction of both into one shared function, dispatched
%   by nargin like PYTHIA/TRACE (#38): each mode's behaviour is preserved
%   exactly as it was, just no longer duplicated in two places that could
%   silently drift apart from each other.
%
%   Inputs
%     rootdir - directory containing metadata.csv (training mode) or
%               metadata_test.csv (evaluation mode), trailing slash
%     opts    - the full opts struct (obj.opts in training mode,
%               trainedModel.opts in evaluation mode)
%     trainedModel - (optional) a previously-built InstanceSpace object's
%               obj.model; presence selects evaluation mode
%
%   Outputs
%     data  - struct: instlabels, S (if a 'source' column is present), X,
%             Y, Xraw, Yraw, algolabels. Training mode also: featlabels
%             (evaluation mode's equivalent is finalised later by the
%             caller, after featsel.idx subsetting -- see extra.featlabelsAll)
%     extra - training mode: struct() (unused)
%             evaluation mode: struct with featlabelsAll (feature labels,
%             'feature_'-prefixed, after opts.selvars.feats filtering but
%             before featsel.idx subsetting), modelalgos, newalgos, nalgos,
%             ninst -- all needed by the caller's downstream
%             PRELIM/normalisation/featsel step, which must treat the
%             first modelalgos columns of data.Y differently from newly
%             appended ones

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

narginchk(2, 3);
isEvalMode = nargin == 3;
extra = struct();

if isEvalMode
    fprintf('[EXPLORE] Loading metadata_test.csv.\n');
    datafile = [rootdir 'metadata_test.csv'];
else
    fprintf('[BUILD] Loading metadata.csv.\n');
    datafile = [rootdir 'metadata.csv'];
end
Xbar = readtable(datafile);
varlabels = Xbar.Properties.VariableNames;
isname   = strcmpi(varlabels, 'instances');
isfeat   = strncmpi(varlabels, 'feature_', 8);
isalgo   = strncmpi(varlabels, 'algo_', 5);
issource = strcmpi(varlabels, 'source');

data = struct();
data.instlabels = Xbar{:,isname};
if isnumeric(data.instlabels)
    data.instlabels = num2cell(data.instlabels);
    data.instlabels = cellfun(@(x) num2str(x), data.instlabels, 'UniformOutput', false);
end
if any(issource)
    data.S = categorical(Xbar{:,issource});
end
data.X = Xbar{:,isfeat};
data.Y = Xbar{:,isalgo};

if isEvalMode
    % ---------------------------------------------------------------
    [ninst, nalgos] = size(data.Y);

    % Mirror the training-mode opts.selvars.algos restriction: without
    % this, data.Y/algolabels keep every raw algo_ column from
    % metadata_test.csv, so an algorithm deliberately excluded from
    % training via opts.selvars.algos looks like a "new" algorithm to
    % the reconciliation step below -- appended as an extra column with
    % real performance data but no trained classifier -- instead of
    % simply being excluded from evaluation the same way it was
    % excluded from training.
    algolabelsAll = varlabels(isalgo);
    if isfield(opts, 'selvars') && isfield(opts.selvars, 'algos')
        isselalgo = false(1, length(algolabelsAll));
        for i = 1:length(opts.selvars.algos)
            isselalgo = isselalgo | strcmp(algolabelsAll, opts.selvars.algos{i});
        end
        data.Y = data.Y(:,isselalgo);
        algolabelsAll = algolabelsAll(isselalgo);
        nalgos = size(data.Y, 2);
    end

    % Mirror the training-mode opts.selvars.feats restriction (applied at
    % build time, before PRELIM computed model.prelim.hibound/lobound/
    % minX/etc, all of which are sized to the RESTRICTED feature count):
    % without this, data.X keeps every raw feature_ column from
    % metadata_test.csv, and the bound-clipping below fails with a
    % dimension mismatch against those restricted-size arrays -- or, if
    % it happened not to error, featsel.idx would silently select the
    % wrong columns positionally instead of the intended named features.
    featlabelsAll = varlabels(isfeat);
    if isfield(opts, 'selvars') && isfield(opts.selvars, 'feats')
        isselfeat = false(1, length(featlabelsAll));
        for i = 1:length(opts.selvars.feats)
            isselfeat = isselfeat | strcmp(featlabelsAll, opts.selvars.feats{i});
        end
        data.X = data.X(:,isselfeat);
        featlabelsAll = featlabelsAll(isselfeat);
    end
    % Validate against trainedModel.prelim.lambdaX (one Box-Cox lambda
    % per feature PRELIM actually fit at training time, i.e. after both
    % opts.selvars.feats AND any opts.prelim.nanThreshold-triggered
    % column drops -- the latter isn't mirrored above, since which
    % columns those were isn't retained anywhere in the model). A count
    % mismatch here means metadata_test.csv doesn't match the training
    % feature set, and would otherwise surface many rows down as an
    % opaque bsxfun dimension-mismatch error.
    if numel(featlabelsAll) ~= numel(trainedModel.prelim.lambdaX)
        error('ISA:INIT:featureCountMismatch', ...
            ['metadata_test.csv has %d feature column(s) after applying ' ...
             'opts.selvars.feats, but the trained model expects %d ' ...
             '(some training features may have been dropped for exceeding ' ...
             'opts.prelim.nanThreshold). Check that metadata_test.csv has the ' ...
             'same feature_ columns as metadata.csv.'], ...
            numel(featlabelsAll), numel(trainedModel.prelim.lambdaX));
    end
    % A count match alone doesn't guarantee the same feature ORDER: every
    % bound-clipping/Box-Cox/normalisation step downstream indexes
    % trainedModel.prelim.hibound/lobound/lambdaX/etc positionally, so if
    % metadata_test.csv lists its feature_ columns in a different order
    % than metadata.csv did, columns would be silently misaligned --
    % normalised/clipped with the wrong feature's parameters -- instead
    % of erroring. Only checkable for models that retained the training
    % order (featsel.labels); older saved models without it just skip
    % this extra check.
    if isfield(trainedModel.featsel, 'labels')
        testLabels = strrep(featlabelsAll, 'feature_', '');
        if ~isequal(testLabels(:), trainedModel.featsel.labels(:))
            error('ISA:INIT:featureOrderMismatch', ...
                ['metadata_test.csv''s feature_ columns are not in the same order ' ...
                 'as the trained model''s (after applying opts.selvars.feats). ' ...
                 'Expected order: %s. Got: %s. Reorder metadata_test.csv''s columns ' ...
                 'to match metadata.csv.'], ...
                strjoin(trainedModel.featsel.labels, ', '), strjoin(testLabels, ', '));
        end
    end

    % Reconcile test algorithms against the trained model's: known
    % algorithms line up by name, unseen ones are appended as new
    % columns (NaN for training-only algorithms).
    data.algolabels = strrep(algolabelsAll, 'algo_', '');
    algoexist = zeros(1, nalgos);
    for ii = 1:nalgos
        aux = find(strcmpi(strtrim(data.algolabels{ii}), strtrim(trainedModel.data.algolabels)));
        if ~isempty(aux)
            algoexist(ii) = aux;
        end
    end
    newalgos  = sum(algoexist == 0);
    modelalgos = length(trainedModel.data.algolabels);
    Yaux  = NaN + ones(ninst, modelalgos+newalgos);
    lblaux = trainedModel.data.algolabels;
    acc = modelalgos + 1;
    for ii = 1:nalgos
        if algoexist(ii) == 0
            Yaux(:,acc) = data.Y(:,ii);
            lblaux(:,acc) = data.algolabels(ii);
            acc = acc + 1;
        else
            Yaux(:,algoexist(ii)) = data.Y(:,ii);
        end
    end
    data.Y = Yaux;
    data.algolabels = lblaux;
    nalgos = size(data.Y, 2);

    data.Xraw = data.X;
    data.Yraw = data.Y;

    extra.featlabelsAll = featlabelsAll;
    extra.modelalgos = modelalgos;
    extra.newalgos = newalgos;
    extra.nalgos = nalgos;
    extra.ninst = ninst;
else
    % ---------------------------------------------------------------
    data.featlabels = varlabels(isfeat);
    if isfield(opts, 'selvars') && isfield(opts.selvars, 'feats')
        msg = 'Using the following features: ';
        isselfeat = false(1, length(data.featlabels));
        for i = 1:length(opts.selvars.feats)
            isselfeat = isselfeat | strcmp(data.featlabels, opts.selvars.feats{i});
            msg = [msg opts.selvars.feats{i} ' ']; %#ok<AGROW>
        end
        fprintf('[BUILD] %s\n', msg);
        data.X = data.X(:,isselfeat);
        data.featlabels = data.featlabels(isselfeat);
    end

    data.algolabels = varlabels(isalgo);
    if isfield(opts, 'selvars') && isfield(opts.selvars, 'algos')
        msg = 'Using the following algorithms: ';
        isselalgo = false(1, length(data.algolabels));
        for i = 1:length(opts.selvars.algos)
            isselalgo = isselalgo | strcmp(data.algolabels, opts.selvars.algos{i});
            msg = [msg opts.selvars.algos{i} ' ']; %#ok<AGROW>
        end
        fprintf('[BUILD] %s\n', msg);
        data.Y = data.Y(:,isselalgo);
        data.algolabels = data.algolabels(isselalgo);
    end

    idx = all(isnan(data.X), 2) | all(isnan(data.Y), 2);
    if any(idx)
        warning('-> There are instances with too many missing values. They are being removed to increase speed.');
        data.X = data.X(~idx,:);
        data.Y = data.Y(~idx,:);
        data.instlabels = data.instlabels(~idx);
        if isfield(data, 'S')
            data.S = data.S(~idx);
        end
    end
    idx = mean(isnan(data.X), 1) >= opts.prelim.nanThreshold;
    if any(idx)
        warning('-> There are features with too many missing values. They are being removed to increase speed.');
        data.X = data.X(:,~idx);
        data.featlabels = data.featlabels(~idx);
    end
    ninst = size(data.X, 1);
    nuinst = size(unique(data.X, 'rows'), 1);
    if nuinst/ninst < 0.5
        warning('-> There are too many repeated instances. It is unlikely that this run will produce good results.');
    end

    data.Xraw = data.X;
    data.Yraw = data.Y;
    data.featlabels = strrep(data.featlabels, 'feature_', '');
    data.algolabels = strrep(data.algolabels, 'algo_', '');
end

end
