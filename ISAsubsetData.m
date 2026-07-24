function data = ISAsubsetData(data, subsetIndex, featIdx)
% ISAsubsetData  Subset rows (and optionally feature columns) of a data struct.
%   data = ISAsubsetData(data, subsetIndex) subsets all row-indexed fields.
%   data = ISAsubsetData(data, subsetIndex, featIdx) also selects feature
%   columns featIdx from data.X (used in the post-SIFTED density path).
if nargin < 3
    data.X = data.X(subsetIndex, :);
else
    data.X = data.X(subsetIndex, featIdx);
    data.featlabels = data.featlabels(featIdx);
end
data.Y            = data.Y(subsetIndex, :);
data.Xraw         = data.Xraw(subsetIndex, :);
data.Yraw         = data.Yraw(subsetIndex, :);
data.Ybin         = data.Ybin(subsetIndex, :);
data.beta         = data.beta(subsetIndex);
data.numGoodAlgos = data.numGoodAlgos(subsetIndex);
data.Ybest        = data.Ybest(subsetIndex);
data.P            = data.P(subsetIndex);
data.instlabels   = data.instlabels(subsetIndex);
if isfield(data, 'S')
    data.S = data.S(subsetIndex);
end
end
