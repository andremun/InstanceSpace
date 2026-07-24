function printOptions(opts, prefix)
% printOptions  Compact one-line-per-setting console dump of an opts struct.
%
%   printOptions(opts)
%
%   Recurses into nested option structs (e.g. opts.pilot, opts.sifted) so
%   every leaf setting prints as "section.field = value", replacing the
%   old fieldnames()+disp() loop that printed each top-level opts section
%   using MATLAB's verbose default struct display (~80+ lines for a full
%   options.json).
if nargin < 2, prefix = ''; end
fields = fieldnames(opts);
for i = 1:numel(fields)
    f = fields{i};
    v = opts.(f);
    key = [prefix f];
    if isstruct(v) && isscalar(v)
        printOptions(v, [key '.']);
    else
        fprintf('  %-28s %s\n', key, formatOptionValue(v));
    end
end
end
% =========================================================================
function s = formatOptionValue(v)
if ischar(v)
    s = ['''' v ''''];
elseif isstring(v) && isscalar(v)
    s = ['''' char(v) ''''];
elseif iscell(v) && all(cellfun(@ischar, v))
    s = ['{' strjoin(v, ', ') '}'];
elseif islogical(v) && isscalar(v)
    s = mat2str(v);
elseif isnumeric(v) && isempty(v)
    s = '[]';
elseif isnumeric(v) && isscalar(v)
    s = num2str(v);
elseif isnumeric(v)
    s = mat2str(v);
else
    s = class(v);
end
end
