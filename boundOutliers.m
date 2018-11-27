function [X, out] = boundOutliers(X, opts)

out.medval = median(X, 1);
out.iqrange = iqr(X, 1);
if opts.flag
    out.hibound = out.medval + 5.*out.iqrange;
    out.lobound = out.medval - 5.*out.iqrange;
    himask = bsxfun(@gt,X,out.hibound);
    lomask = bsxfun(@lt,X,out.lobound);
    X = X.*~(himask | lomask) + bsxfun(@times,himask,out.hibound) + ...
                                bsxfun(@times,lomask,out.lobound);
else
    out.hibound = max(X, [], 1);
    out.lobound = min(X, [], 1);
end

end