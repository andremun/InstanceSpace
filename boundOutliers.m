function [X, out] = boundOutliers(X)

disp('-> Removing extreme outliers from the feature values.');
out.medval = nanmedian(X, 1);
out.iqrange = iqr(X, 1);
out.hibound = out.medval + 5.*out.iqrange;
out.lobound = out.medval - 5.*out.iqrange;
himask = bsxfun(@gt,X,out.hibound);
lomask = bsxfun(@lt,X,out.lobound);
X = X.*~(himask | lomask) + bsxfun(@times,himask,out.hibound) + ...
                            bsxfun(@times,lomask,out.lobound);

end