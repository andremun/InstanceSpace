function [X,Y,Ybest,Ybin,P,numGoodAlgos,beta,out] = PRELIM(X,Y,opts)

Yraw = Y;
nalgos = size(Y,2);
% -------------------------------------------------------------------------
% Determine whether the performance of an algorithm is a cost measure to
% be minimized or a profit measure to be maximized. Moreover, determine
% whether we are using an absolute threshold as good peformance (the
% algorithm has a performance better than the threshold) or a relative
% performance (the algorithm has a performance that is similar that the
% best algorithm minus a percentage).
disp('-------------------------------------------------------------------------');
disp('-> Calculating the binary measure of performance');
msg = '-> An algorithm is good if its performace is ';
if opts.MaxPerf
    Yaux = Y;
    Yaux(isnan(Yaux)) = -Inf;
    [Ybest,P] = max(Yaux,[],2);
    if opts.AbsPerf
        Ybin = Yaux>=opts.epsilon;
        msg = [msg 'higher than ' num2str(opts.epsilon)];
    else
        Ybest(Ybest==0) = eps;
        Y(Y==0) = eps;
        Y = 1-bsxfun(@rdivide,Y,Ybest);
        Ybin = (1-bsxfun(@rdivide,Yaux,Ybest))<=opts.epsilon;
        msg = [msg 'within ' num2str(round(100.*opts.epsilon)) '% of the best.'];
    end
else
    Yaux = Y;
    Yaux(isnan(Yaux)) = Inf;
    [Ybest,P] = min(Yaux,[],2);
    if opts.AbsPerf
        Ybin = Yaux<=opts.epsilon;
        msg = [msg 'less than ' num2str(opts.epsilon)];
    else
        Ybest(Ybest==0) = eps;
        Y(Y==0) = eps;
        Y = bsxfun(@rdivide,Y,Ybest)-1;
        Ybin = (bsxfun(@rdivide,Yaux,Ybest)-1)<=opts.epsilon;
        msg = [msg 'within ' num2str(round(100.*opts.epsilon)) '% of the best.'];
    end
end
disp(msg);
% -------------------------------------------------------------------------
% Testing for ties. If there is a tie in performance, we pick an algorithm
% at random.
bestAlgos = bsxfun(@eq,Yraw,Ybest);
multipleBestAlgos = sum(bestAlgos,2)>1;
aidx = 1:nalgos;
for i=1:size(Y,1)
    if multipleBestAlgos(i)
        aux = aidx(bestAlgos(i,:));
        P(i) = aux(randi(length(aux),1)); % Pick one at random
    end
end
disp(['-> For ' num2str(round(100.*mean(multipleBestAlgos))) '% of the instances there is ' ...
      'more than one best algorithm. Random selection is used to break ties.']);
numGoodAlgos = sum(Ybin,2);
beta = numGoodAlgos>(opts.betaThreshold*nalgos);

disp('=========================================================================');
disp('-> Auto-pre-processing.');
disp('=========================================================================');
if opts.bound
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

if opts.norm
    nfeats = size(X,2);
    nalgos = size(Y,2);
    disp('-> Auto-normalizing the data using Box-Cox and Z transformations.');
    out.minX = min(X,[],1);
    X = bsxfun(@minus,X,out.minX)+1;
    out.lambdaX = zeros(1,nfeats);
    out.muX = zeros(1,nfeats);
    out.sigmaX = zeros(1,nfeats);
    for i=1:nfeats
        aux = X(:,i);
        idx = isnan(aux);
        [aux, out.lambdaX(i)] = boxcox(aux(~idx));
        [aux, out.muX(i), out.sigmaX(i)] = zscore(aux);
        X(~idx,i) = aux;
    end

    out.minY = min(Y(:));
    Y = (Y-out.minY)+eps;
    out.lambdaY = zeros(1,nalgos);
    out.muY = zeros(1,nalgos);
    out.sigmaY = zeros(1,nalgos);
    for i=1:nalgos
        aux = Y(:,i);
        idx = isnan(aux);
        [aux, out.lambdaY(i)] = boxcox(aux(~idx));
        [aux, out.muY(i), out.sigmaY(i)] = zscore(aux);
        Y(~idx,i) = aux;
    end 
end

end
