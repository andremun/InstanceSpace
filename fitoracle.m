function out = fitoracle(X,Y,opts)

% global params
Y = double(Y)+1;
nalgos = size(Y,2);
out.paramgrid = sortrows(2.^(opts.maxcvgrid.*lhsdesign(opts.cvgrid,2) + ...
                             opts.mincvgrid));  % Cross-validation grid
out.cvmcr = NaN.*ones(opts.cvgrid,nalgos);
out.Yhat = NaN.*Y;
out.paramidx = NaN.*ones(1,nalgos);
out.oracle = cell(nalgos,1);

for i=1:nalgos
    for j=1:opts.cvgrid
        out.cvmcr(j,i) = crossval('mcr', X, Y(:,i),...
                                'Kfold', opts.cvfolds,...
                                'Predfun',@(xtrain,ytrain,xtest) svmwrap(xtrain,...
                                                                         ytrain,...
                                                                         xtest,...
                                                                         out.paramgrid(j,:)));
    end
    [~,out.paramidx(i)] = min(out.cvmcr(:,i));
    prior = 1-mean(bsxfun(@eq,Y(:,i),[0 1]));
    prior = 100.*prior./sum(prior);
    command = ['-s 0 -t 2 -q -b 1 -c ' num2str(out.paramgrid(out.paramidx(i),1)) ...
                                ' -g ' num2str(out.paramgrid(out.paramidx(i),2))...
                               ' -w0 ' num2str(prior(1)) ' -w1 ' num2str(prior(2))];
    out.oracle{i} = svmtrain(Y(:,i), X, command); %#ok<SVMTRAIN>
    out.Yhat(:,i) = svmpredict(Y(:,i), X, out.oracle{i}, '-q');
end
out.mcr = mean(Y~=out.Yhat,1);