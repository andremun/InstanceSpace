function ytest = svmwrap(xtrain,ytrain,xtest,params,flag)

% global params
prior = 1-mean(bsxfun(@eq,ytrain,[1 2]));
prior = 100.*prior./sum(prior);
command = ['-s 0 -t 2 -q -b 1 -c ' num2str(params(1)) ' -g ' num2str(params(2))...
                           ' -w1 ' num2str(prior(1)) ' -w2 ' num2str(prior(2))];
[ytest,~,prob] = svmpredict(randi([1 2],size(xtest,1),1), ...
                            xtest, svmtrain(ytrain, xtrain, command), '-q'); %#ok<SVMTRAIN>
if flag
   ytest = [ytest prob]; 
end
