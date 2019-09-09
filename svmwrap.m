function [ytest,svm] = svmwrap(xtrain,ytrain,xtest,params,weight,flag)

% global params
ytrain = double(ytrain);
prior = 1-mean(bsxfun(@eq,ytrain,[1 2]));
prior(2) = (weight>1e-2).*prior(2);
prior = 100.*prior./sum(prior);
command = ['-s 0 -t 2 -q -b 1 -c ' num2str(params(1)) ' -g ' num2str(params(2))...
                           ' -w1 ' num2str(prior(1)) ' -w2 ' num2str(prior(2))];
svm = svmtrain(ytrain, xtrain, command); %#ok<SVMTRAIN>
[ytest,~,prob] = svmpredict(randi([1 2],size(xtest,1),1), xtest, svm, '-q'); 
if flag
   ytest = [ytest prob]; 
end

end
