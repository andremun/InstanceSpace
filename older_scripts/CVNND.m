function CVNND_Single(Net)

%% Descriptions
% CVNND stands for cofficient varioation of the nearest neighbor distances (NNDs).
% For each point in the instance space, this function calculates its distance
% form its nearest neighbor point. Then, claculates the coefficinet variation
% of all such distances. CVNND is used as a measure of distribution
% uniformuity (evenness) for non-classified samples. In fact, uniformity
% can be calculated as 1-CVNND.

% Depending on the flag, we might regard NND for both features and
% algorithm prefromnces (APs). In such cases, NNDs of features and APs are
% calculated with different scales, which concide with the scalses that
% data sets are purified with the corresponding flags in
% purifyInstIS_Cmp(Net, epsilon, flag) function.

% The argument of this function is:
% Net: is the name of metadata file including the instances and their features
% =========================================================================
% Code by: Hossein Alipour
%          School of Mathematics and Statistics
%          The University of Melbourne
%          Australia
%          2020
%          Email: h.alipour@unimelb.edu.au

%  Copyright: Hossein Alipour
%

%% 

TblHeader_CVNND = {'InstNumb'  'CV_All' 'Uniformity_All' };

textHeader_CVNND = strjoin(TblHeader_CVNND, ',');

NewFileName = sprintf('CVNND_%s.csv', Net);
fid = fopen(sprintf('%s',NewFileName),'w');
fprintf(fid,'%s\n',textHeader_CVNND);
fclose(fid);
clear fid*

 Xbar = readtable(sprintf('%s.csv', Net));

%%
    
    
    varlabels = Xbar.Properties.VariableNames;
    isfeat = strncmpi(varlabels,'feature_',8);
    isalgo = strncmpi(varlabels,'algo_',5);
    numFtr = sum(isfeat);
    numAlgo = sum(isalgo);
    X = Xbar{:,isfeat};
    Y = Xbar{:,isalgo};
    
    
    opts.norm.flag = true;
    [X, Y, model.norm] = autoNormalize(X, Y, opts.norm);
    
    
    
    
    %% The loop to calculate CVNND and Uniformity (Evenness);
    
    nearestDist = zeros(length(X),1);
    for ii=1:length(X)
        Xtmp = X;
        Xtmp(ii,:)=[];
        allDiff = bsxfun(@minus,Xtmp,X(ii,:));
        EucDist = cellfun(@norm,num2cell(allDiff,2));
        nearestDist(ii) = min(EucDist);
    end
    
    CV = std(nearestDist)/mean(nearestDist) % coefficient of variation of the nearest neighbor distance
    Uniformity = 1- CV

    

    %% Wrtie data on the table
%     TblHeader_CVNND = {'Epsilon' 'InstNumb' 'NoViSANumb' 'ViSANumb' 'RViSA' 'CV_All' 'Uniformity_All' 'CV_NoViSA' 'Uniformity_NoViSA' };
    Current_data_CVNND = [ length(X), CV, Uniformity];
    fid = fopen(sprintf('%s',NewFileName),'a');
    fprintf(fid,'%f,%f,%f\n', Current_data_CVNND);
    fclose(fid);
    clear fid*

    %%
    clear Xbar;
    clear X;
    clear XbarNoViSA;
    clear XNoViSA;
    clear XbarViSA
end



