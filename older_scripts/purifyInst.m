function purifyInst(Net, epsilon, flag)
%% Description
% This script purifies benchmarks problems based on the similarity among them.
% The similarity of two benchmark problems is defined as the Euclidian
% distnace between their features.
% This script regards the similarity of tow benchmarks based on the
% similarity of the algorithm performances (APs) on them too; the simiarity of
% APs can be calculated as th Euclidian distance among the APs or as the
% similarity of their goodness in saolving the similar benchmarks.
% If two benchmarks are similar, then it checks the similarity among
% algorithm performances on these benchmarks.

% The arguments of this function are:
% Net: is the name of metadata file including the instances and their features
% and APs.
% epsilon: determines the treshhold of the similarity.
% flag: determines which similarity approach must be applied as follows;
% =========================================================================
% 'Ftr': similarity just based on the features
% 'Ftr&AP': both features and APs with Euclidian distance
% 'Ftr&Good': features with Euclidian distnace and APs based on the goodness
% 'Ftr&AP&Good': features with Euclidian dustance and APs
% with both Euclidian distance and goodness criterion
% =========================================================================
% In the cases of goodness criterion, we must define a differen threshold
% for the goddness of APs. This is ture about the
% second case too, but since the number of APs are less than the number of
% features, they are purified with different thresholds simlicitly;

% =========================================================================
% Code by: Hossein Alipour
%          School of Mathematics and Statistics
%          The University of Melbourne
%          Australia
%          2020
%          Email: h.alipour@unimelb.edu.au

%  Copyright: Hossein Alipour
%

%% Main loop

for k = 1:length(epsilon)
    
    Xbar = readtable(sprintf('%s.csv',Net));
    
    varlabels = Xbar.Properties.VariableNames;
    isfeat = strncmpi(varlabels,'feature_',8);
    isalgo = strncmpi(varlabels,'algo_',5);
    numFtr = sum(isfeat);
    numAlgo = sum(isalgo);
    X = Xbar{:,isfeat};
    Y = Xbar{:,isalgo};
    
    toPreclude = false((size(X,1)),1);
    ViSA_idx = false((size(X,1)),1);
    Dissimilar_idx = true((size(X,1)),1);
    Tprcl = 0; % counter of the precluded instances
    
    
    ViSA_D = 0; % A counter of the violation from Euclidian distance for APs
    ViSA_Good = 0; % A counter of the violation from the goodness criterion for APs
    
    
    
    %% The main nested loop to preclude the similar instances according to a
    %  given tolearnce epsilon, preferably from the interval (0, 1).
    Dissimilar_idx(1) = true;
    
    for i=1:size(X,1)
        if (~toPreclude(i))
            for j=i+1:size(X,1)
                if (~toPreclude(j))
                    dist = sqrt(sum((X(i,:) - X(j,:)) .^ 2)); % Euclidean distance
                    if (dist <= epsilon(k))
                        Dissimilar_idx(j) = false;
                        switch flag
                            case 'Ftr'
                                toPreclude(j) = true;
                            case 'Ftr&AP'
                                dist_AP = sqrt(sum((Y(i,:) - Y(j,:)) .^ 2));
%                                 if (dist_AP <= sqrt(numAlgo/numFtr)*epsilon(k)*(1+0.26)) % 0.26 is the value of CPUN (CPU Noise)
                                if (dist_AP <= sqrt(numAlgo/numFtr)*epsilon(k)) 
                                    toPreclude(j) = true;
                                    ViSA_idx(j) = false;
                                else
                                    ViSA_idx(j) = true;
                                end
                            case 'Ftr&Good'
                                if (Ybin(i,:) == Ybin(j,:))
                                    toPreclude(j) = true; 
                                    ViSA_idx(j) = false;
                                else
                                    ViSA_idx(j) = true;
                                end
                            case 'Ftr&AP&Good'
                                if (Ybin(i,:) == Ybin(j,:))
                                    dist_AP = sqrt(sum((Y(i,:) - Y(j,:)) .^ 2));
%                                     if (dist_AP <= sqrt(numAlgo/numFtr)*epsilon(k)*(1+0.26)) % epsilon could be different from that given in the function.
                                    if (dist_AP <= sqrt(numAlgo/numFtr)*epsilon(k))
                                        toPreclude(j) = true;
                                        ViSA_idx(j) = false;
                                    else
                                        ViSA_idx(j) = true;
                                    end
                                else
                                    ViSA_idx(j) = true;
                                end
                            otherwise
                                disp('Invalid flag!')
                        end
                    end
                end
            end
        end
    end
       
    %  Xbar = sortrows(Xbar); % This command is useful if you want to apply a random purification.
    
    PurifiedInst = Xbar;
    PrecludedInst = Xbar;
    DissimilarPurInst = Xbar;
    ViSAPurInst = Xbar;
    
    PurifiedInst(toPreclude,:) = [];
    PrecludedInst(~toPreclude,:) = [];
    DissimilarPurInst(~Dissimilar_idx,:) = [];
    ViSAPurInst(~ViSA_idx,:) = [];
    
    switch flag
        case 'Ftr'
            writetable(PurifiedInst,sprintf('Pur_%s_%s_Dist_%.3f.csv',flag, Net, epsilon(k)));
            writetable(PrecludedInst,sprintf('Prec_%s_%s_Dist_%.3f.csv',flag, Net, epsilon(k)));
        case 'Ftr&AP'
            writetable(PurifiedInst,sprintf('Pur_%s_%s_Dist_%.3f.csv',flag, Net, epsilon(k)));
            writetable(PrecludedInst,sprintf('Prec_%s_%s_Dist_%.3f.csv',flag, Net, epsilon(k)));
        case 'Ftr&Good'
            writetable(PurifiedInst,sprintf('Pur_%s_%s_G_%.2f_Dist_%.3f.csv',flag, Net,  opts.perf.epsilon, epsilon(k)));
            writetable(PrecludedInst,sprintf('Prec_%s_%s_G_%.2f_Dist_%.3f.csv',flag, Net,  opts.perf.epsilon, epsilon(k)));
            writetable(DissimilarPurInst,sprintf('Dissimilar_%s_%s_G_%.2f_Dist_%.3f.csv',flag, Net,  opts.perf.epsilon, epsilon(k)));
            writetable(ViSAPurInst,sprintf('ViSA_%s_%s_G_%.2f_Dist_%.3f.csv',flag, Net,  opts.perf.epsilon, epsilon(k)));
        case 'Ftr&AP&Good'
            writetable(PurifiedInst,sprintf('Pur_%s_%s_G_%.2f_Dist_%.3f.csv',flag, Net,  opts.perf.epsilon, epsilon(k)));
            writetable(PrecludedInst,sprintf('Prec_%s_%s_G_%.2f_Dist_%.3f.csv',flag, Net,  opts.perf.epsilon, epsilon(k)));
            writetable(DissimilarPurInst,sprintf('Dissimilar_%s_%s_G_%.2f_Dist_%.3f.csv',flag, Net,  opts.perf.epsilon, epsilon(k)));
            writetable(ViSAPurInst,sprintf('ViSA_%s_%s_G_%.2f_Dist_%.3f.csv',flag, Net,  opts.perf.epsilon, epsilon(k)));
    end

    %%
    clear Xbar;
    clear X;
    clear Y;
    clear Ybin;
    clear toPreclude;
    fclose('all');
end

