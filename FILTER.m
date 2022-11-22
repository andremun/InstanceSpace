function [subsetIndex,isDissimilar,isVISA] = FILTER(X,Y,Ybin,opts)

[ninst,nalgos] = size(Y);
nfeats = size(X,2);

subsetIndex = false(ninst,1);
isDissimilar = true(ninst,1);
isVISA = false(ninst,1);
gamma = sqrt(nalgos/nfeats)*opts.mindistance;

for ii=1:ninst
    if ~subsetIndex(ii)
        for jj=ii+1:ninst
            if ~subsetIndex(jj)
                Dx = pdist2(X(ii,:),X(jj,:));
                Dy = pdist2(Y(ii,:),Y(jj,:));
                Db = all(Ybin(ii,:) & Ybin(jj,:));
                if  Dx <= opts.mindistance
                    isDissimilar(jj) = false;
                    switch opts.type
                        case 'Ftr'
                            subsetIndex(jj) = true;
                        case 'Ftr&AP'
                            if Dy <= gamma
                                subsetIndex(jj) = true;
                                isVISA(jj) = false;
                            else
                                isVISA(jj) = true;
                            end
                        case 'Ftr&Good'
                            if Db
                                subsetIndex(jj) = true;
                                isVISA(jj) = false;
                            else
                                isVISA(jj) = true;
                            end
                        case 'Ftr&AP&Good'
                            if Db
                                if Dy <= gamma
                                    subsetIndex(jj) = true;
                                    isVISA(jj) = false;
                                else
                                    isVISA(jj) = true;
                                end
                            else
                                isVISA(jj) = true;
                            end
                        otherwise
                            disp('Invalid flag!')
                    end
                end
            end
        end
    end
end

% Assess the uniformity of the data
D = squareform(pdist(X(~subsetIndex,:)));
ninst = size(D,1);
D(eye(ninst,'logical')) = NaN;
nearest = min(D,[],2,'omitnan');
model.data.unif = 1-(std(nearest)./mean(nearest));
disp(['Uniformity of the instance subset: ' num2str(model.data.unif,4)]);

end