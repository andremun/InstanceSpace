function exampleWeb(rootdir)

try
    model = trainIS(rootdir);
    testIS(rootdir);
catch ME
    disp('EOF:ERROR');
    rethrow(ME)
end
end
