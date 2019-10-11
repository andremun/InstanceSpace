function exampleWeb(rootdir)

try
    model = trainIS(rootdir);
catch ME
    disp('EOF:ERROR');
    rethrow(ME)
end
end
