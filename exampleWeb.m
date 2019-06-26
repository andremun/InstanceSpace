function exampleWeb(rootdir)

try
    matilda(rootdir);
catch ME
    disp('EOF:ERROR');
    rethrow(ME)
end
end
