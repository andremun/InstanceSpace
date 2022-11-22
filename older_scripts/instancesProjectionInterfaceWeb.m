function instancesProjectionInterfaceWeb(rootdir)
    try
        testIS(rootdir);
    catch ME
        disp('EOF:ERROR');
        rethrow(ME)
    end
end