function instancesProjectionInterfaceWeb(rootdir)
% function instancesProjectionInterfaceWeb()
% rootdir = '/servers/matilda-users/neelofar/data_projection-black_box_single_objective_1_7/';
    try
        testIS(rootdir);
    catch ME
        disp('EOF:ERROR');
        rethrow(ME)
    end
end