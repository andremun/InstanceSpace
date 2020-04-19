try
    rootdir = 'E:/InstanceSpace_Classification/MATILDA_trial/';
    testingDir = strcat(rootdir, 'test/');
    status = mkdir(testingDir);
    if status == 0
        error('Can not create sub-directory to store testing results. ');
    else
        source_metadata = fullfile(rootdir, 'metadata_test.csv');
        source_model = fullfile(rootdir, 'model.mat');
        status = copyfile(source_metadata, testingDir);
        if status == 0
               error('Error: Please place metadata_test inside %s', rootdir );
        end
         status = copyfile(source_model, testingDir);
        if status == 0
               error('Error: could not find model.mat inside %s', rootdir );
        end

           source_options = fullfile(rootdir, 'options.json');
           status = copyfile(source_options, testingDir);
           if status == 0
               error('Error: could not find options.json inside %s', rootdir );
           end

        testIS(testingDir);  
    end
catch ME
    disp('EOF:ERROR');
    rethrow(ME)
end