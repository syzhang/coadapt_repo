function result_mult = run_cv_multiple_02()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cross validation main function
% - calls 'run_cv_session.m' 
% - runs multiple subjects/iterations to produce overall stats
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % subject list
    % sjList = 1:19;
    sjList = 1;

    % number of CV iterations
    numIteration = 5;

    result_mult = [];

    for dd = 1:2

        for ddd = 1:2
            PAR.ta = num2str(dd);
            PAR.te = num2str(ddd);

            for sj = 1:length(sjList)
                PAR.name = [num2str(sjList(sj), '%02d')];
                PAR.scanidx_suff = '_rest2test3'; % scan idx file suffix (restNtestN)

                for ii = 1:numIteration
                    % calls CV function
                    result = run_cv_session(PAR);
                    
                    % recording results
                    res_size = size(result);
                    roiNum = res_size(2);
                    result = [result, repmat(sj, size(result, 1), 1)];
                    result = [result, repmat(ii, size(result, 1), 1)];
                    result_mult = [result_mult; result];
                end % end ii

            end % end sj

        end % end ddd

    end % end dd

    % saving results
    SAVEDIR_RESUlTS = fullfile('..', 'res', 'cv_results');
    if ~exist(SAVEDIR_RESUlTS); mkdir(SAVEDIR_RESUlTS); end

    % writing header to csv
    cHeader = {'ROI', 'Accuracy', 'AUC', 'Sensitivity', 'Specificity', '#Voxels', 'Session', 'Train (day)', 'Test (day)', 'sj', 'iter'};
    commaHeader = [cHeader;repmat({','},1,numel(cHeader))]; % insert commaas
    commaHeader = commaHeader(:)';
    textHeader = cell2mat(commaHeader); %cHeader in text with commas

    % write header to file
    fid = fopen([SAVEDIR_RESUlTS, '/result_mult.txt'],'w'); 
    fprintf(fid,'%s\n',textHeader)
    fclose(fid)

    % write data to end of file
    dlmwrite([SAVEDIR_RESUlTS, '/result_mult.txt'], result_mult, '-append');

end % end function
