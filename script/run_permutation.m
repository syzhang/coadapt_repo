function result_mult = run_permutation()
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % run permutation test
    % - runs multiple subjects/iterations to produce statistical significance of cv results
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        addpath('../SLR1.51/')
        PAR.TOPDIR = fileparts(pwd);

        % subject list
        sjList = 1:19;
    
        % number of CV iterations
        numIteration = 1;
    
        result_mult = [];

        % saving results
        SAVEDIR_RESUlTS = fullfile('..', 'res', 'permutation_results');
        if ~exist(SAVEDIR_RESUlTS); mkdir(SAVEDIR_RESUlTS); end
    
        % writing header to csv
        cHeader = {'cv_iter', 'test_acc', 'p05', 'Train (day)', 'Test (day)', 'sj'};
        commaHeader = [cHeader;repmat({','},1,numel(cHeader))]; % insert commaas
        commaHeader = commaHeader(:)';
        textHeader = cell2mat(commaHeader); %cHeader in text with commas

        % write header to file
        fid = fopen([SAVEDIR_RESUlTS, '/permutation_results.txt'],'w'); 
        fprintf(fid,'%s\n',textHeader);
        fclose(fid);

        for dd = 1:2    
            for ddd = 1:2
                PAR.ta = num2str(dd);
                PAR.te = num2str(ddd);
    
                for sj = 1:length(sjList)
                    % which roi to use
                    PAR.roiuse = 1; % insula only
                    PAR.roiName = 'Insula';
            
                    % file directory
                    PAR.name = [num2str(sjList(sj), '%02d')];
                    PAR.scanidx_suff = '_rest2test3'; % scan idx file suffix (restNtestN)
                    FEATUREDIR = [PAR.TOPDIR, '/dat/', PAR.name, filesep, 'features_scantrend/'];
                    PARADIR = [PAR.TOPDIR, '/dat/', PAR.name, filesep, 'params/'];
            
                    % load data
                    if PAR.ta == '1' && PAR.te == '1' % training with day 1, testing with day 1
                        PAR.tadir = [FEATUREDIR, 'tmpl_d1_d1/'];
                        PAR.tedir = PAR.tadir;
                    elseif PAR.ta == '1' && PAR.te == '2' % training with day 1, testing with day 2
                        PAR.tadir = [FEATUREDIR, 'tmpl_d1_d1/'];
                        PAR.tedir = [FEATUREDIR, 'tmpl_d1_d2/'];
                    elseif PAR.ta == '2' && PAR.te == '2' % training with day 2, testing with day 2
                        PAR.tadir = [FEATUREDIR, 'tmpl_d2_d2/'];
                        PAR.tedir = PAR.tadir;
                    elseif PAR.ta == '2' && PAR.te == '1' % training with day 2, testing with day 1
                        PAR.tadir = [FEATUREDIR, 'tmpl_d2_d2/'];
                        PAR.tedir = [FEATUREDIR, 'tmpl_d2_d1/'];
                    end

                    result = run_permutation_slr(PAR);
                    % result = run_permutation_svm(PAR);

                    % recording results
                    res_size = size(result);
                    roiNum = res_size(2);
                    result = [result, repmat(dd, size(result, 1), 1)];
                    result = [result, repmat(ddd, size(result, 1), 1)];
                    result = [result, repmat(sj, size(result, 1), 1)];
                    % result_mult = [result_mult; result];

                    % write data to end of file
                    % dlmwrite([SAVEDIR_RESUlTS, '/permutation_results.txt'], result_mult, '-append');
                    dlmwrite([SAVEDIR_RESUlTS, '/permutation_results.txt'], result, '-append');

                end % end sj
    
            end % end ddd
    
        end % end dd
    

    
    end % end function
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function DAm_rec = run_permutation_slr(PAR)
%% classification:
    p = [0.05];%, 0.01]; % p value
    permutations = 1/min(p);
    display(['-> Subject: ' num2str(PAR.name)])

    % setting parameters and load data files
    DAT_ta = load([PAR.tadir, num2str(PAR.roiuse, '%03d'), filesep, 'feature_', PAR.roiName]); % load 'feature', 'label'
    DAT_te = load([PAR.tedir, num2str(PAR.roiuse, '%03d'), filesep, 'feature_', PAR.roiName]); % load 'feature', 'label'

    
    nfold = 10;
    % balance data first
    [X_train, Y_train] = feat_unroll_balance(DAT_ta);
    [X_test, Y_test] = feat_unroll_balance(DAT_te);
    % not balancing
    % [X_train, Y_train] = feat_unroll(DAT_ta);
    % [X_test, Y_test] = feat_unroll(DAT_te);
    if length(Y_train)<=length(Y_test)
        cv = cvpartition(Y_train, 'KFold', nfold);
    else
        cv = cvpartition(Y_test, 'KFold', nfold);
    end % if

    DAm_rec = [];
    % load feature
    for f=1:nfold
        if PAR.tadir == PAR.tedir
            X_ta = X_train(cv.training(f), :);
            Y_ta = Y_train(cv.training(f));
            X_te = X_train(cv.test(f), :);
            Y_te = Y_train(cv.test(f));
        else
            X_ta = X_train(cv.training(f), :);
            Y_ta = Y_train(cv.training(f));
            X_te = X_test(cv.test(f), :);
            Y_te = Y_test(cv.test(f));
        end % if
        % calculate 
        percent = SLR_nonCV(X_ta, Y_ta, X_te, Y_te);
    
        for k=1:permutations
            display(['  -> Permutations :' num2str(k)])
            rng_idx = randi([6000 8000],1,1);
            rng(k,'twister')
            
            y_perm = Y_te(randperm(length(Y_te)));
            percent_perm(k) = SLR_nonCV(X_ta, Y_ta, X_te, y_perm);
        end
    
        ChanceLevel = 100/length(unique(Y_te));
        PercDatasetsM = mean(percent_perm);
    
        CurPercPerm = sort(percent_perm);
        [u,~] = find(percent < CurPercPerm);
        for k=1:length(p)
            LegStr{k} = ['p<' num2str(p(k))];
            Nperm = p(k)*permutations;
            DA(k) = CurPercPerm(end-Nperm+1);
        end
    
        % DAm = squeeze(mean(DA,1))
        DAm = [f, percent, squeeze(mean(DA,1))]
        % fprintf('\nStatistical chance level: %.3f.\n', DAm)
        DAm_rec = [DAm_rec; DAm];
    end % f

end % end function

function DAm_rec = run_permutation_svm(PAR)
    %% classification:
        p = [0.05];%, 0.01]; % p value
        permutations = 1/min(p);
        display(['-> Subject: ' num2str(PAR.name)])
    
        % setting parameters and load data files
        DAT_ta = load([PAR.tadir, num2str(PAR.roiuse, '%03d'), filesep, 'feature_', PAR.roiName]); % load 'feature', 'label'
        DAT_te = load([PAR.tedir, num2str(PAR.roiuse, '%03d'), filesep, 'feature_', PAR.roiName]); % load 'feature', 'label'
        
        nfold = 10;
        % balance data first
        [X_train, Y_train] = feat_unroll_balance(DAT_ta);
        [X_test, Y_test] = feat_unroll_balance(DAT_te);
        % not balancing
        % [X_train, Y_train] = feat_unroll(DAT_ta);
        % [X_test, Y_test] = feat_unroll(DAT_te);
        cv = cvpartition(Y_train, 'KFold', nfold);
        
        DAm_rec = [];
        % load feature
        for f=1:nfold
            if PAR.tadir == PAR.tedir
                X_ta = X_train(cv.training(f), :);
                Y_ta = Y_train(cv.training(f));
                X_te = X_train(cv.test(f), :);
                Y_te = Y_train(cv.test(f));
            else
                X_ta = X_train(cv.training(f), :);
                Y_ta = Y_train(cv.training(f));
                X_te = X_test(cv.test(f), :);
                Y_te = Y_test(cv.test(f));
            end % if
            % calculate 
            percent = SVM_nonCV(X_ta, Y_ta, X_te, Y_te);
        
            for k=1:permutations
                display(['  -> Permutations :' num2str(k)])
                rng_idx = randi([6000 8000],1,1);
                rng(k,'twister')
                
                y_perm = Y_te(randperm(length(Y_te)));
                percent_perm(k) = SVM_nonCV(X_ta, Y_ta, X_te, y_perm);
            end
        
            ChanceLevel = 100/length(unique(Y_te));
            PercDatasetsM = mean(percent_perm);
        
            CurPercPerm = sort(percent_perm);
            [u,~] = find(percent < CurPercPerm);
            for k=1:length(p)
                LegStr{k} = ['p<' num2str(p(k))];
                Nperm = p(k)*permutations;
                DA(k) = CurPercPerm(end-Nperm+1);
            end
        
            % DAm = squeeze(mean(DA,1))
            DAm = [f, percent, squeeze(mean(DA,1))]
            % fprintf('\nStatistical chance level: %.3f.\n', DAm)
            DAm_rec = [DAm_rec; DAm];
        end % f
    
    end % end function    

function acc2 = SVM_nonCV(X_ta, Y_ta, X_te, Y_te)
    % train svm
    SVMModel = fitcsvm(X_ta, Y_ta, 'ClassNames', [-1, 1]);
    [Y_pred, score] = predict(SVMModel, X_te);
    errTable_te = confusionmat(Y_te, Y_pred, 'Order', [-1, 1]);

    % accuracy
    acc2 = calc_percor(errTable_te) ./ 100;
end % function

function [X, Y] = feat_unroll_balance(DAT)
    % balancing high/low pain trials in each session (run)

    feat_all = [];
    label_all = [];
    for csess = 1:length(DAT.feature)
        dum_feat = DAT.feature{csess};
        dum_label = DAT.label(:, csess);

        % concat all session
        feat_all = [feat_all; dum_feat];
        label_all = [label_all; dum_label];
    end % i
    % balance
    h_num = sum(label_all == 1);
    l_num = sum(label_all == -1);
    min_num = min([h_num, l_num]);
    max_num = max([h_num, l_num]);
    
    % randomly choose samples with replacement (upsampling)
    h_idx = find(label_all == 1);
    l_idx = find(label_all == -1);
    % h_chosen = datasample(h_idx, min_num, 'Replace', false);
    % l_chosen = datasample(l_idx, min_num, 'Replace', false);
    h_chosen = datasample(h_idx, max_num, 'Replace', true);
    l_chosen = datasample(l_idx, max_num, 'Replace', true);
    chosen_idx = [h_chosen; l_chosen];
    chosen_idx = chosen_idx(randperm(length(chosen_idx)));

    % chosen trials in session
    chosen_feat = feat_all(chosen_idx, :);
    chosen_label = label_all(chosen_idx, :);

    X = chosen_feat;
    Y = chosen_label;
end % end function

function [X, Y] = feat_unroll(DAT)
    % unroll feature and label
    % data struct
    feat_all = [];
    label_all = [];
    for csess = 1:length(DAT.feature)
        dum_feat = DAT.feature{csess};
        dum_label = DAT.label(:, csess);

        % concat all session
        feat_all = [feat_all; dum_feat];
        label_all = [label_all; dum_label];
    end % i
    X = feat_all;
    Y = label_all;
    
end % end function
    
function acc2 = SLR_nonCV(X_ta, Y_ta, X_te, Y_te)
    % train SLR decoder
    [ww_f, ~, ~, errTable_te, ~, ~, ~, ~, ~, t_test_est] = biclsfy_slrvar(X_ta, Y_ta, X_te, Y_te, ...
        'nlearn', 300, 'mean_mode', 'none', 'scale_mode', 'none', 'invhessian', 0, 'usebias', 1);
    % slr_results = sum(diag(errTable_te));
    % W_data = ww_f;
    % performance = sum(slr_results) / size(X_te, 1);
    acc2 = sum(t_test_est==Y_te)/length(Y_te);
    fprintf('\nCalc_label == pred_label perc: %.3f.\n', acc2)

    % saving weights with roi name
    % save([SUBSLRDIR ['Selected_W_nonCV_', PAR.roiName]], 'W_data');

end % end function
