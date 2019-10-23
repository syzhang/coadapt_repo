function result = run_cv_session(PAR)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cross validation fucntions
% - check decoder performance by running CV 
% - run through all combinations of training/testing days
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    PAR.TOPDIR = fileparts(pwd);
    result = [];
    numROIs = 10; % number of ROI masks

    for ii = 1:numROIs
        % which roi to use
        PAR.roiuse = ii;

        % file directory
        FEATUREDIR = [PAR.TOPDIR, '/dat/', PAR.name, filesep, 'features_scantrend/'];
        PARADIR = [PAR.TOPDIR, '/dat/', PAR.name, filesep, 'params/'];

        % load data
        if PAR.ta == '1' && PAR.te == '1' % training with day 1, testing with day 1
            PAR.tadir = [FEATUREDIR, 'tmpl_d1_d1/'];
            PAR.tedir = PAR.tadir;
            PAR.scan_info_ta = [PARADIR, 'scan_idx/scan_idx', PAR.scanidx_suff, '.mat'];
            PAR.scan_info_te = PAR.scan_info_ta;
        elseif PAR.ta == '1' && PAR.te == '2' % training with day 1, testing with day 2
            PAR.tadir = [FEATUREDIR, 'tmpl_d1_d1/'];
            PAR.tedir = [FEATUREDIR, 'tmpl_d1_d2/'];
            PAR.scan_info_ta = [PARADIR, 'scan_idx/scan_idx', PAR.scanidx_suff, '.mat'];
            PAR.scan_info_te = [PARADIR, 'scan_idx_rt/scan_idx', PAR.scanidx_suff, '.mat'];
        elseif PAR.ta == '2' && PAR.te == '2' % training with day 2, testing with day 2
            PAR.tadir = [FEATUREDIR, 'tmpl_d2_d2/'];
            PAR.tedir = PAR.tadir;
            PAR.scan_info_ta = [PARADIR, 'scan_idx_rt/scan_idx', PAR.scanidx_suff, '.mat'];
            PAR.scan_info_te = PAR.scan_info_ta;
        elseif PAR.ta == '2' && PAR.te == '1' % training with day 2, testing with day 1
            PAR.tadir = [FEATUREDIR, 'tmpl_d2_d2/'];
            PAR.tedir = [FEATUREDIR, 'tmpl_d2_d1/'];
            PAR.scan_info_ta = [PARADIR, 'scan_idx_rt/scan_idx', PAR.scanidx_suff, '.mat'];
            PAR.scan_info_te = [PARADIR, 'scan_idx/scan_idx', PAR.scanidx_suff, '.mat'];
        end

        % load features
        featFile_ta = dir([PAR.tadir, num2str(PAR.roiuse, '%03d'), filesep, 'feature_*.mat']);
        featFile_te = dir([PAR.tedir, num2str(PAR.roiuse, '%03d'), filesep, 'feature_*.mat']);
        temp = strsplit(featFile_ta.name(1:end - 4), '_');
        PAR.roiName = temp{end}

        % run main check function
        [acc, auc, sen, spe, nv] = SVM_param_balance_tate(PAR);

        % recording results
        result_roi = repmat(ii, size(acc, 1), size(acc, 2));
        result_acc = acc;
        result_auc = auc;
        result_sen = sen;
        result_spe = spe;
        result_nv = nv;
        result_sess = [1:6]';
        result_ta = repmat(str2num(PAR.ta), size(acc, 1), size(acc, 2));
        result_te = repmat(str2num(PAR.te), size(acc, 1), size(acc, 2));
        % result = ['ROI','Accuracy','AUC','Sensitivity','Specificity','#Voxels','Session','Train (day)','Test (day)'];
        result_all = [result_roi, result_acc, result_auc, result_sen, result_spe, result_nv, result_sess, result_ta, result_te];
        result = [result; result_all];
    
    end % end ii

end % end function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [acc, auc, sen, spe, nv] = SVM_10folds(DATA_ta, DATA_te)
    % 10 fold CV with SVM

    X_ta = DATA_ta.X;
    Y_ta = DATA_ta.Y;
    X_te = DATA_te.X;
    Y_te = DATA_te.Y;

    % balance number of trials/dataset
    ta_p = sum(Y_ta == 1);
    te_p = sum(Y_te == 1);
    ta_n = sum(Y_ta == -1);
    te_n = sum(Y_te == -1);
    n_trial = min([ta_p, te_p, ta_n, te_n]);

    % resample trials
    idx_ta_p = randsample(find(Y_ta == 1), n_trial);
    idx_ta_n = randsample(find(Y_ta == -1), n_trial);
    rX_ta = [X_ta(idx_ta_p, :); X_ta(idx_ta_n, :)];
    rY_ta = [Y_ta(idx_ta_p, :); Y_ta(idx_ta_n, :)];
    idx_te_p = randsample(find(Y_te == 1), n_trial);
    idx_te_n = randsample(find(Y_te == -1), n_trial);
    rX_te = [X_te(idx_te_p, :); X_te(idx_te_n, :)];
    rY_te = [Y_te(idx_te_p, :); Y_te(idx_te_n, :)];

    % 10-fold CV partition
    nfold = 10;
    cv = cvpartition(rY_ta, 'KFold', nfold);

    for csample = 1:nfold

        X_train = rX_ta(cv.training(csample), :);
        Y_train = rY_ta(cv.training(csample));
        X_test = rX_te(cv.test(csample), :);
        Y_test = rY_te(cv.test(csample));

        % train SVM (matching searchlight method)
        SVMModel = fitcsvm(X_train, Y_train, 'ClassNames', [-1, 1]);
        [Y_pred, score] = predict(SVMModel, X_test);
        errTable_te = confusionmat(Y_test, Y_pred, 'Order', [-1, 1]);

        % accuracy
        acc(csample, 1) = calc_percor(errTable_te) ./ 100;
        % AUC
        [~, ~, ~, auc(csample, 1), ~] = perfcurve(Y_test, score(:, 2), max(Y_te));

        % sensitivity/specificity (errTable_te, pos/row1=L, neg/row2=H)
        sen(csample, 1) = errTable_te(2, 2) / sum(errTable_te(2, :));
        spe(csample, 1) = errTable_te(1, 1) / sum(errTable_te(1, :));

        % number of voxels
        nv(csample, 1) = size(X_train, 2);
 
    end % end csample

end % end function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [acc, auc, sen, spe, nv] = SVM_leaveRun(DATA_ta, DATA_te)
    % leave one session out CV with SVM

    X_ta = DATA_ta.X;
    Y_ta = DATA_ta.Y;
    sess_ta = DATA_ta.session;
    X_te = DATA_te.X;
    Y_te = DATA_te.Y;
    sess_te = DATA_te.session;
    nrun = max(sess_ta);

    for csample = 1:nrun

        X_train = X_ta(sess_ta ~= csample, :);
        Y_train = Y_ta(sess_ta ~= csample);
        X_test = X_te(sess_te == csample, :);
        Y_test = Y_te(sess_te == csample);

        % train svm
        SVMModel = fitcsvm(X_train, Y_train, 'ClassNames', [-1, 1]);
        [Y_pred, score] = predict(SVMModel, X_test);
        errTable_te = confusionmat(Y_test, Y_pred, 'Order', [-1, 1]);

        % accuracy
        acc(csample, 1) = calc_percor(errTable_te) ./ 100;
        % AUC
        [~, ~, ~, auc(csample, 1), ~] = perfcurve(Y_test, score(:, 2), max(Y_te));

        % sensitivity/specificity (errTable_te, pos/row1=L, neg/row2=H)
        sen(csample, 1) = errTable_te(2, 2) / sum(errTable_te(2, :));
        spe(csample, 1) = errTable_te(1, 1) / sum(errTable_te(1, :));
        
        % number of voxels
        nv(csample, 1) = size(X_train, 2);
    
    end % end csample

end % end function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DATA = balance_data(DAT)
    % balancing high/low pain trials in each session (run)

    % initialise
    feat_all = [];
    label_all = [];
    sess_all = [];

    ntrial = size(DAT.label);
    nsess = ntrial(2);

    for csess = 1:nsess

        dum_feat = DAT.feature{csess};
        dum_label = DAT.label(:, csess);
        dum_sess = repmat(csess, length(dum_label), 1);
        
        % balancing number of samples
        h_num = sum(dum_label == 1);
        l_num = sum(dum_label == -1);
        min_num = min([h_num, l_num]);
        
        % randomly choose samples
        h_idx = find(dum_label == 1);
        l_idx = find(dum_label == -1);
        h_chosen = datasample(h_idx, min_num);
        l_chosen = datasample(l_idx, min_num);
        chosen_idx = [h_chosen; l_chosen];
        chosen_idx = chosen_idx(randperm(length(chosen_idx)));

        % chosen trials in session
        chosen_feat = dum_feat(chosen_idx, :);
        chosen_label = dum_label(chosen_idx, :);
        chosen_sess = dum_sess(chosen_idx, :);

        % concat all session
        feat_all = [feat_all; chosen_feat];
        label_all = [label_all; chosen_label];
        sess_all = [sess_all; chosen_sess];
    
    end % end csess

    DATA.X = feat_all;
    DATA.Y = label_all;
    DATA.session = sess_all;

end % end function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [acc, auc, sen, spe, nv] = SVM_param_balance_tate(PAR)
    % main function for CV with SVM

    % setting parameters and load data files
    DAT_ta = load([PAR.tadir, num2str(PAR.roiuse, '%03d'), filesep, 'feature_', PAR.roiName]); % load 'feature', 'label'
    DAT_te = load([PAR.tedir, num2str(PAR.roiuse, '%03d'), filesep, 'feature_', PAR.roiName]); % load 'feature', 'label'

    % balance data first
    DATA_ta = balance_data(DAT_ta);
    DATA_te = balance_data(DAT_te);

    % perform leave run out CV
    [acc, auc, sen, spe, nv] = SVM_leaveRun(DATA_ta, DATA_te);

    % alternatively, run 10 fold CV with balance data
    % [acc,auc,sen,spe,nv] = SVM_10folds(DATA_ta,DATA_te);
end
