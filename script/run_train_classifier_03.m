function [] = run_train_classifier_03()
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Classifier training
    % - train SLR classifiers
    % - prepare files needed to run day 2 experiment
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % subject list
    % sjList = 1:19;
    sjList = 1;

    % add SLR toolbox to path
    PAR.TOPDIR = fileparts(pwd);
    toolboxDir = fullfile(PAR.TOPDIR, 'SLR1.51');
    addpath(toolboxDir);

    numROIs = 1;

    for ii = 1:length(sjList)

        PAR.name = [num2str(sjList(ii), '%02d')];
        PAR.scanidx_suff = '_rest2test3'; %scan idx file suffix (restNtestN)

        for jj = 1:numROIs
            % which roi to use
            PAR.roiuse = jj;

            % feature directory
            FEATUREDIR = [PAR.TOPDIR, '/dat/', PAR.name, filesep, 'features_scantrend/'];
            PAR.FEATDIR = [FEATUREDIR, 'tmpl_d1_d1/']; % day 1 data aligned to day 1 reference
            featFile = dir([PAR.FEATDIR, num2str(PAR.roiuse, '%03d'), filesep, 'feature_*.mat']);
            featFile
            temp = strsplit(featFile.name(1:end - 4), '_');
            PAR.roiName = temp{end};

            % run SLR and saving weights
            SLR_nonCV_balance(PAR);
            make_roi_template(PAR);
            make_roi_weight(PAR);

        end % end jj

        % move reference
        move_dicom(PAR);

    end % end ii

end % end function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = SLR_nonCV_balance(PAR)
    % train SLR with balanced data within each run

    % setting parameters and load data
    TOPDIR = PAR.TOPDIR;
    DATDIR = [TOPDIR, '/dat/'];
    name = PAR.name;
    scanidx_suff = PAR.scanidx_suff; % sacnidx suffix
    PARADIR = [DATDIR, name, '/params/'];
    FEATDIR = PAR.FEATDIR;

    SLRDIR = [TOPDIR, '/SLR/'];
    scan_info = load([PARADIR, 'scan_idx_d1/scan_idx', scanidx_suff, '.mat']);

    useroi = PAR.roiuse; % ROI to use
    trial_label = scan_info.trial_label;
    nsess = size(trial_label, 2);
    SUBSLRDIR = [SLRDIR, name, filesep];
    if ~exist(SUBSLRDIR); mkdir(SUBSLRDIR); end

    % load feature
    featFile = dir([FEATDIR, num2str(useroi, '%03d'), filesep, 'feature_*.mat']);
    DAT = load([FEATDIR, num2str(useroi, '%03d'), filesep, featFile.name]); % load 'feature', 'label'

    % initialise
    feat_all = [];
    label_all = [];
    sess_all = [];

    % balance data
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

    % train SLR decoder
    [ww_f, ~, ~, errTable_te] = biclsfy_slrvar(DATA.X, DATA.Y, DATA.X, DATA.Y, ...
        'nlearn', 300, 'mean_mode', 'none', 'scale_mode', 'none', 'invhessian', 0, 'usebias', 1);
    slr_results = sum(diag(errTable_te));
    W_data = ww_f;
    performance = sum(slr_results) / size(DATA.X, 1);

    % check calc_label
    correct_count = 0;

    for jj = 1:size(DATA.X, 1)
        [a, label_calc] = calc_label([DATA.X(jj, :), 1], ww_f);
        pred_out = 0;

        if label_calc > 0.5
            pred_out = 1;
        else
            pred_out = -1;
        end

        correct_count = correct_count + (DATA.Y(jj) == pred_out);
    end

    acc2 = correct_count / jj;
    fprintf('\nCalc_label == pred_label perc: %.3f.\n', acc2)

    % saving weights with roi name
    save([SUBSLRDIR ['Selected_W_nonCV_', PAR.roiName]], 'W_data');

end % end function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = make_roi_template(PAR)
    % make ROI mask template for coordinate extraction later

    % parameters
    TOPDIR = PAR.TOPDIR;
    DATDIR = [TOPDIR, '/dat/'];
    name = PAR.name;
    useroi = PAR.roiuse; % ROI to use
    scanidx_suff = PAR.scanidx_suff; % sacnidx suffix
    PARADIR = [DATDIR, name, '/params/'];
    FEATDIR = PAR.FEATDIR;

    scan_info = load([PARADIR, 'scan_idx_d1/scan_idx', scanidx_suff, '.mat']);
    SAVEDIR = [DATDIR, name, '/templates/'];
    if (~exist(SAVEDIR)), mkdir(SAVEDIR); end

    trial_label = scan_info.trial_label;
    nsess = size(trial_label, 2);

    load([PARADIR, 'roi_signal_scantrend/tmpl_d1_d1/', num2str(useroi, '%03d'), filesep, 'ROISIG_rs.mat']); % load raw ROI signal
    load([FEATDIR, num2str(useroi, '%03d'), filesep, 'feature_', PAR.roiName]); % load 'feature', 'label', 'ROILIST'

    usescan = find(scan_info.scan_idx(:, 1:nsess) ~= 0);
    template_sig = mean(ROISIG_rs(usescan, :), 1)';
    ROISIG = spm_read_vols(spm_vol(ROILIST(useroi, :)));
    FNIIDIR = fullfile(DATDIR, name, '1stNii', 'd1');
    dum = dir(fullfile(FNIIDIR, '*.nii'));
    dummy = spm_vol(fullfile(FNIIDIR, dum(1).name));

    % excluding ROISIG zero voxels (can make classifier training useless)
    if size(template_sig) ~= size(ROISIG(ROISIG ~= 0))
        ROISIG(ROISIG ~= 0) = template_sig(template_sig ~= 0);
    else
        ROISIG(ROISIG ~= 0) = template_sig;
    end

    % saving templates
    dummy.fname = [SAVEDIR, '/template.nii'];
    spm_write_vol(dummy, ROISIG)

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = make_roi_weight(PAR)
    %% Setting parameters and Load data files
    TOPDIR = PAR.TOPDIR;
    DATDIR = [TOPDIR, '/dat/'];
    name = PAR.name;

    SLRDIR = [TOPDIR, '/SLR/'];
    OUTDIR = [TOPDIR, '/res/weight_txt/', name, filesep];
    TEMPSIG = spm_read_vols(spm_vol([DATDIR, name, '/templates/template.nii']));
    load([DATDIR, name, filesep, 'params/roi_signal_scantrend/tmpl_d1_d1/ROILIST'])% load ROILIST

    SUBSLRDIR = [SLRDIR, filesep, name, filesep];
    useroi = PAR.roiuse;
    if (~exist(OUTDIR)), mkdir(OUTDIR); end

    % make weight file
    load([SUBSLRDIR, 'Selected_W_nonCV_', PAR.roiName, '.mat']); % load 'W_data'
    ROISIG = spm_read_vols(spm_vol(ROILIST(useroi, :)));

    % convert to index
    idx_temp = find(W_data(1:end - 1, 1) ~= 0);
    W_data_use = W_data(idx_temp);
    W_data_use(end + 1) = W_data(end, 1);
    mni_idx = find(ROISIG(:, :, 1:end) ~= 0);
    mni_idx_w = mni_idx(idx_temp);
    [i, j, k] = ind2sub(size(ROISIG), mni_idx_w);
    TEMPSIG_w = TEMPSIG(mni_idx_w);
    i(end + 1) = 0; j(end + 1) = 0; k(end + 1) = 0; TEMPSIG_w(end + 1) = 0;
    idx_ROI = [i, j, k, W_data_use, TEMPSIG_w];

    % saving weight txt
    fileID = fopen([OUTDIR, 'ROI_', PAR.roiName, '.txt'], 'w');
    fprintf(fileID, '%d %d %d\n', size(ROISIG, 1), size(ROISIG, 2), size(ROISIG, 3));
    fprintf(fileID, '%d\n', size(idx_ROI, 1) - 1);
    fprintf(fileID, '%d %d %d %2.6f %5.6f\n', idx_ROI');
    fclose(fileID);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = move_dicom(PAR)
    % moving reference dicom for day 2 use

    % setting parameters and load data files
    TOPDIR = PAR.TOPDIR;
    DATDIR = [TOPDIR, '/dat/'];
    name = PAR.name;
    OUTDIR = [TOPDIR, '/res/weight_txt/', name, filesep];
    DCMDIR = [DATDIR, name, '/1stEPI/d1/'];
    DCMfile = dir([DCMDIR, '/*.dcm']);

    % copy dicom over to result dir
    copyfile([DCMDIR, DCMfile.name], [OUTDIR, '001_000004_000001.dcm'])

end
