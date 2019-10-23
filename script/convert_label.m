% convert proprietary experiment result files from .txt to .mat (extracting label)

clear
addpath('../../coadaptask/coadapt_experiment/experiment/')
sjList = 1:19;

for d = 1:2

    for sj = 1:length(sjList)

        subject = ['adpt', num2str(sjList(sj), '%02d')]; % input txt data
        subject_rt = [num2str(sjList(sj), '%02d')];
        coadaptaskDir = '../../coadaptask';

        TOPDIR = fullfile('..', 'dat', subject_rt);

        if d == 1
            % day 1 (train)
            txt_dir = fullfile(coadaptaskDir, 'DATA', subject, 'output_sham');
            SAVEDIR = fullfile(TOPDIR, 'params', 'scan_idx_d1');
        elseif d == 2
            % day 2 (rt)
            txt_dir = fullfile(coadaptaskDir,'DATA',subject,'output_adpt');
            SAVEDIR = fullfile(TOPDIR,'params','scan_idx_d2');
        end
        if (~exist(SAVEDIR)), mkdir(SAVEDIR); end

        % sbj information
        behavFiles = dir(fullfile(txt_dir, '*.txt'));
        totalSessions = length(behavFiles);
        totalScans = 188 - 3; % in each session

        numFirstBlankTR = 3; %dummy scan
        numBaselineTR = 5; %dummy scan
        % rest=2, test=3 TR
        numBoldDelayTR = 2; % pretend bold delay and move scans by
        numTestTR = 3; % how many TR to use for decoder construction
        numITITR = 1; % how many TR for ITI
        numTrialTR = numBoldDelayTR + numTestTR + numITITR; % TR in one trial
        TR = 2; %sec

        % initialise scan info
        scan_idx = zeros(totalScans, totalSessions);
        trial_label = [];

        % converting for SLR training
        for ii = 1:totalSessions

            % load task conditions from behavioural files
            % (this is a proprietary function that reads ascii txt files produced in the experiment, not included in this repo)
            stim_data = load_ascii_data(txt_dir, behavFiles(ii).name);

            % getting shock time and conditions
            shockLevel = stim_data.data.painlevel;
            numTrials = length(shockLevel);
            fprintf('session %d, high perc=%.3f, #high shock=%d, #low shock=%d.\n', ii, sum(shockLevel == 1) ./ length(shockLevel), sum(shockLevel == 1), sum(shockLevel == -1))

            % baseline 5 scans labelled as 99
            baseline_start = 1;
            baseline_end = numBaselineTR;
            scan_idx(baseline_start:baseline_end, ii) = 99 * ones(numBaselineTR, 1);

            % produce scan index and trial labels
            for jj = 1:numTrials
                trialStart = numBaselineTR + numBoldDelayTR + numTrialTR * (jj - 1) + 1; %trial start scan
                trialEnd = trialStart + numTestTR - 1; % trial end scan
                scan_idx(trialStart:trialEnd, ii) = jj * ones(numTestTR, 1);
            end

            trial_label(:, ii) = shockLevel; % high shock = 1, low shock = -1

        end % end ii

        % saving label files
        saveName = ['scan_idx_', 'rest', num2str(numBoldDelayTR), 'test', num2str(numTestTR), '.mat'];
        save(fullfile(SAVEDIR, saveName), 'trial_label', 'scan_idx')

    end % end sj

end % end d
