function [] = run_make_feature_01(sj)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Feature extraction
% Functions to extract features for decoder training
% - extract temporal signals from ROI masks
%   - progressive detrend (as in real-time)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % subject list
    % sjList = 1:19;
    % sjList = 1;

    % for ii = 1:length(sjList)
    for ii = sj

        for dd = 1:2 % days

            % file directory
            % PAR.name   = [num2str(sjList(ii),'%02d')];
            PAR.name   = [num2str(ii,'%02d')];
            PAR.day = num2str(dd);
            PAR.TOPDIR = fileparts(pwd);

            PAR.nscan  = 188-3; % total already excluded dummy
            PAR.scanidx_suff = '_rest2test3'; % scan idx file suffix (restNtestN)
            PAR.DATDIR = [PAR.TOPDIR,'/dat/'];
            PAR.SUBDIR = [PAR.DATDIR,PAR.name,filesep];
            PAR.EPIDIR = [PAR.SUBDIR,'/EPI/'];
            PAR.T1DIR  = [PAR.SUBDIR,'/t1/'];
            PAR.ROIDIR = [PAR.SUBDIR,'/ROI/d',PAR.day,'/'];
            PAR.PARADIR= [PAR.SUBDIR,'/params/'];
            PAR.FEATDIR= [PAR.SUBDIR,'/features_scantrend/'];
            
            % progressive detrend (as in real-time)
            make_roi_signal_scantrend(PAR);

            if ( ~exist( PAR.PARADIR ) ), mkdir( PAR.PARADIR ); end
            if ( ~exist( PAR.FEATDIR ) ), mkdir( PAR.FEATDIR ); end

        end % end dd

        % remove EPI folder to save discspace
        rmdir(PAR.EPIDIR, 's')

    end % end ii

end % end function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = make_roi_signal_scantrend(PAR)
    % extract ROI signal with progressive detrend

    SUBDIR      = PAR.SUBDIR;
    
    for dd = 1:2

        dir_name = ['tmpl_d',PAR.day,'_d',num2str(dd)];
        EPIDIR     = [SUBDIR,'EPI/',dir_name,'/'];
        SAVEDIR     = [PAR.PARADIR,'/roi_signal_scantrend/',dir_name,'/']; 
        scanidx_suff = PAR.scanidx_suff;
    
        % load label info
        if dd == 1
            scan_info = load([PAR.PARADIR,'scan_idx_d1/scan_idx',scanidx_suff,'.mat']); % for d1 data
        elseif dd == 2
            scan_info = load([PAR.PARADIR,'scan_idx_d2/scan_idx',scanidx_suff,'.mat']); % for d2 data
        end

        scan_idx = scan_info.scan_idx;
        trial_label = scan_info.trial_label;
        label = trial_label;
        nsess = size(trial_label,2);
        nscan = PAR.nscan; % all scans in session minus dummy
        if ( ~exist( SAVEDIR ) ), mkdir( SAVEDIR ); end
    
        % list realigned EPI files
        EPILIST    = spm_select('FPlist', EPIDIR, '^r.*nii');

        % list subject space ROIs
        ROIDIR      = [SUBDIR,'ROI/d',PAR.day,'/'];
        ROILIST     = (spm_select('FPlist', ROIDIR, '^rw.*nii'));

        % extracting realigned EPI signals
        EPISIG     = spm_read_vols(spm_vol(EPILIST));
        for cscan = 1:size(EPISIG,4)
            temp    = reshape(EPISIG(:,:,:,cscan),1,[]);
            EPISIG_rs(cscan,:) = temp;
        end
    
        % extracting preproc ROI signals
        ROI_temp    = spm_read_vols(spm_vol(ROILIST));

        for croi = 1:size(ROI_temp,4)

            SAVEDIR_croi        = [SAVEDIR,num2str(croi,'%03d'),filesep]; 
            if ( ~exist( SAVEDIR_croi ) ), mkdir( SAVEDIR_croi ); end
            SUBSAVEDIR_croi = [PAR.FEATDIR,dir_name,'/',num2str(croi,'%03d'),filesep];
            if ( ~exist( SUBSAVEDIR_croi ) ), mkdir( SUBSAVEDIR_croi ); end
    
            dum_roimsk          = ROI_temp(:,:,:,croi) > 0;
            ROISIG_rs           = EPISIG_rs(:,dum_roimsk);
            ROISIG_dt           = [];
            roi_vol             = [];
            feature             = [];
            end_scan = find(diff(scan_idx(:,1))<0);
            end_scan(1) = []; % first end from baseline scans

            for csess = 1:nsess

                if csess == 1
                    session_start       = 1;
                    session_end         = nscan;
                else
                    session_start = session_end + 1;
                    session_end   = session_start + nscan - 1;
                end
                baseline_scans = find(scan_idx(:,csess) == 99); % baseline scans at beginning of session
                end_idx = end_scan + (csess-1)*nscan;
                ctrial = 1;
                
                for sc = session_start:session_end

                    % session_scans = session_start:session_start+sc;
                    roi_vol = ROISIG_rs(session_start:sc,:); % all voxels' signal in ROI
                    roi_vox_num = size(roi_vol,2); % number of voxels in ROI
                    if isempty( find(isnan(roi_vol)) )	% roi_vol nan check
                        mean_vol = mean(roi_vol);
                            ROISIG_dt{csess}(1:size(roi_vol,1),:) =...
                        detrend(roi_vol,'linear') + ones(size(roi_vol,1),1)*mean_vol;
                    else					
                        for ii=1:roi_vox_num
                            p = ~isnan( roi_vol(:,ii) );
                            ROISIG_dt{csess}(p,ii) =...
                        detrend(roi_vol(p,ii),'linear') + mean(roi_vol(p,ii));
                        end
                    end

                    % extract each trial test scans
                    if sc == session_start + 4
                        baseline_data = ROISIG_dt{csess}(baseline_scans,:); 
                        roi_vox_num = size(ROISIG_dt{csess},2);
            
                        % baseline data processing
                        if isempty( find(isnan(baseline_data)) )
                            roi_baseline_mean = mean(baseline_data, 1);
                            roi_baseline_std = std(baseline_data, 0, 1);
                        else
                            roi_baseline_mean = nan(1, roi_vox_num);
                            roi_baseline_std = nan(1, roi_vox_num);
                            for ii  =1:roi_vox_num
                                p = ~isnan( baseline_data(:,ii) );
                                roi_baseline_mean(ii) = mean(baseline_data(p,ii));
                                roi_baseline_std(ii) = std(baseline_data(p,ii), 0);
                            end
                        end
                    
                    elseif length(end_idx)>0 && sc == end_idx(1) && sc > (session_start+5)
                        % only use test scans not immediately after shock
                        test_scans = [sc-session_start-2+1:sc-session_start+1];
                        test_vol = ROISIG_dt{csess}(test_scans,:);
                        test_scan_num = length(test_scans);
            
                        % baseline processing
                        baseline_mean = ones(test_scan_num,1)*roi_baseline_mean;
                        baseline_std = ones(test_scan_num,1)*roi_baseline_std;

                        if isempty( find( roi_baseline_std == 0.0 ) )
                            z_score = (test_vol - baseline_mean)./baseline_std;
                        else
                            z_score = nan(test_scan_num, roi_vox_num);
                            p = roi_baseline_std ~= 0.0;
                            z_score(:,p) = (test_vol(:,p) - baseline_mean(:,p))./baseline_std(:,p);
                        end
                        % check nan
                        if isempty( find(isnan(z_score)) ) &&... 
                            isempty( find(isinf(z_score)) )	
                            z_score_mean = mean(z_score, 1);
                        else
                            z_score_mean = nan(1, roi_vox_num);
                            for ii=1:roi_vox_num
                                p = ~isnan( z_score(:,ii) ) & ~isinf( z_score(:,ii) );
                                z_score_mean(ii) = mean(z_score(p,ii));
                            end
                        end	
            
                    % pulling all trials in session into cell
                    feature{csess}(ctrial,:) = z_score_mean;

                    ctrial = ctrial + 1;
                    end_idx(1) = [];

                    else
                        continue

                    end % end if sc

                end % end for sc

            end % end csess
            
            % saving features with ROI names (last part before extension)
            [a,b,c] = fileparts(ROILIST(croi,:));
            temp = strsplit(b,'_');
            roiName = temp{end}
            save([SUBSAVEDIR_croi,'feature_',roiName],'feature','label','ROILIST')
            clear feature

            % saving raw ROI signal for template construction
            save([SAVEDIR_croi,filesep,'ROISIG_rs'],'ROISIG_rs','ROILIST','-v7.3')

        end % end croi

        % saving a list of ROIs in text
        save([SAVEDIR,filesep,'ROILIST'],'ROILIST','-v7.3')

    end % end dd

end % end function
