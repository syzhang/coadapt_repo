function [] = run_preprocessing_00()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preprocessing
% Functions to prepare dicom files and ROIs
% - import dicom files
% - remove dummy scans
% - copy reference scan
% - realign to reference scan
% - normalise ROI masks to subject space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % subject list
    sjList = 1:19;

    % SPM12 directory
    tempDir = userpath;
    PAR.spmdir = fullfile(tempDir, 'spm12');

    % set number of dummy scans
    PAR.numdum = 3; 

    for ii = 1:length(sjList)

        PAR.name = [num2str(sjList(ii), '%02d')];
        PAR.TOPDIR = fileparts(pwd);

        % defining data directory
        PAR.SUBDIR = fullfile(PAR.TOPDIR, 'dat', PAR.name, filesep);
        PAR.T1DIR = fullfile(PAR.SUBDIR, 't1', filesep);
        PAR.ROIDIR = fullfile(PAR.TOPDIR, 'dat', 'MNI_ROI_to_use', filesep);

        for dd = 1:2 % days
            PAR.day = num2str(dd);

            % preprocessing functions start here
            dicomimport(PAR); % dicom import
            movedummy(PAR); % remove dummy
            create_first(PAR); % copy reference

            % register to 1st EPI reference
            if ( ~exist( fullfile(PAR.T1DIR, 'mmprage.nii') ) )
                run_deformation_1(PAR); % only do this once to save time
            end
            run_deformation_ROI(PAR); % coreg ROIs to 1st EPI

        end % end dd

        for ddd = 1:2 % days
            PAR.day = num2str(ddd);
            run_realign(PAR); % realign
        end % end ddd

    end % end ii

end % end function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = dicomimport(PAR)
    % import dicom files for EPIs and T1
    
    name = PAR.name;
    addpath(genpath(PAR.spmdir))
    spm('Defaults', 'fMRI');
    spm_jobman('initcfg');
    
    % file directory
    TOPDIR = PAR.TOPDIR;
    SUBDIR = PAR.SUBDIR;
    BACKUP_epi = [SUBDIR, filesep, 'd', PAR.day, '/'];
    BACKUP_t1 = [SUBDIR, filesep, 't1/'];
    epi_dir = [SUBDIR, filesep, 'EPI/d', PAR.day, '/'];

    backup_epi_dir = dir(BACKUP_epi);
    backup_epi_dir(1:2) = [];
    t1_dir = dir(BACKUP_t1);
    t1_dir = [BACKUP_t1, t1_dir(3).name];
    j = 1;

    T1DIR = PAR.T1DIR;

    if (~exist(epi_dir)), mkdir(epi_dir); end
    if (~exist(T1DIR)), mkdir(T1DIR); end

    % EPI dicom convert
    if length(backup_epi_dir) < 20 % unlikely to have more than 20 sessions
        % finding dcm files in subdirectories (when collected in sessions for sham)
        for cdir = 1:length(backup_epi_dir)
            dum_dir = [BACKUP_epi backup_epi_dir(cdir).name];
            dcm_files_epi = spm_select('FPList', dum_dir, '^*.dcm$')

            matlabbatch{j}.spm.util.dicom.data = cellstr(dcm_files_epi);
            matlabbatch{j}.spm.util.dicom.root = 'flat';
            matlabbatch{j}.spm.util.dicom.outdir = cellstr(epi_dir);
            matlabbatch{j}.spm.util.dicom.convopts.format = 'nii';
            matlabbatch{j}.spm.util.dicom.convopts.icedims = 0;
            j = j + 1;
        end

    else
        % finding dcm files
        dcm_files_epi = spm_select('FPList', BACKUP_epi, '^*.dcm$');

        matlabbatch{j}.spm.util.dicom.data = cellstr(dcm_files_epi);
        matlabbatch{j}.spm.util.dicom.root = 'flat';
        matlabbatch{j}.spm.util.dicom.outdir = cellstr(epi_dir);
        matlabbatch{j}.spm.util.dicom.convopts.format = 'nii';
        matlabbatch{j}.spm.util.dicom.convopts.icedims = 0;
        j = j + 1;
    end

    % T1 dicom convert
    % finding t1 dcm files
    dcm_files_t1 = spm_select('FPList', t1_dir, '^*.dcm$')

    matlabbatch{j}.spm.util.dicom.data = cellstr(dcm_files_t1);
    matlabbatch{j}.spm.util.dicom.root = 'flat';
    matlabbatch{j}.spm.util.dicom.outdir = cellstr(T1DIR);
    matlabbatch{j}.spm.util.dicom.convopts.format = 'nii';
    matlabbatch{j}.spm.util.dicom.convopts.icedims = 0;
    j = j + 1;

    spm_jobman('run', matlabbatch);

    if (~exist(fullfile(PAR.T1DIR, 'mprage.nii')))
        t1filename = dir([T1DIR, filesep, '*.nii']);
        movefile([T1DIR, filesep, t1filename(1).name], [T1DIR, filesep, 'mprage.nii'])
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = movedummy(PAR)
    % moving dummy scans from EPI dir

    % file directory
    SUBDIR = PAR.SUBDIR;
    ndum = PAR.numdum;
    EPIDIR = [SUBDIR, filesep, 'EPI/d', PAR.day, '/'];
    DUMDIR = [SUBDIR, filesep, 'dummy/d', PAR.day, '/'];
    if (~exist(DUMDIR)), mkdir(DUMDIR); end

    % moving front dummy
    for cdum = 1:ndum
        filename = [EPIDIR, sprintf('*%06d-01.nii', cdum)];

        if ~isempty(dir(filename))
            movefile(filename, DUMDIR)
        end

    end

    dumScan = [189:190]; % number of extra dummy
    % moving end dummy
    for cdum = 1:length(dumScan)
        filename = [EPIDIR, sprintf('*%06d-01.nii', dumScan(cdum))];

        if ~isempty(dir(filename))
            movefile(filename, DUMDIR)
        end

    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = create_first(PAR)
    % create 1stEPI and 1stNii by copying from data dir (after removing dummy)

    % file directory
    SUBDIR = PAR.SUBDIR;
    EPIDIR = [SUBDIR, filesep, 'd', PAR.day, '/'];
    FEPIDIR = [SUBDIR, filesep, '1stEPI/d', PAR.day, '/'];
    FNIIDIR = [SUBDIR, filesep, '1stNii/d', PAR.day, '/'];
    if (~exist(FEPIDIR)), mkdir(FEPIDIR); end
    if (~exist(FNIIDIR)), mkdir(FNIIDIR); end

    % copy first EPI
    epi_files = dir(fullfile(EPIDIR));
    fepi_file_name = epi_files(3 + 3).name; % 3 for list to start in matlab (first two are .. and ., then 3 dummy scans)
    copyfile(fullfile(EPIDIR, fepi_file_name), fullfile(FEPIDIR, fepi_file_name));

    % convert 1st EPI to nii
    dicom_file_name = dir(FEPIDIR);
    dicom_file_name(1:2) = [];
    FEPI_file_name = fullfile(FEPIDIR, dicom_file_name.name);
    hdr = spm_dicom_headers(FEPI_file_name);
    opts     = 'all';
    root_dir = 'flat';
    format   = spm_get_defaults('images.format');
    out_dir  = FNIIDIR;	% output dir
    nifti = spm_dicom_convert(hdr, opts, root_dir, format, out_dir);
    templ_nifti_fname = nifti.files{1};

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = run_realign(PAR)
    % realign to first EPI

    SUBDIR = PAR.SUBDIR;
    FNIIDIR = [SUBDIR, filesep, '1stNii/d', PAR.day, '/'];

    for dd = 1:2
        % get all epi in epi folder
        epi = cellstr(spm_select('FPList', [SUBDIR, 'EPI/d', num2str(dd), '/'], '^f.*.nii'));
        % 1st EPI nii (reference)
        templ_nifti_fname = cellstr(spm_select('FPList', FNIIDIR, '^f.*.nii'));

        % reslice parameters
        reslice_flags = spm_get_defaults('realign.write');
        reslice_flags.interp = 1; % interp : the B-spline interpolation method.
        reslice_flags.mask = 1; % 1: mask output images
        reslice_flags.prefix = 'r'; % prefix for resliced image [default:'r']
        reslice_flags.which = 1; % 1: don't reslice the first image.
        reslice_flags.mean = 0; % 0: not write mean image
        reslice_flags.wrap = [0, 0, 0]; %representing wrapping in each of the dimensions.(*)

        % matching realign/reslice with decnef code
        % Realign & Reslice
        for i = 1:length(epi)
            nifti_fname = epi{i};
            TempfNameArray = ...
                strvcat(templ_nifti_fname{1}, nifti_fname);
            spm_realign(TempfNameArray);
            spm_reslice(TempfNameArray, reslice_flags);
        end

        % move data
        outdir = [SUBDIR, 'EPI/tmpl_d', PAR.day, '_d', num2str(dd), '/'];
        if (~exist(outdir)), mkdir(outdir); end
        movefile([SUBDIR, 'EPI/d', num2str(dd), '/r*.nii'], outdir)
        movefile([SUBDIR, '1stNii/d', PAR.day, '/rp*.txt'], outdir)
    
    end % end dd

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = run_deformation_1(PAR)
    % coregistration and segmentation of t1 

    % Coregister t1 to spm template t1
    single_subject_t1 = fullfile(PAR.spmdir, 'canonical', 'single_subj_T1.nii');
    current_subject_t1 = fullfile(PAR.T1DIR, 'mprage.nii');
    matlabbatch{1}.spm.spatial.coreg.estimate.ref = cellstr(single_subject_t1);
    matlabbatch{1}.spm.spatial.coreg.estimate.source = cellstr([current_subject_t1]);

    % segmentation of coregistered t1
    matlabbatch{2}.spm.spatial.preproc.channel.vols = cellstr(current_subject_t1);
    matlabbatch{2}.spm.spatial.preproc.channel.write = [0 1];
    matlabbatch{2}.spm.spatial.preproc.warp.write = [1 1];

    [~, prov] = spm_jobman('run', matlabbatch);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = run_deformation_ROI(PAR)
    % normalise ROI masks to subject space

    clear matlabbatch
    % imcalc on segmentation outpus
    c1_seg = spm_select('FPListRec', PAR.T1DIR, '^c1.*\.nii$'); % GS
    c2_seg = spm_select('FPListRec', PAR.T1DIR, '^c2.*\.nii$'); % WM
    c3_seg = spm_select('FPListRec', PAR.T1DIR, '^c3.*\.nii$'); % CSF
    ms_seg = spm_select('FPListRec', PAR.T1DIR, '^mm.*\.nii$'); % bias corrected
    imcalc_input = {ms_seg; c1_seg; c2_seg; c3_seg};

    % produce MS
    matlabbatch{1}.spm.util.imcalc.input = imcalc_input;
    matlabbatch{1}.spm.util.imcalc.output = fullfile(PAR.T1DIR, 'mprage_brain_spm.img');
    matlabbatch{1}.spm.util.imcalc.expression = 'i1.*((i2+i3+i4)>0.5)';

    % normalise ROIs to subject space
    roi_list = spm_select('FPList', [PAR.ROIDIR], '^*.nii');
    MNI_sub_file = spm_select('FPListRec', PAR.T1DIR, '^iy.*\.nii$');
    matlabbatch{2}.spm.spatial.normalise.write.subj.def = cellstr(MNI_sub_file);
    matlabbatch{2}.spm.spatial.normalise.write.subj.resample = cellstr(roi_list);
    matlabbatch{2}.spm.spatial.normalise.write.woptions.vox = [2 2 2];
    matlabbatch{2}.spm.spatial.normalise.write.woptions.prefix = 'w';

    [~, prov] = spm_jobman('run', matlabbatch);

    clear matlabbatch

    % coreg to 1st EPI not mean from realign
    SUBDIR = PAR.SUBDIR;
    BACKUP_1stnii = [SUBDIR, filesep, '1stNii/d', PAR.day, '/'];
    dum = dir([BACKUP_1stnii, '/*.nii']);
    mean_file = [BACKUP_1stnii, filesep, dum.name];
    roi_list = cellstr(spm_select('FPList', [PAR.ROIDIR], 'w.*.nii'));
    matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {mean_file};
    matlabbatch{1}.spm.spatial.coreg.estwrite.source = {[PAR.T1DIR, '/mprage_brain_spm.img']};
    matlabbatch{1}.spm.spatial.coreg.estwrite.other = roi_list;
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';
    [~, prov] = spm_jobman('run', matlabbatch);

    clear matlabbatch

    % using imcalc to get rid of bounding box of ROI
    rw_roi_warped = cellstr(spm_select('FPListRec', [PAR.ROIDIR], '^rw.*\.nii$'));

    for ii = 1:length(rw_roi_warped)
        [a, b, c] = fileparts(rw_roi_warped{ii});
        matlabbatch{ii}.spm.util.imcalc.input = cellstr(rw_roi_warped{ii});
        matlabbatch{ii}.spm.util.imcalc.output = fullfile([PAR.ROIDIR], [b, '.nii']);
        matlabbatch{ii}.spm.util.imcalc.expression = 'i1>0.9';
    end

    [~, prov] = spm_jobman('run', matlabbatch);

    % move rw ROI files to subject's own ROI dir
    SUBROIDIR = [SUBDIR, filesep, 'ROI/d', PAR.day, '/'];
    if (~exist(SUBROIDIR)), mkdir(SUBROIDIR); end
    movefile([PAR.ROIDIR, 'rw*.nii'], SUBROIDIR)
    delete(fullfile(PAR.ROIDIR, 'w*.nii'))

end
