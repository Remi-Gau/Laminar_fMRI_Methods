function bold_profiles_vol
clc; clear;

StartDir = fullfile(pwd, '..');
cd (StartDir)
addpath(genpath(fullfile(StartDir, 'code', 'subfun')))

NbLayers = 6;

NbVoxPerSubROI = 2000;

CondNames = {...
    'AStimL','AStimR';...
    'VStimL','VStimR';...
    'TStimL','TStimR'};

% --------------------------------------------------------- %
%                            ROIs                           %
% --------------------------------------------------------- %
Mask_Ori.ROI(1) = struct('name', 'V1_L_thres', 'fname', 'SubjName_lcr_V1_Pmap_Ret_thres_10_data.nii');
Mask_Ori.ROI(end+1) = struct('name', 'V1_R_thres', 'fname', 'SubjName_rcr_V1_Pmap_Ret_thres_10_data.nii');


SubLs = dir('sub*');
NbSub = numel(SubLs);


for iSub = 1:NbSub % for each subject
    
    fprintf('\n\nProcessing %s\n', SubLs(iSub).name)
    
    Mask = Mask_Ori;
    
    for iROI =1:length(Mask.ROI)
        Mask.ROI(iROI).fname = strrep(Mask_Ori.ROI(iROI).fname,'SubjName',SubLs(iSub).name);
    end
    
    SubDir = fullfile(StartDir, SubLs(iSub).name);
    RoiFolder = fullfile(SubDir, 'roi', 'vol');
    AnalysisFolder = fullfile(SubDir, 'ffx_nat', 'vol');
    
    SaveDir = fullfile(SubDir, 'results', 'profiles');
    [~,~,~] = mkdir(SaveDir);
    
    % Gets the number of each beta images and the numbers of the beta of
    % interest
    load(fullfile(SubDir, 'ffx_nat','SPM.mat'))
    [BetaOfInterest, BetaNames] = GetBOI(SPM,CondNames);
    
    for i=1:size(BetaNames,1)
        if BetaNames(i,6)==' '
            tmp(i,1:6) = BetaNames(i,7:12);
        else
            tmp(i,1:6) = BetaNames(i,8:13);
        end
    end
    BetaNames = tmp;
    clear SPM
    
    % Defines what condition each line of the feature matrix corresponds to
    CondLines = nan(20,numel(CondNames));
    FilesList = {};
    iRow = 1;
    
    for iCond=1:numel(CondNames)
        
        tmp=BetaNames(BetaOfInterest,1:length(CondNames{iCond}));
        tmp = BetaOfInterest(strcmp(CondNames{iCond}, cellstr(tmp)));
        
        for i=1:length(tmp)
            
            CondLines(i,iCond) = iRow;
            
            FilesList{end+1,1} = fullfile(AnalysisFolder, ...
                sprintf('r%s_beta-%04d.nii', SubLs(iSub).name, tmp(i)));
            
            iRow = iRow+1;
        end
        
        clear tmp
    end
    
    % Loads which runs happened on which day to set up the CVs
    load(fullfile(StartDir, 'RunsPerSes.mat'))
    Idx = ismember({RunPerSes.Subject}, SubLs(iSub).name);
    RunPerSes = RunPerSes(Idx).RunsPerSes;
    CVs = {...
        1:RunPerSes(1), ...
        RunPerSes(1)+1:RunPerSes(1)+RunPerSes(2),...
        RunPerSes(1)+RunPerSes(2)+1:sum(RunPerSes)};
    clear Idx
    
    % Leave one day out
    LODO = [2,3;1,3;1,2];
    
    
    %% Gets global mask from GLM and ROI masks for the data
    fprintf(' Reading masks\n')
    
    if ~exist(fullfile(AnalysisFolder, ['r' SubLs(iSub).name '_GLM_mask.nii']), 'file')
        try
            gunzip(fullfile(AnalysisFolder, ['r' SubLs(iSub).name '_GLM_mask.nii.gz']))
        catch
            error('The GLM mask file %s is missing.', ['r' SubLs(iSub).name '_GLM_mask.nii'])
        end
    end
    Mask.global.hdr = spm_vol(fullfile(AnalysisFolder, ['r' SubLs(iSub).name '_GLM_mask.nii']));
    Mask.global.img = logical(spm_read_vols(Mask.global.hdr));
    
    for i=1:length(Mask.ROI)
        Mask.ROI(i).hdr = spm_vol(fullfile(RoiFolder, Mask.ROI(i).fname));
    end
    
    hdr = cat(1, Mask.ROI.hdr);
    sts = spm_check_orientations([Mask.global.hdr; hdr]);
    if sts ~= 1
        error('Images not in same space!');
    end
    
    clear sts hdr i
    
    % Create mask in XYZ format (both world and voxel coordinates)
    [X, Y, Z] = ind2sub(size(Mask.global.img), find(Mask.global.img));
    Mask.global.XYZ = [X'; Y'; Z']; % XYZ format
    clear X Y Z
    Mask.global.size = size(Mask.global.XYZ, 2);
    Mask.global.XYZmm = Mask.global.hdr.mat(1:3,:) ...
        * [Mask.global.XYZ; ones(1, Mask.global.size)]; % voxel to world transformation
    
    % Combine masks
    xY.def = 'mask';
    for i=1:length(Mask.ROI)
        xY.spec = fullfile(RoiFolder, Mask.ROI(i).fname);
        [xY, Mask.ROI(i).XYZmm, j] = spm_ROI(xY, Mask.global.XYZmm);
        Mask.ROI(i).XYZ = Mask.global.XYZ(:,j);
        Mask.ROI(i).size = size(Mask.ROI(i).XYZ, 2);
        A = spm_read_vols(Mask.ROI(i).hdr);
        A(isnan(A)) = 0;
        Mask.ROI(i).size(2) = sum(logical(A(:)));
    end
    
    
    clear xY j i A
    
    
    %% Gets Layer labels
    fprintf(' Reading layer labels\n')
    
    LayerLabelsFile = dir(fullfile(SubDir, 'anat', ...
        ['sub-*_MP2RAGE_T1map_Layers-' sprintf('%02.0f', NbLayers) '.nii']));
    
    % Unzip the file if necessary
    if ~isempty(LayerLabelsFile)
        LayerLabelsHdr = spm_vol(fullfile(SubDir, 'anat', ...
            LayerLabelsFile.name));
    else
        try
            LayerLabelsFile = dir(fullfile(SubDir, 'anat', ...
                ['sub-*_MP2RAGE_T1map_Layers-' sprintf('%02.0f', NbLayers) '.nii.gz']));
            gunzip(fullfile(SubDir, 'anat', 'cbs', ...
                LayerLabelsFile.name));
            LayerLabelsHdr = spm_vol(fullfile(SubDir, 'anat', ...
                LayerLabelsFile.name(1:end-3)));
        catch
            error(['The layer label file ' LayerLabels 'is missing.'])
        end
    end
    
    sts = spm_check_orientations([Mask.global.hdr; LayerLabelsHdr]);
    if sts ~= 1
        error('Images not in same space!');
    end
    clear sts
    
    for i=1:length(Mask.ROI)
        LayerLabels{i} = spm_get_data(LayerLabelsHdr, Mask.ROI(i).XYZ); %#ok<*AGROW>
    end
    
    
    %% Mask each image by each ROI and create a features set (images x voxel)
    fprintf('\n Get features\n')
    
    for i=1:length(Mask.ROI)
        
        FeatureFile = fullfile(AnalysisFolder, ['BOLD_' Mask.ROI(i).name '_data_l-' num2str(NbLayers) '.mat']);
        
        % Load the feature file if it exists
        if exist(FeatureFile, 'file')
            load(FeatureFile, 'Features', 'MaskSave', 'FilesListSave')
            
            % Make sure that we have the right ROI
            if ~isequal(MaskSave, Mask.ROI(i))
                NeedFeat(i) = true;
            end
            
            % Make sure that the right features were extracted
            if ~isequal(FilesListSave, FilesList)
                NeedFeat(i) = true;
            end
            
            FeaturesAll{i} = Features{1};
            NeedFeat(i) = false;
            
            % Otherwise flag this ROI to feature extraction
        else
            NeedFeat(i) = true;
        end
        
        clear FilesListSave MaskSave Features
        
        if NeedFeat(i)
            error('There is a problem with the data from this file \n%s',...
                FeatureFile)
        end
    end
    
    Features = FeaturesAll;
    
    clear FilesList FeaturesAll
    
    
    %% Averages across blocks and voxels for each ROI
    for iROI=1:numel(Mask.ROI)
        
        fprintf('\n Processing %s \t', Mask.ROI(iROI).name)
        Features_ROI = Features{iROI};
        LayerLabels_ROI = LayerLabels{iROI};
        
        % Only keep voxels that belong to a layer
        Vox2Keep = LayerLabels_ROI>0;
        Features_ROI = Features_ROI(:,Vox2Keep);
        LayerLabels_ROI = LayerLabels_ROI(Vox2Keep);
        
        clear Data_ROI Vox2Keep SubROI DeactVox ActVox MostAct
        
        Data_ROI.info.name = Mask.ROI(iROI).name;
        Data_ROI.info.fname = Mask.ROI(iROI).fname;
        Data_ROI.info.size = Mask.ROI(iROI).size;
        Data_ROI.info.vox_per_layer = tabulate(LayerLabels_ROI);
        
        disp(Data_ROI.info.vox_per_layer(:,2)')
        
        %% Get subsets of voxels
        
        for iCond = 1:numel(CondNames)
            
            Img2Sel = CondLines(:,iCond);
            Img2Sel(isnan(Img2Sel)) = [];
            
            % If we Z-scored based on the mean activation across days
            SubROI(iCond).Zscore(1,:) = zscore(mean(Features_ROI(Img2Sel,:),1));
            
            % If we Z-scored based on the mean activation for each day
            for iDay = 1:3
                if  iDay==3  && iSub==5 && (iCond==1 || iCond==4)
                    tmp = CVs{iDay};
                    tmp(tmp==20) = [];
                    SubROI(iCond).Zscore(iDay+1,:) = zscore(mean(Features_ROI(Img2Sel(tmp),:),1));
                    clear tmp
                else
                    SubROI(iCond).Zscore(iDay+1,:) = zscore(mean(Features_ROI(Img2Sel(CVs{iDay}),:),1));
                end
            end
        end
        
        
        Data_ROI.Act = struc('LayerMean', [], ...
            'LayerMedian', [], ...
            'NbVoxel', [] );
        Data_ROI.Deact = struc('LayerMean', [], ...
            'LayerMedian', [], ...
            'NbVoxel', [] );
        Data_ROI.MostAct = struc('LayerMean', [], ...
            'LayerMedian', [], ...
            'NbVoxel', [] );
        
        
        %% For each Condition
        for iCond = 1:numel(CondNames)

            for iLayer = 1:NbLayers % Averages over voxels of a given layer
                
                for iCond2 = 1:numel(CondNames)
                    
                    Img2Sel = CondLines(:,iCond2);
                    Img2Sel(isnan(Img2Sel)) = [];
                    
                    CVs = {...
                        1:RunPerSes(1), ...
                        RunPerSes(1)+1:RunPerSes(1)+RunPerSes(2),...
                        RunPerSes(1)+RunPerSes(2)+1:sum(RunPerSes)};
                    
                    if  iSub==5 && (iCond2==1 || iCond2==4)
                        CVs{3}(end) = [];
                    end
                    
                    %% Sub-division: activated/deactivated based on all days
                    % Activation
                    ActVox = all([...
                        LayerLabels_ROI==iLayer ;...
                        SubROI(iCond).Zscore(1,:)>0]);
                    
                    Data_ROI.Act.LayerMean(iLayer,iCond,iCond2,1) = ...
                        nanmean(nanmean(Features_ROI(Img2Sel,ActVox),2));
                    Data_ROI.Act.LayerMedian(iLayer,iCond,iCond2,1) = ...
                        nanmean(nanmedian(Features_ROI(Img2Sel,ActVox),2));
                    Data_ROI.Act.NbVoxel(iLayer,iCond,iCond2,1) = sum(ActVox);
                    
                    % Deactivation
                    DeactVox = all([...
                        LayerLabels_ROI==iLayer ;...
                        SubROI(iCond).Zscore(1,:)<0]);

                    Data_ROI.Deact.LayerMean(iLayer,iCond,iCond2,1) = ...
                        nanmean(nanmean(Features_ROI(Img2Sel,DeactVox),2));
                    Data_ROI.Deact.LayerMedian(iLayer,iCond,iCond2,1) = ...
                        nanmean(nanmedian(Features_ROI(Img2Sel,DeactVox),2));
                    Data_ROI.Deact.NbVoxel(iLayer,iCond,iCond2,1) = sum(DeactVox);
                    
                    
                    %% Sub-division: activated/deactivated with "Cross-validation"
                    % "Cross-validation" subdivision based on e.g day 1 and applied
                    % to mean over day 2 and 3
                    
                    for iCV = 1:numel(CVs)
                        % Activation
                        ActVox = all([...
                            LayerLabels_ROI==iLayer ;...
                            SubROI(iCond).Zscore(iCV+1,:)>0]);
                        
                        % Runs to 'generalize' to
                        TestRuns = [ CVs{LODO(iCV,1)} CVs{LODO(iCV,2)} ];
                        
                        Data_ROI.Act.LayerMean(iLayer,iCond,iCond2,iCV+1) = ...
                            nanmean(nanmean(Features_ROI(Img2Sel(TestRuns),ActVox),2));
                        Data_ROI.Act.LayerMedian(iLayer,iCond,iCond2,iCV+1) = ...
                            nanmean(nanmedian(Features_ROI(Img2Sel(TestRuns),ActVox),2));
                        Data_ROI.Act.NbVoxel(iLayer,iCond,iCond2,iCV+1) = sum(ActVox);
                        
                        % Deactivation
                        DeactVox = all([...
                            LayerLabels_ROI==iLayer ;...
                            SubROI(iCond).Zscore(iCV+1,:)<0]);
                        
                        Data_ROI.Deact.LayerMean(iLayer,iCond,iCond2,iCV+1) = ...
                            nanmean(nanmean(Features_ROI(Img2Sel(TestRuns),DeactVox),2));
                        Data_ROI.Deact.LayerMedian(iLayer,iCond,iCond2,iCV+1) = ...
                            nanmean(nanmedian(Features_ROI(Img2Sel(TestRuns),DeactVox),2));
                        Data_ROI.Deact.NbVoxel(iLayer,iCond,iCond2,iCV+1) = sum(DeactVox);
                    end
                    
                    %% Most activated voxels based on all days
                    % Identify voxels in that layer among the most activated voxels
                    MostAct = false(size(SubROI(iCond).Zscore(1,:)));
                    [~,I] = sort(SubROI(iCond).Zscore(1,:));
                    MostAct(I(1:NbVoxPerSubROI)) = true;
                    MostAct = all([...
                        LayerLabels_ROI==iLayer ;...
                        MostAct]);
                    
                    Data_ROI.MostAct.LayerMean(iLayer,iCond,iCond2,1) = ...
                        nanmean(nanmean(Features_ROI(Img2Sel,MostAct),2));
                    Data_ROI.MostAct.LayerMedian(iLayer,iCond,iCond2,1) = ...
                        nanmean(nanmedian(Features_ROI(Img2Sel,MostAct),2));
                    Data_ROI.MostAct.NbVoxel(iLayer,iCond,iCond2,1) = sum(MostAct);
                    
                    
                    %% Most activated voxels with "Cross-validation"
                    % "Cross-validation" subdivision based on day 1 and applied
                    % to mean over day 2 and 3
                    for iCV = 1:numel(CVs)
                        % Identify voxels in that layer among the most activated voxels
                        MostAct = false(size(SubROI(iCond).Zscore(iCV+1,:)));
                        [~,I] = sort(SubROI(iCond).Zscore(iCV+1,:));
                        MostAct(I(1:NbVoxPerSubROI)) = true;
                        MostAct = all([...
                            LayerLabels_ROI==iLayer ;...
                            MostAct]);
                        
                        % Runs to 'generalize' to
                         TestRuns = [ CVs{LODO(iCV,1)} CVs{LODO(iCV,2)} ];
                        
                        Data_ROI.MostAct.LayerMean(iLayer,iCond,iCond2,iCV+1) = ...
                            nanmean(nanmean(Features_ROI(Img2Sel(TestRuns),MostAct),2));
                        Data_ROI.MostAct.LayerMedian(iLayer,iCond,iCond2,iCV+1) = ...
                            nanmean(nanmedian(Features_ROI(Img2Sel(TestRuns),MostAct),2));
                        Data_ROI.MostAct.NbVoxel(iLayer,iCond,iCond2,iCV+1) = sum(MostAct);
                    end
                end
            end
        end
        
        save(fullfile(SaveDir, strcat('Data_', Mask.ROI(iROI).name, '_l-', num2str(NbLayers), '.mat')), 'Data_ROI')
        
        
        
    end % iSubROI=1:numel(SVM(iSVM).ROI)
    
    fprintf('\n')
    
end % for iSub = 1:NbSub



end