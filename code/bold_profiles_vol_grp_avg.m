clear
clc


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
Mask.ROI(1) = struct('name', 'V1_L_thres', 'fname', 'SubjName_lcr_V1_Pmap_Ret_thres_10_data.nii');
Mask.ROI(end+1) = struct('name', 'V1_R_thres', 'fname', 'SubjName_rcr_V1_Pmap_Ret_thres_10_data.nii');


SubLs = dir('sub*');
NbSub = numel(SubLs);


for iSub = 1:NbSub % for each subject
    
    fprintf('\n\nProcessing %s', SubLs(iSub).name)
    
    SubDir = fullfile(StartDir, SubLs(iSub).name);
    SaveDir = fullfile(SubDir, 'results', 'profiles');

    for iROI=1:numel(Mask.ROI)
        load(fullfile(SaveDir, strcat('Data_', Mask.ROI(iROI).name, '_l-', num2str(NbLayers), '.mat')), 'Data_ROI')
        
        AllSubjects(iROI).Act.LayerMean(:,:,:,:,iSub) = Data_ROI.Act.LayerMean;
        AllSubjects(iROI).Act.LayerMedian(:,:,:,:,iSub) = Data_ROI.Act.LayerMedian;
        AllSubjects(iROI).Act.NbVoxel(:,:,:,:,iSub) = Data_ROI.Act.NbVoxel;
        
        AllSubjects(iROI).Deact.LayerMean(:,:,:,:,iSub) = Data_ROI.Deact.LayerMean;
        AllSubjects(iROI).Deact.LayerMedian(:,:,:,:,iSub) = Data_ROI.Deact.LayerMedian;
        AllSubjects(iROI).Deact.NbVoxel(:,:,:,:,iSub) = Data_ROI.Deact.NbVoxel;
        
        AllSubjects(iROI).MostAct.LayerMean(:,:,:,:,iSub) = Data_ROI.MostAct.LayerMean;
        AllSubjects(iROI).MostAct.LayerMedian(:,:,:,:,iSub) = Data_ROI.MostAct.LayerMedian;
        AllSubjects(iROI).MostAct.NbVoxel(:,:,:,:,iSub) = Data_ROI.MostAct.NbVoxel;
        
    end
end

%%
close all
for iROI=1%:numel(Mask.ROI)
    
    figure(1)
    iSubplot = 1;
    for iCond1=1:numel(CondNames)
            for iCond2=1:numel(CondNames)
                Data2Plot = squeeze(AllSubjects(iROI).Deact.LayerMean(:,iCond1,iCond2,1,:));
                subplot(6,6,iSubplot)
                hold on
                grid on
                plot(1:6,mean(Data2Plot,2))
%                 errorbar(1:6,mean(Data2Plot,2),nanstd(Data2Plot,[],2))
%                 plot(repmat((1:6)',1,10), Data2Plot, 'color', [.5 .5 .5])
                iSubplot = iSubplot + 1;
                
                axis([0.5 6.5 -5 5])
            end
    end
end
