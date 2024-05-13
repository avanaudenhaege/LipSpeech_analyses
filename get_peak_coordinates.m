% % % (C) Copyright 2020 CPP ROI developers
% % %Below is the original code from the demo within CPP_Roi
% % % Shows how to get the peak coordinate within a ROI
% % roiImage = extractRoiFromAtlas(pwd, 'wang', 'V1v', 'L');
% % 
% % % Data to read the maximum from
% % gunzip(fullfile(pwd, 'inputs', '*.gz'));
% % dataImage = fullfile(pwd, 'inputs', 'TStatistic.nii');
% % 
% % % If there is no value above a certain threshold the function will return NaN
% % threshold = 1;
% % 
% % % The image and the ROI must have the same dimension if we want to use the threshold option
% % reslicedImages = resliceRoiImages(dataImage, roiImage);
% % 
% % % Get to work.
% % [worldCoord, voxelCoord, maxVal] = getPeakCoordinates(dataImage, reslicedImages, threshold);

%%spmT number in my visual localizer
% 013 {T} : face_gt_others
% 014 {T} : word_gt_others
% 015 {T} : houses_gt_others

%%spmT number in my phono localizer
% 012 {T} : Syl_gt_scr
% 013 {T} : Scr_gt_syl

%%spmT number in my TVSA localizer
% 011 {T} : VL_gt_NL


%% Define dir&paths
this_dir = fileparts(mfilename('fullpath'));
mvpa_dir = fullfile(this_dir, '..');
lipspeech_dir = fullfile(mvpa_dir, '..');
bids_dir = fullfile(lipspeech_dir, 'Lip_BIDS');

%% Define subjects and ROIs to work on
subject_label = {'04', '05', '06', '07', '08', '09', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24', '26', '27'}; 
roi_label = {'ppaR'};%, 'ppaL', 'vwfa', 'ffaR', 'phonoL', 'phonoR', 'TVSA'};%, 'vwfa', 'ffaR', 'ffaL', 'ppaR', 'ppaL', 'phonoL', 'phonoR'

%%!!!!!!!!!!! Don't forget to load bidspm before starting :) 

%% Get the peak
for r=1:length(roi_label)
    roi = roi_label{r};
    
    % Look for the good Group ROI to look for the peak in it. 
    P = fullfile(this_dir,'group-ROIs');
    F = sprintf('%s*', roi);
    D = dir(fullfile(P,F)); 
    N = D.name; 
    %Get the path and name of this GroupROI to use in the function GetPeakCoordinates
    roiImage = fullfile(P, N);
    

    for s=1:length(subject_label)
        sub = subject_label{s};
        %this is the path to the spmT contrast I want to use to build my ROI (e.g.
        %for VWFA the contrast words>otherCategories
        if strcmp(roi, 'vwfa')
            dataImage= fullfile(bids_dir, '/derivatives/bidspm-stats', ...
            strcat('sub-', sub, '/task-VisLoc_space-IXI549Space_FWHM-6/spmT_0015.nii'));%spmT_0013 pour sub-04 - changed manually in the file
        elseif strncmp(roi, 'ffa', 3)
            dataImage= fullfile(bids_dir, '/derivatives/bidspm-stats', ...
            strcat('sub-', sub, '/task-VisLoc_space-IXI549Space_FWHM-6/spmT_0014.nii'));%spmT_0012 ''
        elseif strncmp(roi, 'ppa', 3)
            dataImage= fullfile(bids_dir, '/derivatives/bidspm-stats', ...
            strcat('sub-', sub, '/task-VisLoc_space-IXI549Space_FWHM-6/spmT_0016.nii'));%spmT_0014 ''
        elseif strncmp(roi, 'phono', 4)
            dataImage= fullfile(bids_dir, '/derivatives/bidspm-stats', ...  
            strcat('sub-', sub, '/task-PhonoLoc_space-IXI549Space_FWHM-6/spmT_0012.nii'));
        elseif strncmp(roi, 'TVSA', 4)
            %for TVSA, I have a few sub that did not run this localizer.
            dataImgpath = strcat(bids_dir, '/derivatives/bidspm-stats/sub-', sub, '/task-TVSALoc_space-IXI549Space_FWHM-6/spmT_0011.nii');
            %Check if the subject did TVSA localizer
            if exist(dataImgpath, 'file')
                dataImage = fullfile(dataImgpath); 
            %if yes, take their image (spmT_0011.nii)  
            %For the one who did not, use the group peak coordinate [-57 -31 3] from the group spmT file. 
            else 
                dataImage = fullfile(bids_dir, '/derivatives/bidspm-stats/derivatives/bidspm-groupStats/task-TVSALoc_Stats-GroupLevel_sub-ALL/spmT_0001.nii');
            end
            
            %dataImage= fullfile(bids_dir, '/derivatives/bidspm-stats', ...
            %strcat('sub-', sub, '/task-TVSALoc_space-IXI549Space_FWHM-6/spmT_0011.nii'));

                
        end 
        % If there is no value above a certain threshold the function will return NaN
        threshold = 2;

         % The image and the ROI must have the same dimension if we want to use the threshold option
        reslicedImage = resliceRoiImages(dataImage, roiImage);

        % Get to work.
        [worldCoord, voxelCoord, maxVal] = getPeakCoordinates(dataImage,reslicedImage, threshold);
        
        subID(s,1)={sub};
        Coord(s,:)=worldCoord;
        maxValues(s,1)=maxVal;
        
        output_name = strcat(roi, '-IndividualPeaksNEW'); 
        save(output_name,'subID', 'Coord','maxValues'); %, '-append'); % --> this is if you want to add new data to existing file
    end
end 