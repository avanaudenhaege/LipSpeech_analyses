% (C) Copyright 2020 CPP ROI developers

%% Script function 
% create ROIs from an individual peak coordinate
% within the localizer mask (following boundaries)
% and extract func data

%This script works with input :
%     - Individual Peak coordinates (from script get_peak_coordinates)
%     - dataImage = image of reference, one of the 3D beta/tmap volume of your
%     decoding task
%     - LocMask = the localizer binary mask that was saved manually in spm,
%     from the localizer data (at a liberal treshold p<0.001unc, in order to have
%     enough voxels). 


clear;
clc;

%% ASSUMPTION
%
% This assumes that the 2 images are in the same space (MNI, individual)
% but they might not necessarily have the same resolution.
%
% In SPM lingo this means they are coregistered but not necessarily resliced.
%

%% IMPORTANT: for saving ROIs
%
% If you want to save the ROI you are creating, you must make sure that the ROI
% image you are using DOES have the same resolution as the image you will
% sample.
%
% You can use the resliceRoiImages for that.
% 

%% Define dir&paths
this_dir = fileparts(mfilename('fullpath'));
mvpa_dir = fullfile(this_dir, '..');
lipspeech_dir = fullfile(mvpa_dir, '..');
bids_dir = fullfile(lipspeech_dir, 'Lip_BIDS');
output_dir = fullfile(this_dir, 'expand-withinmask');

%% Define subjects and ROIs to work on
subject_label = {'04', '05', '06', '07', '08', '09', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24', '26', '27'}; 
roi_label = {'ppaR'};%, 'ppaL', 'phonoR', 'phonoL', 'vwfa', 'ffaR', 

val = 'beta'; %or 'tmap' ?

%% script

for r=1:length(roi_label)
    roi = roi_label{r};
    load(strcat(roi, '-IndividualPeaks.mat')); %mat file with individual peaks of VWFA
    
    opt.save.roi = true;
    if opt.save.roi
      opt.reslice.do = true;
    else
      opt.reslice.do = false;
    end

    for s=1:length(subject_label)
        sub = subject_label{s};
        sub_name = strcat('sub-',sub);
        sub_dir = fullfile(output_dir, roi, sub_name);
        
        if ~exist(sub_dir, 'dir')
            mkdir (sub_dir)
        end 
        opt.outputDir = fullfile(sub_dir); % if this is empty new masks are saved in the current directory.
                         
                if strcmp(val, 'beta')
                    file = 'con';
                elseif strcmp(val, 'tmap')
                    file = 'spmT';
                end 

                dataImage= fullfile(bids_dir, 'derivatives', 'bidspm-stats', sub_name, ...
                    '/task-MVPAAud_space-IXI549Space_FWHM-2_node-MVPAAudCons', ...
                    strcat('/', file, '_0001.nii'));
%% open if I want to make a ROI > 300vx (PPA NeuroSynth not enough vx)                 
% % %                 if strncmp(roi, 'ppa', 3)
% % %                     %for ppa, not have vx in NeuroSynth Mask (around 200vx)
% % %                     %so I use the group localizer mask
% % %                     P = fullfile(mvpa_dir,'ROIs', 'group-ROIs');
% % %                     F = sprintf('%s*', roi);
% % %                     D = dir(fullfile(P,F)); 
% % %                     N = D.name; 
% % %                 else
% % %                     % Look for the good Neurosynth mask to grow the sphere in it. 
% % %                     P = fullfile(mvpa_dir,'ROIs', 'NeuroSynth');
% % %                     F = sprintf('%s*', roi);
% % %                     D = dir(fullfile(P,F)); 
% % %                     N = D.name; 
% % %                 end
%%                 
                % Look for the good Neurosynth mask to grow the sphere in it. 
                P = fullfile(mvpa_dir,'ROIs', 'NeuroSynth');
                F = sprintf('%s*', roi);
                D = dir(fullfile(P,F)); 
                N = D.name; 
                %Get the path and name of this mask to use in the function getDataFromExpansion
                LocMask = fullfile(P, N);

                opt.location = Coord(s,:);
                opt.roiname = roi;
                

                data_expand = getDataFromExpansion(opt,dataImage,LocMask);

    end
end 


function data_expand = getDataFromExpansion(opt, dataImage, LocMask)

%   % X Y Z coordinates of left VWFA
%   %If I want to use canonical coordinates instead of individual peak !
%   location = [-45 -57 -12];

%   % radius in millimeters
%   radius = 8;
% 
%   sphere.location = location;
%   sphere.radius = radius;
  
  sphere = struct; 
  sphere.location = opt.location; %for individual peak 
  sphere.radius = 1;
  sphere.maxNbVoxels = 200;

  specification = struct('mask1', LocMask, 'mask2', sphere); 

  mask = createRoi('expand', specification, dataImage, opt.outputDir, opt.save.roi);
  
  lookfor = dir(fullfile(opt.outputDir, ['label-', mask.label, '_', opt.roiname, '*', '.nii']));
  maskpath = fullfile(lookfor.folder, lookfor.name);  
  
 %%%%%%%%%%% à retirer si les deux lignes au dessus fonctionnent
% % % % % %   if strcmp(opt.roiname, 'vwfa')
% % % % % %       lookfor = dir(fullfile(opt.outputDir, ['label-', mask.label, '_vwfa', '*', '.nii']));
% % % % % %       maskpath = fullfile(lookfor.folder, lookfor.name);
% % % % % %       %maskpath = fullfile(opt.outputDir, strcat('label-', mask.label, '_rvwfa-neurosynthCluster_mask.nii'));
% % % % % %   elseif strncmp(opt.roiname, 'ffa', 3)
% % % % % %       maskpath = fullfile(opt.outputDir, strcat('label-', mask.label, '_all-cluster_faces-gtothers_pt001-unc_mask.nii'));
% % % % % %   elseif strncmp(opt.roiname, 'ppa', 3)
% % % % % %       maskpath = fullfile(opt.outputDir, strcat('label-', mask.label, '_all-cluster_houses-gtothers_pt001-unc_mask.nii'));
% % % % % %   elseif strncmp(opt.roiname, 'phono', 5)
% % % % % %       maskpath = fullfile(opt.outputDir, strcat('label-', mask.label, '_all-cluster_syll-gtscr_pt001-unc_mask.nii'));
% % % % % %   end
% % % % % %   
    disp(opt.roiname)
    disp(mask.label)
  
    
    
  %%%% PROBLEM : FIND HOW TO MODIFY NAME OF OUTPUT !!!! %%%%
  % such as newname of the mask = fullfile(opt.outputDir, strcat(opt.roiname, '_label-', mask.label, '_mask.nii'))
  %%%% bcs here it gives an output but I don't know if it is the left or
  %%%% the right roi and it does not match with the spherical rois
  
  
  %need to reslice the mask, because after createRoi, it is saved outside
  %of space a brain space, just a mask of [x,y,z] voxels.
  mask = resliceRoiImages(dataImage,maskpath);

  data_expand = spm_summarise(dataImage, mask);
  
  
  % equivalent to
  % b = spm_summarise(dataImage, ...
  %                   struct( ...
  %                          'def', 'sphere', ...
  %                          'spec', radius, ...
  %                          'xyz', location'));

end


% 
% 
% 
% 
% % ROI: use the probability map for visual motion from Neurosynth
% %   link: https://neurosynth.org/analyses/terms/visual%20motion/
% %
% % Data: tmap of Listening > Baseline from the MoAE demo
% 
% 
% 
% run ../../initCppRoi;
% 
% 
% 
% %%
% zMap = fullfile(pwd, 'inputs', 'visual motion_association-test_z_FDR_0.01.nii');
% dataImage = fullfile(pwd, 'inputs', 'TStatistic.nii');
% 
% opt.unzip.do = true;
% opt.save.roi = true;
% opt.outputDir = []; % if this is empty new masks are saved in the current directory.
% if opt.save.roi
%   opt.reslice.do = true;
% else
%   opt.reslice.do = false;
% end
% 
% % all of these functions can be found below and show you how to create ROIs and
% % / or ROIs to extract data from an image.
% %
% [roiName, zMap] = prepareDataAndROI(opt, dataImage, zMap);
% 
% data_mask = getDataFromMask(dataImage,  roiName);
% data_sphere = getDataFromSphere(opt, dataImage);
% data_intersection = getDataFromIntersection(opt, dataImage,  roiName);
% data_expand = getDataFromExpansion(opt, dataImage,  roiName);
% 
% %% Mini functions
% 
% % only to show how each case works
% 
% function data_mask = getDataFromMask(dataImage, roiName)
% 
%   data_mask = spm_summarise(dataImage, roiName);
% 
% end
% 
% 
% 
% function data_intersection = getDataFromIntersection(opt, dataImage,  roiName)
%   %
%   % Gets the voxels at the intersection of:
%   % - a binary mask and user defined sphere
%   % - TODO: 2 binary masks
%   %
% 
%   % X Y Z coordinates of right V5 in millimeters
%   location = [44 -67 0];
% 
%   sphere.location = location;
%   sphere.radius = 5;
% 
%   specification  = struct( ...
%                           'mask1', roiName, ...
%                           'mask2', sphere);
% 
%   mask = createRoi('intersection', specification, dataImage, opt.outputDir, opt.save.roi);
% 
%   data_intersection = spm_summarise(dataImage, mask.roi.XYZmm);
% 
% end
% 
% function data_expand = getDataFromExpansion(opt, dataImage,  roiName)
%   %
%   % will expand a ROI with a series of expanding spheres but within the
%   % constrains of a binary mask
%   %
%   % the expansion stops once the number of voxels goes above a user defined
%   % threshold.
%   %
% 
%   % X Y Z coordinates of right V5 in millimeters
%   location = [44 -67 0];
% 
%   sphere.location = location;
%   sphere.radius = 1; % starting radius
%   sphere.maxNbVoxels = 50;
% 
%   specification  = struct('mask1', roiName, ...
%                           'mask2', sphere);
% 
%   mask = createRoi('expand', specification, dataImage, opt.outputDir, opt.save.roi);
% 
%   data_expand = spm_summarise(dataImage, mask.roi.XYZmm);
% 
% end
% 
% %% HELPER FUNCTION
% 
% function [roiName, zMap] = prepareDataAndROI(opt, dataImage, zMap)
% 
%   if opt.unzip.do
%     gunzip(fullfile('inputs', '*.gz'));
%   end
% 
%   % give the neurosynth map a name that is more bids friendly
%   %
%   % space-MNI_label-neurosynthKeyWordsUsed_probseg.nii
%   %
%   zMap = renameNeuroSynth(zMap);
% 
%   if opt.reslice.do
%     % If needed reslice probability map to have same resolution as the data image
%     %
%     % resliceImg won't do anything if the 2 images have the same resolution
%     %
%     % if you read the data with spm_summarise,
%     % then the 2 images do not need the same resolution.
%     zMap = resliceRoiImages(dataImage, zMap);
%   end
% 
%   % Threshold probability map into a binary mask
%   % to keep only values above a certain threshold
%   threshold = 10;
%   roiName = thresholdToMask(zMap, threshold);
%   roiName = removePrefix(roiName, spm_get_defaults('realign.write.prefix'));
% 
% end
