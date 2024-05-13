
%% Demo: fMRI searchlights LDA classifier
%
%To run MVPA searchlight with LDA/SVM classifier, crossvalidation method
%
% #   For CoSMoMVPA's copyright information and license terms,   #
% #   see the COPYING file distributed with CoSMoMVPA.           #


clear all;


%subject definition (see function in same dir)
global sub
sub =sub_data;
no_sub= [27]; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load data and brain mask to run searchlight within 
%config=cosmo_config();

%% Define dir&paths
mvpa_dir = fileparts(mfilename('fullpath'));
mask_dir = fullfile(mvpa_dir, 'ROIs','Masks');
output_dir = fullfile(mvpa_dir, 'Searchlight_CrossM');

lipspeech_dir = fullfile(mvpa_dir, '..'); %root dir
bids_dir = fullfile(lipspeech_dir, 'Lip_BIDS');
stats_dir = fullfile(bids_dir, 'derivatives', 'bidspm-stats');
preproc_dir = fullfile(bids_dir, 'derivatives', 'bidspm-preproc');

%add libsvm in matlab path
addpath('/Applications/libsvm/matlab');

%% Define subjects, ROIs, tasks to work on

sub_all=sub_data; %call function in same dir named sub_data.m - with all info on each subject
nsub_all = length(sub_all);

algo = 'svm'; %'lda' or 'svm'
val = 'beta'; %or 'beta' ?

%%Define paths 
%%%%%study_path= cd;%one folder back 

%tmap to use as data 
model_list = {'Cons'}; %'Cons', 'Vowels' or 'Speak'
modality_label = {'Aud', 'Vis'}; %order of cond in 4D file --> Auditory - Visual 
decoding_cond = {'both'}; %, 'trainA_testV', 'trainV_testA'}; % 

% % task_list = {'Aud', 'Vis'}; %Vis or Aud
% % model_list = {'Cons'}; %'Cons', 'Vowels' or 'Speak'

for d = 1:length(decoding_cond)
    decodingCondition = decoding_cond{d};
    if strcmp(decodingCondition, 'trainV_testA')
        test = 1;
    elseif strcmp(decodingCondition, 'trainA_testV')
        test = 2;
    elseif strcmp(decodingCondition, 'both')
        test = [];
    end 
    
    for f = 1:length(model_list)
        model = model_list{f};


            % reset citation list
            %cosmo_check_external('-tic');


        for s = no_sub

            %exctract the subname to select the correct data 4D file 
            sub_name=(sub(s).id);
            
            working_on=strcat(sub_name,'_CrossMsearchlight_100vx_', algo, '_', decodingCondition, '_', model, '_', val);
            disp(working_on)

            %%% MASK WHOLE BRAIN 
            mask_fn=fullfile(mvpa_dir, 'ROIs','sub-all_T1-mask', 'mean', 'binary-rmean-22sub.nii');%mean of all anat masks - already binarized and resliced
            % mask_fn=fullfile(mask_dir, 'fullbrainMASK_mni.nii');% whole brain mask excluding cerebellum
            %%% MASK SUB
            %KO !! needs to binarize the mask after reslicing !! don't use
            %following lines of code without changing it. 
% %                 %see if mask image already resliced or not
% %                 lookfor = dir(fullfile(preproc_dir, sub_name, 'ses-01', 'anat', strcat('r', sub_name, '_ses-01_space-IXI549Space_res-r1pt0_desc-brain_mask.nii')));
% %                 if isfile(fullfile(preproc_dir, sub_name, 'ses-01', 'anat', strcat('r', sub_name, '_ses-01_space-IXI549Space_res-r1pt0_desc-brain_mask.nii')))
% %                     %if yes, take this resliced image
% %                     mask_fn = fullfile(lookfor.folder, lookfor.name);
% %                 else 
% %                     %if no, reslice roi image
% %                     imageToCheck = fullfile(preproc_dir, sub_name, 'ses-01', 'anat', strcat(sub_name, '_ses-01_space-IXI549Space_res-r1pt0_desc-brain_mask.nii'));              
% %                     referenceImage = fullfile(mvpa_dir, 'ROIs/Masks/reference_image_for_reslice.nii'); % choose any image from the scanner that will be used with the mask.
% %                     mask_fn = resliceRoiImages(referenceImage, imageToCheck, 0);
% %                 end
            
            %find the data of this subject
                   
            data_img=fullfile(stats_dir, 'mvpaVol', sub_name, ...
                strcat(sub_name, '_task-CrossM_space-IXI549Space_FWHM-2_node-CrossM', model, '_desc-4D_', val, '.nii'));
            
            %% LDA classifier searchlight analysis %%
            % This analysis identified brain regions where the categories can be
            % distinguished using a a Linear Discriminant Analysis (LDA) classifier.

            %prepare the targets - first for one modality (1M) then double
            %this vector bcs I have 2modalities in the same 4D file
            targets1M=repmat(1:3,1,str2num(sub_all(s).Nrun))'; %there are 3 consonants (mean value for the 9 iteration of the cons), and each is repeated in each run. 
            targets1M= sort(targets1M); %comment this line if want to see a "random" decoding - labels will not correspond anymore
            targets = cat(1, targets1M, targets1M);

            chunks1M = repmat(1:str2num(sub_all(s).Nrun), 1, 3)';
            chunks = cat(1, chunks1M, chunks1M);

            modality = repmat(1:2, 1, str2num(sub_all(s).Nrun)*3)'; %There are 2 modalities (Aud AND Vis), and we have to repeat this *3(for each cons) and *Nrun for each run (usually 19 or 20) = usually 114 or 120
            modality = sort(modality);

            ds = cosmo_fmri_dataset(data_img, ... 
                                            'mask', mask_fn,...
                                            'targets',targets, ...
                                            'chunks',chunks); %uses target and chunk structure that we just did
            ds.sa.modality = modality; %add modality in the ds info
                
                
            %set targets
            targets=repmat(1:3,1,str2num(sub(s).Nrun))'; %there are 3 consonants (mean value for the 9 iteration of the cons), and each is repeated in each run. 
            targets= sort(targets); %comment this line if want to see a "random" decoding - labels will not correspond anymore

            %set chunks, corresponds to the number of runs
            chunks = repmat(1:str2num(sub(s).Nrun), 1, 3)';

            %remove constant features (due to liberal masking)
            ds=cosmo_remove_useless_data(ds);
            
            % IMPORTANT FOR CROSSMODAL DECODING : Demean  every pattern to remove univariate effect differences
            meanPattern = mean(ds.samples,2);  % get the mean for every pattern
            meanPattern = repmat(meanPattern,1,size(ds.samples,2)); % make a matrix with repmat
            ds.samples  = ds.samples - meanPattern; % remove the mean from every every point in each pattern

            % Slice the dataset according to modality
            modIndex = (ds.sa.modality == 1) | (ds.sa.modality==2);
            ds = cosmo_slice(ds, modIndex); 

            % print dataset
            fprintf('Dataset input:\n');
            cosmo_disp(ds);


            % Use the cosmo_cross_validation_measure and set its parameters
            % (classifier and partitions) in a measure_args struct.
            measure = @cosmo_crossvalidation_measure;
            measure_args = struct();

            % Define which classifier to use, using a function handle.
            % Alternatives are @cosmo_classify_{svm,matlabsvm,libsvm,nn,naive_bayes}
            if strcmp(algo, 'lda')
                measure_args.classifier = @cosmo_classify_lda;
            elseif strcmp(algo, 'svm')
                measure_args.classifier = @cosmo_classify_libsvm;
            end 

            %%  Set partition scheme. odd_even is fast; for publication-quality analysis nfold_partitioner is recommended.
            %I use nfold partitions here

            % Other alternatives are:
            % - cosmo_nfold_partitioner    (take-one-chunk-out crossvalidation) %% for within modality 
            % - cosmo_nchoosek_partitioner (take-K-chunks-out  " "). %%for crossmodal decoding
            % - cosmo_oddeven_partitioner(ds_per_run); %%faster
            measure_args.partitions = cosmo_nchoosek_partitioner(ds, 1, 'modality', test); % this is for training in one modality, test on the other etc.
            measure_args.normalization = 'zscore';
            % print measure and arguments
            fprintf('Searchlight measure:\n');
            cosmo_disp(measure);
            fprintf('Searchlight measure arguments:\n');
            cosmo_disp(measure_args);

            %% Define the size of the voxels' sphere you want to take for each searchlight

            % Here I define a neighborhood with approximately 100 voxels in each searchlight.
            nvoxels_per_searchlight=100;
            nbrhood=cosmo_spherical_neighborhood(ds,...
                                    'count',nvoxels_per_searchlight);


            % Run the searchlight
            sl_results = cosmo_searchlight(ds,nbrhood,measure,measure_args);

            % print output dataset
            fprintf('Dataset output:\n');
            cosmo_disp(sl_results);

            % Plot the output
            cosmo_plot_slices(sl_results);

            % Define the name of the file to be saved in the output location
            output_fn = fullfile(output_dir, strcat(working_on, '.nii'));

            % Store results to disc
            cosmo_map2fmri(sl_results, output_fn);

            % Show citation information
            cosmo_check_external('-cite');
        end
    end 
end 