clear all

%% Define dir&paths
mvpa_dir = fileparts(mfilename('fullpath'));
output_dir = fullfile(mvpa_dir, 'Decoding_ROIs_CrossM');

lipspeech_dir = fullfile(mvpa_dir, '..'); %root dir
bids_dir = fullfile(lipspeech_dir, 'Lip_BIDS');
stats_dir = fullfile(bids_dir, 'derivatives', 'bidspm-stats');

sub_all=sub_data; %call function in same dir named sub_data.m - with all info on each subject
nsub_all = length(sub_all);

%% Define subjects, ROIs, tasks to work on
%CHANGE HERE BEFORE RUNNING 

sub_included = [4:24 26 27]; %which subjects to analyze ? 
sub_list = sub_all(sub_included);
nsub = length(sub_list);

model_label = {'Cons'}; %'Cons','Speak','Vowels','Trialbytrial'
roi_label = {'vwfa', 'AVoverlap'};%, 'vwfa', 'TVSA'};%'ffa', 'ppaR', 'ppaL', 

modality_label = {'Aud', 'Vis'}; %order of cond in 4D file --> 1 = Auditory - 2 = Visual 
decoding_cond = {'both', 'trainA_testV', 'trainV_testA'};  %'trainA_testV', 'trainV_testA', 'both'

algo = 'svm'; %ideally the same as in the searchlight decoding
val = 'beta'; %or 'beta'

numFeatures=200; % in decoding, you can take the whole ROI, or a fixed number of voxels. (decoding is affected by the size of the ROI). It will choose the 120 most informative voxels for decoding. 

doPermut = 1; %if you don't want to run permutation part, change to 0
nbIter = 100; % number of iterations for the permutation part

% ext = '.nii';
% %ext = '.img';

%% Load data

    
for m=1:length(model_label)
    model = model_label{m};

    %%% preallocate for saving later
    MEAN_confusion_matrix= zeros(3,3); %if you want something else than binary 
    %decoding, you have to change the matrix (if 6 conditions, matrix 6 by 6)
    
    Acc_allROIs = zeros(nsub_all, length(roi_label));
    nFeatures_allROIs = zeros(nsub_all, length(roi_label));          
    nullAcc_allROIs = zeros(nbIter, length(roi_label), nsub_all); 
    p_allROIs = zeros(nsub_all, length(roi_label));

    
    for d = 1:length(decoding_cond) %see the types in decoding_cond list
        decodingCondition = decoding_cond{d};
        if strcmp(decodingCondition, 'trainV_testA')
            test = 1;
        elseif strcmp(decodingCondition, 'trainA_testV')
            test = 2;
        elseif strcmp(decodingCondition, 'both')
            test = [];
        end 
            

        
        for r=1:length(roi_label)
            roi = roi_label{r};

            for s = sub_included 
                sub_name= (sub_all(s).id);
 
                %% Pick ROI to work in
                roi_dir = fullfile(mvpa_dir, 'ROIs'); %reset
                if strcmp(roi, 'TVSA')
                    roi_dir = strcat(roi_dir, '/expand-sphere/', roi);
                    % Pick up the ROI of the subject 
                    P = fullfile(roi_dir, sub_name);
                    F = 'rlabel-sphere10*.nii';
                    D = dir(fullfile(P,F));
                elseif strcmp(roi, 'AVoverlap')
                    if strcmp(sub_name, 'sub-04') || strcmp(sub_name, 'sub-20') || strcmp(sub_name, 'sub-23')
                        % Pick the group overlap ROI 
                        P = fullfile(roi_dir, 'Group-AVphono');
                        F = 'AV-overlap*.nii';
                        D = dir(fullfile(P,F));
                    else
                        roi_dir = strcat(roi_dir, '/expand-sphere/', roi);
                        % Pick the ROI of the subject 
                        P = fullfile(roi_dir, sub_name);
                        F = 'rlabel-sphere10*.nii';
                        D = dir(fullfile(P,F));
                    end
                else
                    roi_dir = strcat(roi_dir, '/expand-withinmask200vx/', roi);
                    % Pick up the ROI of the subject 
                    P = fullfile(roi_dir, sub_name);
                    F = sprintf('rlabel*%s*.nii', roi);
                    D = dir(fullfile(P,F));
                end
                 
                %Get the path and name of this mask to use in the function getDataFromExpansion
                roi_img = fullfile(D.folder, D.name);

                working_on=strcat('best', num2str(numFeatures), ...
                    '_', algo, '_MVPACrossM-', decodingCondition, '_', model); %we usually use 2-3 types of classifiers (lda, svm -which takes more time, not ideal for searchlight-). It changes the way the algorithm works, but we don't go in the detail of that. 

                disp(strcat(sub_name, '___DECODING IN:', roi, '___FOR:', working_on));

                %load 4D file ---> KO wrong 4D !!! 
                data_img=fullfile(stats_dir, 'mvpaVol', sub_name, ...
                    strcat(sub_name, '_task-CrossM_space-IXI549Space_FWHM-2_node-CrossM', model, '_desc-4D_', val, '.nii'));
                data_imgA=fullfile(stats_dir, sub_name, ...
                    strcat('task-MVPAAud_space-IXI549Space_FWHM-2_node-MVPAAud', model), strcat(sub_name, '_task-MVPAAud_space-IXI549Space_desc-4D_', val, '.nii'));
                data_imgV=fullfile(stats_dir, sub_name, ...
                    strcat('task-MVPAVis_space-IXI549Space_FWHM-2_node-MVPAVis', model), strcat(sub_name, '_task-MVPAVis_space-IXI549Space_desc-4D_', val, '.nii'));
                
                    
                
                
                %% set your dataset
                %START WITH COSMO FUNCTIONS
                % !!!! very important part of the script !! 
                % is based on the tsv file that you get in the folder BIDS
                % derivatives/stats/sub-XX/ ...._labelfold.tsv
                % it gives the order of the chunks and targets and modality for the matrix after. 
                % 
                %prepare the targets - first for one modality (1M) then double
                %this vector bcs I have 2modalities in the same 4D file
                targets1M=repmat(1:3,1,str2num(sub_all(s).Nrun))'; %there are 3 consonants (mean value for the 9 iteration of the cons), and each is repeated in each run. 
                targets1M= sort(targets1M); %comment this line if want to see a "random" decoding - labels will not correspond anymore
                targets = cat(1, targets1M, targets1M);

                chunks1M = repmat(1:str2num(sub_all(s).Nrun), 1, 3)';
                chunks2M = repmat((1:str2num(sub_all(s).Nrun))+str2num(sub_all(s).Nrun), 1, 3)';
                chunks = cat(1, chunks1M, chunks2M); %changed here last time !!!! 
                
                modality = repmat(1:2, 1, str2num(sub_all(s).Nrun)*3)'; %There are 2 modalities (Aud AND Vis), and we have to repeat this *3(for each cons) and *Nrun for each run (usually 19 or 20) = usually 114 or 120
                modality = sort(modality);
                
                dsA = cosmo_fmri_dataset(data_imgA, 'mask', roi_img); 
                dsV = cosmo_fmri_dataset(data_imgV, 'mask', roi_img); 

                ds = cosmo_fmri_dataset(data_img, ... 
                                                'mask', roi_img,...
                                                'targets',targets, ...
                                                'chunks',chunks); %uses target and chunk structure that we just did
                                            
                ds.samples = [dsA.samples; dsV.samples]; % concatenate auditory and visual runs               
                                         
                ds.sa.modality = modality; %add modality in the ds info

                % remove constant features (due to liberal masking)
                ds=cosmo_remove_useless_data(ds); %removes voxels that have constant values (e.g. outside the brain, white matter, csf etc)
                
                % IMPORTANT FOR CROSSMODAL DECODING : Demean  every pattern to remove univariate effect differences
                meanPattern = mean(ds.samples,2);  % get the mean for every pattern
                meanPattern = repmat(meanPattern,1,size(ds.samples,2)); % make a matrix with repmat
                ds.samples  = ds.samples - meanPattern; % remove the mean from every every point in each pattern
                
                % Slice the dataset according to modality
                modIndex = (ds.sa.modality == 1) | (ds.sa.modality==2);
                ds = cosmo_slice(ds, modIndex); 
                %what is this step for ?????????????????????????????
                
                
                %% Choose cosmo param and run the decoding
                measure= @cosmo_crossvalidation_measure;

                % Make a struct containing the arguments for the measure:
                args=struct();
                
                %choose classifier to use. Default is lda
                if strcmp(algo, 'svm') 
                    args.child_classifier = @cosmo_classify_svm;
                else 
                    args.child_classifier = @cosmo_classify_lda; 
                end 
                
                args.output='predictions'; %output results : how they are classified or not
                args.normalization = 'zscore';%'demean'; %we decided in the lab to normalize in zscore and not demean. 
                args.feature_selector=@cosmo_anova_feature_selector; %it will select the most informative voxels (features = voxels)

                %if not enough voxels (=features), the script gives an error. so here under, if the ROI is smaller than it will take all the voxels.     
                n=size(ds.samples);
                n_features=n(2);
                
                if  n_features<numFeatures
                    args.feature_selection_ratio_to_keep =  n_features;
                else
                    args.feature_selection_ratio_to_keep = numFeatures;% thats for the feature selection
                end
                
                disp('Num features used :'); %to have a feedback in case there is a bug 
                disp(args.feature_selection_ratio_to_keep);
                
                %args.max_feature_count=6000; %automatycallly Cosmo set the limit at 5000,
                %open if your ROI is bigger
                
                %Divide (partition) the runs for test and training : cross validation
                %partitions = cosmo_nfold_partitioner(ds); %this is the code for within modality decoding : 1-against all, run 2-against all, run 3-against all etc. and then mean the 12 results. 
                partitions = cosmo_nchoosek_partitioner(ds, 1, 'modality', test); % this is for training in one modality, test on the other etc. 
                
                % TO SEE THE CROSSVALIDATION SCHEME : cosmo_disp(partitions);
                % Apply the measure to ds, with args as second argument. Assign the result
                % to the variable 'ds_accuracy'.
                [pred, accuracy] = cosmo_crossvalidate(ds,@cosmo_classify_meta_feature_selection,partitions,args); 
                disp(strcat ('Accuracy:',num2str(accuracy)));
                
                %% Run permutation part (from Jacek https://github.com/JacMatu/ReadSpeech_MVPA/blob/main/code/src/cosmo-mpva/cosmomvpaRoiCrossValidation_ReadSpeech.m)
                
                %% PERMUTATION PART
                if doPermut  == 1

                    % allocate space for permuted accuracies
                    nullAcc = zeros(nbIter, 1);

                    % make a copy of the dataset for null distribution
                    ds0 = ds;

                    % for *nbIter* iterations, reshuffle the labels and compute accuracy
                    for k = 1:nbIter
                        % shuffle with function cosmo_randomize_targets
                         ds0.sa.targets = cosmo_randomize_targets(ds);
  
                        % OR 
                        % manually randomize the targets (because of cross-modal error)
                        % In every modality separately and in every chunk,
                        % randomize the labels
%                           for iChunk=1:max(ds.sa.chunks)
%                               for iTestModality = 1:max(ds.sa.modality)
%                                   ds0.sa.targets(ds.sa.chunks==iChunk & ds.sa.modality==iTestModality) = Shuffle(ds.sa.targets(ds.sa.chunks==iChunk & ds.sa.modality==iTestModality));
%                               end
%                           end
%   
                          %let's start random decoding
                          [~, nullAcc(k)] = cosmo_crossvalidate(ds0, ...
                                                         @cosmo_meta_feature_selection_classifier, ...
                                                         partitions, args);
                    end

                    % sum(A<B) calculates how many times B (vector) is greater than A (number)
                    p = sum(accuracy < nullAcc) / nbIter;
                    fprintf('%d permutations: accuracy=%.3f, p=%.4f\n', nbIter, accuracy, p);
  
                
                end

                
                %% Create confusion matrix

                cosmo_warning('off');
                
                %Create the confusion matrix by comparing targets of the
                %test modality (half of targets) and the predictions from
                %cosmo.
                pred = reshape(pred', 1, []);
                pred = pred(~isnan(pred))';
                if strcmp(decodingCondition, 'both')
                    targ_compare = targets;
                else
                    targ_compare = targets(1:str2num(sub_all(s).Nrun)*3);
                end 
                confusion_matrix=cosmo_confusion_matrix(targ_compare,pred); 
                %we take the ds_accuracy vector and put the results in a vector: how
                %many targets of condition 1 where classified as condition 1, as
                %cdition 2 and as condition 3 ... then same for targets of cdtion 2 and
                %3. 
                
                
                %%% how to calculate accuracy manually, from the confusion matrix :
                %%% uncoment if you want to check : both accuracies should be the same.
%                 sum_diag=sum(diag(confusion_matrix));
%                 sum_total=sum(confusion_matrix(:));
%                 OTHERaccuracy=sum_diag/sum_total; %general accuracy of decoding is calculated like that !! sum_diag/sum_total
%                 disp(strcat ('Accuracy:',num2str(OTHERaccuracy)));
                
                %% save results
                Acc_allROIs(s,r)=accuracy;
                nFeatures_allROIs(s,r) = args.feature_selection_ratio_to_keep; 
                if doPermut ==1
                    nullAcc_allROIs(:,r, s) = nullAcc;
                    p_allROIs(s,r) = p;
                end


                %This is wrong, needs to be changed  !!! 
                MEAN_confusion_matrix=MEAN_confusion_matrix + confusion_matrix; %add the new confusion matrix to the other ones. It is a sum, so after that you need to divide by the number of subjects to get the mean. 

            end 
        end 
        save(strcat(output_dir, '/', model, '/CrossM_decAccuracies_vwfa-AVoverlap_', algo, '_', num2str(numFeatures), 'vx_MVPA_', model, '_', val, '_', decodingCondition), ...
            'Acc_allROIs', 'nFeatures_allROIs', 'MEAN_confusion_matrix', 'roi_label', 'nullAcc_allROIs', 'p_allROIs');
%             csvwrite(strcat(output_dir, '/', model, '/decoding_', working_on),Acc_all)
    end 
end
