clear all

%% Define dir&paths
mvpa_dir = fileparts(mfilename('fullpath')); 

lipspeech_dir = fullfile(mvpa_dir, '..'); %root dir
bids_dir = fullfile(lipspeech_dir, 'Lip_BIDS');
stats_dir = fullfile(bids_dir, 'derivatives', 'bidspm-stats');
output_dir = fullfile(mvpa_dir, 'Decoding_ROIs');

%% Define subjects, ROIs, tasks to work on

sub_all = sub_data;%call function in same dir named sub_data.m - with all info on each subject
nsub_all = length(sub_all);



sub_no = [4:24 26 27]; %which subjects to analyze ? 
sub_list = sub_all(sub_no);
nsub = length(sub_list);

task_label = {'Vis', 'Aud'}; %'Aud', 'Vis'
model_label = {'Cons'}; %'Cons','Speak','Vowels','Trialbytrial'
roi_label = {'vwfa', 'ppaR', 'ppaL', 'ffaR', 'ffaL', 'phonoR', 'phonoL', 'TVSA', 'AVoverlap'}; 

algo = 'svm'; %'lda' or 'svm'
val = 'beta'; %or 'beta' ?

numFeatures=200; % in decoding, you can take the whole ROI, or a fixed number of voxels. (decoding is affected by the size of the ROI). It will choose the 120 most informative voxels for decoding. 

doPermut = 1; %if you don't want to run permutation part, change to 0
nbIter = 100; % number of iterations for the permutation part

% ext = '.nii';
% %ext = '.img';

%% Load data


%for the moment, one ROI at the time
%but the script is made to run several ROIs with ROI = {'...', '...', '...'};

%data (tmap) to use ? in the directory 4D-files 


%%%To run all the possible pair of binary classification in one loop
for t=1:length(task_label)  
    task = task_label{t};
    
    for m=1:length(model_label)
        model = model_label{m};
        
        %%% preallocate for saving later
        MEAN_confusion_matrix= zeros(3,3); %if you want something else than binary 
        %decoding, you have to change the matrix (if 6 conditions, matrix 6 by 6)

        Acc_allROIs = zeros(nsub_all, length(roi_label));
        nFeatures_allROIs = zeros(nsub_all, length(roi_label));          
        nullAcc_allROIs = zeros(nbIter, length(roi_label), nsub_all); 
        p_allROIs = zeros(nsub_all, length(roi_label));
        

        for r=1:length(roi_label)
            roi = roi_label{r};
            
            
            for s = sub_no 
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
                elseif strcmp(roi, 'ffaL')    
                    roi_dir = '/Volumes/T7 Shield/DATA_HD/LipSpeech/LipSpeech_mvpa_ROI_svm/ROIs/Group-ffaL-sphereExclCerebellum'; 
                    P = fullfile(roi_dir); 
                    F = 'ffaL-MinusGroupvwfa-306voxels.nii'; 
                    D = dir(fullfile(P,F)); 
                else
                    roi_dir = strcat(roi_dir, '/expand-withinmask200vx/', roi);
                    % Pick up the ROI of the subject 
                    P = fullfile(roi_dir, sub_name);
                    F = sprintf('rlabel*%s*.nii', roi);
                    D = dir(fullfile(P,F));
                end
                 
                %Get the path and name of this mask to use in the function getDataFromExpansion
                roi_img = fullfile(D.folder, D.name);
                
                
% %                 %see if roi image already resliced or not
% %                 lookfor = dir(fullfile(roi_dir, sub_name,['r', roi, '*']));
% %                 if isfile(fullfile(roi_dir, sub_name, lookfor.name))
% %                     %if yes, take this resliced image
% %                     roi_img = fullfile(roi_dir, sub_name, lookfor.name);
% %                 else 
% %                     %if no, reslice roi image
% %                     roifile = dir(fullfile(roi_dir, sub_name,[roi, '*']));
% %                     imageToCheck = fullfile(roi_dir, sub_name, roifile(2).name);              
% %                     referenceImage = fullfile(mvpa_dir, 'ROIs/Masks/reference_image_for_reslice.nii'); % choose any image from the scanner that will be used with the mask.
% %                     roi_img = resliceRoiImages(referenceImage, imageToCheck, 0);
% %                 end 
                
                
                working_on=strcat('best', num2str(numFeatures), ...
                    '_', algo, '_MVPA', task, '_', model); %we usually use 2-3 types of classifiers (lda, svm -which takes more time, not ideal for searchlight-). It changes the way the algorithm works, but we don't go in the detail of that. 
        
                disp(strcat(sub_name, '___DECODING IN:', roi, '___FOR:', working_on));
                
                %load 4D file
                data_img=fullfile(stats_dir, sub_name, ...
                    strcat('/task-MVPA', task, '_space-IXI549Space_FWHM-2_node-MVPA', task, model), ...
                    strcat(sub_name, '_task-MVPA', task, '_space-IXI549Space_desc-4D_', val, '.nii'));
                

                % !!!! very important part of the script !! 
                % is based on the tsv file that you get in the folder BIDS
                % derivatives/stats/sub-XX/ ...._labelfold.tsv
                % it gives the order of the chunks and targets for the matrix after. 
                % 
                %prepare the targets 
                targets=repmat(1:3,1,str2num(sub_all(s).Nrun))'; %there are 3 consonants (mean value for the 9 iteration of the cons), and each is repeated in each run. 
                targets= sort(targets); %comment this line if want to see a "random" decoding - labels will not correspond anymore
                chunks = repmat(1:str2num(sub_all(s).Nrun), 1, 3)';
                
                ds = cosmo_fmri_dataset(data_img, ... 
                                                'mask', roi_img,...
                                                'targets',targets, ...
                                                'chunks',chunks); %uses target and chunk structure that we just did


                % remove constant features (due to liberal masking)
                ds=cosmo_remove_useless_data(ds); %removes voxels that have constant values (e.g. outside the brain, white matter, csf etc)

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

             % if not enough voxels (=features), the script gives an error. so here under, if the ROI is smaller than it will take all the voxels.     
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
                partitions = cosmo_nfold_partitioner(ds); %how to divide the runs ??? we can do even and odd but it gives only 1 result. But what we prefer to do is to do run 1-against all, run 2-against all, run 3-against all etc. and then mean the 12 results. 

                % Apply the measure to ds, with args as second argument. Assign the result
                % to the variable 'ds_accuracy'.
                ds_accuracy = cosmo_crossvalidate(ds,@cosmo_meta_feature_selection_classifier,partitions,args);
                ds_accuracy = reshape(ds_accuracy', 1, []);
                ds_accuracy = ds_accuracy(~isnan(ds_accuracy))';
                %ds_accuracy is a list of each target (in the order defined by
                %label_fold.tsv) and the result of the classifier for each target. For
                %example, target one has been classified as a condition "2"
                %(uncorrect, it is a condition 1). 

                cosmo_warning('off');
                %if sub(isub).Nrun=='4'
                %ds_accuracy =[ds_accuracy(1:12,1);ds_accuracy(13:24,2);ds_accuracy(25:36,3);ds_accuracy(37:48,4)];
                %else
                %ds_accuracy =[ds_accuracy(1:12,1);ds_accuracy(13:24,2);ds_accuracy(25:36,3);ds_accuracy(37:48,4);ds_accuracy(49:60,5)];
                %end
                confusion_matrix=cosmo_confusion_matrix(targets,ds_accuracy); 
                %we take the ds_accuracy vector and put the results in a vector: how
                %many targets of condition 1 where classified as condition 1, as
                %cdition 2 and as condition 3 ... then same for targets of cdtion 2 and
                %3. 


                %confusion_matrix=cosmo_confusion_matrix(ds_accuracy.sa.targets,ds_accuracy.samples);
                sum_diag=sum(diag(confusion_matrix));
                sum_total=sum(confusion_matrix(:));
                accuracy=sum_diag/sum_total; %general accuracy of decoding is calculated like that !! sum_diag/sum_total
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
  
                          %let's start random decoding
                          [~, nullAcc(k)] = cosmo_crossvalidate(ds0, ...
                                                         @cosmo_meta_feature_selection_classifier, ...
                                                         partitions, args);
                    end

                    % sum(A<B) calculates how many times B (vector) is greater than A (number)
                    p = sum(accuracy < nullAcc) / nbIter;
                    fprintf('%d permutations: accuracy=%.3f, p=%.4f\n', nbIter, accuracy, p);
  

                
                end

                %% save results
                Acc_allROIs(s,r)=accuracy;
                nFeatures_allROIs(s,r) = args.feature_selection_ratio_to_keep; 
                if doPermut ==1
                    nullAcc_allROIs(:,r, s) = nullAcc;
                    p_allROIs(s,r) = p;
                end
                
                

                % print classification accuracy in terminal window
            %    fprintf('%s\n',desc);
                % Show the result
                %fprintf('\nOutput dataset (with classification accuracy)\n');
                % Show the contents of 'ds_accuracy' using 'cosmo_disp'
                %cosmo_disp(ds_accuracy);

                MEAN_confusion_matrix=MEAN_confusion_matrix + confusion_matrix; %add the new confusion matrix to the other ones. It is a sum, so after that you need to divide by the number of subjects to get the mean. 

                
            end 
            save(strcat(output_dir, '/', model, '/decAccuracies_allROIs_', algo, '_', num2str(numFeatures), 'vx_MVPA', task, '_', model, '_', val), ...
                'Acc_allROIs', 'nFeatures_allROIs', 'MEAN_confusion_matrix', 'roi_label', 'nullAcc_allROIs', 'p_allROIs');%, 
%             csvwrite(strcat(output_dir, '/', model, '/decoding_', working_on),Acc_all)
        end 
    end

end%%for iroi