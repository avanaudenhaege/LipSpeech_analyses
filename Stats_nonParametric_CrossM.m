%% statistical significance in mvpa
% non-parametric technique by combining permutations and bootstrapping

% step1: 
% For each subject, the labels of the different conditions (eg. motion_vertical and motion_horizontal) were permuted,
% and the same decoding analysis was performed. 
% The previous step was repeated 100 times for each subject.

% DONE in our decoding scripts

% step2:
% A bootstrap procedure was applied in order to obtain a group-level null distribution 
% that was representative of the whole group. 
% From each subject’s null distribution, one value was randomly chosen (with replacement) 
% and averaged across all participants. 
% This step was repeated 100,000 times resulting in a group-level null distribution of 100,000 values. 

% step3:
% The statistical significance of the MVPA results was estimated by comparing the observed result 
% to the group-level null distribution. This was done by calculating the proportion of observations 
% in the null distribution that had a classification accuracy higher than the one obtained in the real test.

% step4:
% To account for the multiple comparisons, all p values were corrected using false discovery rate (FDR) correction 

%% set which file, condition and roi label are we testing

% load the .mat file with decoding accuracies
decodTitle= 'trainA_testV'; %'trainV_testA' 'both'


accu_file = ['/Volumes/T7 Shield/DATA_HD/LipSpeech/LipSpeech_mvpa_ROI_svm/Decoding_ROIs_CrossM/Cons/VWFA-indivlpSTS/CrossM_decAccuracies_vwfa-AVoverlap_svm_200vx_MVPA_Cons_beta_', decodTitle, '.mat'];
load(accu_file)



% decodTitle= 'visDirSel';
% decodingConditionList = {'handDown_pinkyThumb_vs_handDown_fingerWrist',...
%     'handUp_pinkyThumb_vs_handUp_fingerWrist',...
%       'visual_vertical_vs_visual_horizontal'};
% decodingCondition='visual_vertical_vs_visual_horizontal';

% decodTitle= 'crossMod_Ext';
% decodingConditionList = {'trainVisual_testTactile','trainTactile_testVision','both'};
% decodingCondition= 'trainVisual_testTactile';


%roiList={'lSTG', 'rSTG', 'lPCG', 'rPCG', 'lMTG', 'rAST'};
roiList = roi_label;

%sub definition
sub_all=sub_data; %load sub informations
subList = [4:24 26 27]; %which subjects to analyze ? 
%sub_info = sub_all(subList);

         
val='beta';%'tmap', 'beta'
smooth='2';
voxNb='200'; %%?? need to keep these 2 lines ?

% number of iterations for group level null distribution
nbIter = 100000;

groupNullDistr=zeros(length(roiList),nbIter);
subAccu=zeros(length(subList),length(roiList));
subSamp = zeros(length(subList), nbIter);

for iRoi=1:length(roiList)
    roiLabel=roiList(iRoi);
    disp(roiList(iRoi))
    
    accuracy = Acc_allROIs(:,iRoi);

    %% STEP 1: DONE

    %% STEP 2: create group null distribution
    timeStart=datestr(now,'HH:MM')

    for iIter = 1:nbIter
    %     disp(iIter)
        %for iAccu=1:length(Acc_allROIs)

            for iSub=1:length(subList)
                
                s = subList(iSub);
                
                %if strcmp(char({accuracy(iAccu).subID}.'),char(subID))==1
                    %check if all the parameters and conditions match             
                    %if strcmp(char({accuracy(iAccu).image}.'), val)==1 && strcmp(num2str([accuracy(iAccu).ffxSmooth].'),smooth)==1 && strcmp(num2str([accuracy(iAccu).choosenVoxNb].'),voxNb)==1   
                        %if strcmp(string({accuracy(iAccu).decodingCondition}.'),decodingCondition)==1  
                            %if strcmp(string({accuracy(iAccu).mask}.'),roiLabel)==1

                                %read the subject level permutations for this roi = nullAcc_allROIs(:,iRoi,s);
                                %pick one decoding accuracy randomly with replacement
                                
                                subSamp(iSub, iIter) = datasample(nullAcc_allROIs(:, iRoi, s),1); 

                            %end

                        %end
                    %end
                %end
            end
        %end
    end

    timeEnd=datestr(now,'HH:MM')
    groupNullDistr(iRoi,:) = mean(subSamp);

    %% STEP 3: check where does the avg accu of the group falls in the group level null ditribution
    % calculate the proportion of values in the group null ditribution which are above the actual decoding
    % accuracy for a one-tailed test. accordingly change for two-tailed test.
    % p = sum(accuracy < acc0) / nbIter; %from Ceren
    % pValue = (sum(Pooled_Perm>accuracy)+1)/(1+NrPermutations); % from Mohamed
    
    %get the actual decoding accuracy 
    subAccu = Acc_allROIs(subList, :);

    subObsPVal(iRoi) = sum(mean(subAccu(:,iRoi))<groupNullDistr(iRoi,:))/nbIter;

end

%% STEP 4: correct the obtained p-value 

% function mafdr([vector of pvalues], BHFDR, 'true') % from Stefania
fdrCorrPVal= mafdr(subObsPVal, 'BHFDR', 'true')
% fdrCorrPValBasic=mafdr(subObsPVal);

%% save the outout

% set output folder/name
pathOutput='/Volumes/T7 Shield/DATA_HD/LipSpeech/LipSpeech_mvpa_ROI_svm/Decoding_ROIs_CrossM/Cons/VWFA-indivlpSTS/stats';
savefileMat = fullfile(pathOutput, ...
                     ['stats', '_CrossM_2ROIs_svm_200vx_Cons_', decodTitle, '_beta.mat']);
                 
% set structure array for keeping the results
% mvpaStats = struct( ...
%             'decodTitle', [], ...
%             'decodCondition', [], ...
%             'roiList', [], ...
%             'groupNullDis', [], ...
%             'obsPVal', [], ...
%             'fdrCorPVal', []);

%store output
mvpaStats.decodTitle = decodTitle;
mvpaStats.imageType = val;
mvpaStats.roiList = roiList; % this tells the order of corresponding p-values
mvpaStats.groupNullDistr = groupNullDistr; % the rows are in the order of Roi list
mvpaStats.subAccu = subAccu;
mvpaStats.obsPVal = subObsPVal; % in the order of roi list
mvpaStats.fdrCorPVal = fdrCorrPVal;

% mat file
save(savefileMat, 'mvpaStats');