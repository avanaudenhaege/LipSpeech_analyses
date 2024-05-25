%% MULTIDIMENSIONAL SCALING
% Visual representation of relationships between stimuli categories
%
% Adapted from scripts of Filippo
% see https://cosmomvpa.org/_static/publish/demo_fmri_distatis.html for
% detailed explanations 


%% Define dir&paths 

clear;
clc;


mvpa_dir = fileparts(mfilename('fullpath'));
roi_dir = fullfile(mvpa_dir, 'ROIs', 'crossM-clusters');
output_dir = fullfile(mvpa_dir, 'MDS');

lipspeech_dir = fullfile(mvpa_dir, '..'); %root dir
bids_dir = fullfile(lipspeech_dir, 'Lip_BIDS');
stats_dir = fullfile(bids_dir, 'derivatives', 'bidspm-stats');


%% Define subjects, ROIs, tasks to work on
%subjects to analyze
sub_no = [4:24 26 27];

sub_all = sub_data; 
sub_list = sub_all(sub_no);


task_label = {'CrossM'}; %'Aud', 'Vis'
model_label = {'Cons', 'Trialbytrial'}; %'Cons','Speak','Vowels','Trialbytrial'
roi_label = {'rightSTS', 'leftSTS', 'leftMotor'};

algo = 'svm'; %'lda' or 'svm'
val = 'beta'; %or 'beta' ?

%% choose : RDMs already computed (1) or not (0) ? 

rdm_computed = 1; 

if rdm_computed == 0
    %% Set dataset
    % distatis requires a stack of RDMs, one for each subject.
    % Steps:
    % - organize dataset
    % - average different targets
    % - get RDM
    %
    % After distatis, unlfatten function provides data to plot.
    % (still unclear how)

    % To set cosmo mvpa structure, we need:
    % - labels
    % - decoding conditions
    % - indices
    % - pairs
    % they change based on our analyses


    % Extract and stack single subject's RDMs
    for t=1:length(task_label)  
        task = task_label{t};

        %to change when I will automatize the script
        model = model_label{1};
        roi = roi_label{3}; 


        for s = 1:length(sub_list)

            subID= (sub_list(s).id);

            %% Pick ROI to work in
            F = sprintf('ClusterMask_%s*.nii', roi);
            D = dir(fullfile(roi_dir,F));
            roi_img = fullfile(D.folder, D.name);

            %% get 4D file     
            if strcmp(model, 'Trialbytrial')
                %do nothing
            else
                data_img=fullfile(stats_dir, 'mvpaVol', subID, ...
                        strcat(subID, '_task-CrossM_space-IXI549Space_FWHM-2_node-CrossM', model, '_desc-4D_', val, '.nii'));
            end


            %% prepare the structure of the dataset 
            if strcmp(model, 'Trialbytrial')
                %for each stimuli
                targets = repmat(1:54, 1, str2double(sub_list(s).Nrun)); %[1-9 = S1]; [1,2,3,10,11,12,19,20,21 = a]; [1,4,7,10,13,16,19,22,25 = f]
                targets = sort(targets)';           
                chunks1M = repmat(1:str2double(sub_list(s).Nrun), 1, 27)';
                chunks2M = repmat((1:str2double(sub_list(s).Nrun))+str2double(sub_list(s).Nrun), 1, 27)';
                chunks = cat(1, chunks1M, chunks2M);
    %           labels = {'S1fa', 'S1la', 'S1pa', 'S1fi', 'S1li', 'S1pi', 'S1fe', 'S1le', 'S1pe', ...
    %                         'S2fa', 'S2la', 'S2pa', 'S2fi', 'S2li', 'S2pi', 'S2fe', 'S2le', 'S2pe', ...
    %                         'S3fa', 'S3la', 'S3pa', 'S3fi', 'S3li', 'S3pi', 'S3fe', 'S3le', 'S3pe'}';

            else 
                %for cons/speak/vow

                targets=repmat(1:6,1,str2double(sub_list(s).Nrun))'; %there are 3 consonants (mean value for the 9 iteration of the cons), and each is repeated in each run. 
                targets= sort(targets); 

                chunks1M = repmat(1:str2double(sub_list(s).Nrun), 1, 3)';
                chunks2M = repmat((1:str2double(sub_list(s).Nrun))+str2double(sub_list(s).Nrun), 1, 3)';
                chunks = cat(1, chunks1M, chunks2M); %changed here last time !!!! 

                modality = repmat(1:2, 1, str2double(sub_list(s).Nrun)*3)'; %There are 2 modalities (Aud AND Vis), and we have to repeat this *3(for each cons) and *Nrun for each run (usually 19 or 20) = usually 114 or 120
                modality = sort(modality);

                %labels = string(1:str2num(sub_list(s).Nrun));
            end

    data_img1 = strcat(stats_dir, '/', subID, '/task-MVPAAud_space-IXI549Space_FWHM-2_node-MVPAAud', model, '/', subID, '_task-MVPAAud_space-IXI549Space_desc-4D_beta.nii');
    data_img2 = strcat(stats_dir, '/', subID, '/task-MVPAVis_space-IXI549Space_FWHM-2_node-MVPAVis', model, '/', subID, '_task-MVPAVis_space-IXI549Space_desc-4D_beta.nii');
    ds1 = cosmo_fmri_dataset(data_img1, 'mask', roi_img); 
    ds2 = cosmo_fmri_dataset(data_img2, 'mask', roi_img);

    if strcmp(model, 'Trialbytrial')
        ds = ds1;
        ds.samples = [ds1.samples; ds2.samples];
        ds.sa.targets = targets;
        ds.sa.chunks = chunks;
    else
        ds = cosmo_fmri_dataset(data_img, ...
                                            'mask', roi_img,...
                                            'targets',targets, ...
                                            'chunks',chunks);
        ds.samples = [ds1.samples; ds2.samples];
    end

            % Getting rid off zeros
            zeroMask = all(ds.samples == 0, 1);
            ds = cosmo_slice(ds, ~zeroMask, 2);

            % remove constant features
            ds = cosmo_remove_useless_data(ds);

            % Compute average for each unique target (e.g. for each consonant audF, audL, audP, visF, visL, visP)
            % Dataset should have one sample for each target
            ds_mean = cosmo_fx(ds, @(x)mean(x,1), 'targets');
            %ds_mean = cosmo_fx(ds, @(x)mean(x,1), 'chunks'); 



            ds_rdm = cosmo_dissimilarity_matrix_measure(ds_mean);

            % simple sanity check to ensure all attributes are set properly
            cosmo_check_dataset(ds);

            % set chunks (one chunk per subject)
            ds_rdm.sa.chunks = s * ones(size(ds_rdm.samples, 1), 1);
            ds_rdms{s} = ds_rdm;

        end
    end 

elseif rdm_computed == 1
    %% load RDMs if already loaded
    roi = roi_label{1};
    rdm_path = ['/Volumes/T7 Shield/DATA_HD/LipSpeech/LipSpeech_mvpa_SL_svm/MDS/allRDMs_ds_Cons_', roi, '.mat'];
    load(rdm_path)
end 

%% Stack all the RDMs in a single dataset and extract values for MDS plot
% still untouched
allRDMs_ds = cosmo_stack(ds_rdms);

% Run DISTATIS
distatis = cosmo_distatis(allRDMs_ds);

% Compute compromise distance matrix
[compromise_matrix, dim_labels, values] = cosmo_unflatten(distatis, 1);


%% Plot multidimensional scaling: like a scatter-plot but fancy

% Labels
if strcmp(model, 'Trialbytrial')
    labels = {'S1fa_A', 'S1la_A', 'S1pa_A', 'S1fi_A', 'S1li_A', 'S1pi_A', 'S1fe_A', 'S1le_A', 'S1pe_A', ...
                'S2fa_A', 'S2la_A', 'S2pa_A', 'S2fi_A', 'S2li_A', 'S2pi_A', 'S2fe_A', 'S2le_A', 'S2pe_A', ...
                'S3fa_A', 'S3la_A', 'S3pa_A', 'S3fi_A', 'S3li_A', 'S3pi_A', 'S3fe_A', 'S3le_A', 'S3pe_A', ...
                'S1fa_V', 'S1la_V', 'S1pa_V', 'S1fi_V', 'S1li_V', 'S1pi_V', 'S1fe_V', 'S1le_V', 'S1pe_V', ...
                'S2fa_V', 'S2la_V', 'S2pa_V', 'S2fi_V', 'S2li_V', 'S2pi_V', 'S2fe_V', 'S2le_V', 'S2pe_V', ...
                'S3fa_V', 'S3la_V', 'S3pa_V', 'S3fi_V', 'S3li_V', 'S3pi_V', 'S3fe_V', 'S3le_V', 'S3pe_V'}';
else
    %labels = {'audF', 'audL', 'audP', 'visF', 'visL', 'visP'};
    labels = {'F', 'L', 'P', 'F', 'L', 'P'};
end

n_labels = numel(labels);

% % Make figure - compromise matrix
figure();
imagesc(compromise_matrix)
% 
% % Set title and labels on both axes
str = sprintf('DSM ROI: %s', roi);
title(str)
set(gca, 'YTick', 1:n_labels, 'YTickLabel', labels);
set(gca, 'XTick', 1:n_labels, 'XTickLabel', labels);
colorbar % add side legend

% skip if stats toolbox is not present
if cosmo_check_external('@stats',false)

    % Make (tree-like) figure - dendogram 
    figure();
    hclus = linkage(compromise_matrix);
    dendrogram(hclus, 'labels', labels);
    str = sprintf('dendrogram ROI: %s', roi);
    title(str)

    % Make figure - bidimensional plot
    figure();
    F = cmdscale(squareform(compromise_matrix)); % actual mds function
    hold on 
    % Plot Aud data (in BLUE)
    plot(F(1:(n_labels/2),1), F(1:n_labels/2,2), ...
        'linestyle','none', ...
        'marker','o', ...
        'MarkerSize',15,...
        'MarkerEdgeColor',[64/256 141/256 243/256]) % to fill in the dot : 'MarkerFaceColor',[64/256 141/256 243/256]

    % Plot Vis data (in RED)
    plot(F((n_labels/2)+1:n_labels,1), F((n_labels/2)+1:n_labels,2), ...
        'linestyle','none', ...
        'marker','o', ...
        'MarkerSize',15,...
        'MarkerEdgeColor',[241/256 75/256 57/256]) % to fill in the dot : 'MarkerFaceColor',[241/256 75/256 57/256]
    
    % Plot labels for auditory stimuli (first 3) in red
    text(F(1:(n_labels/2),1), F(1:(n_labels)/2,2), labels(1:(n_labels/2)), ...
        'FontSize',14, ...
        'VerticalAlignment','middle', ...
        'HorizontalAlignment','center', ...
        'FontName', 'Avenir', ...
        'Color', [64/256 141/256 243/256]);

    % Plot labels for visual stimuli (last 3) in blue
    text(F((n_labels/2)+1:n_labels,1), F((n_labels/2)+1:n_labels,2), labels((n_labels/2)+1:n_labels), ...
        'FontSize',14, ...
        'VerticalAlignment','middle', ...
        'HorizontalAlignment','center', ...
        'FontName', 'Avenir', ...
        'Color', [241/256 75/256 57/256]);
    
    
    %text(F(:,1), F(:,2), labels,'FontSize',14,'VerticalAlignment','bottom','HorizontalAlignment','left'); %'FontWeight','bold','Color', colors(i,:)    
    str = sprintf('2D MDS plot ROI: %s', roi);
    title(str, 'FontName', 'Avenir')
    
%     % Do not show axis at all
%     axis off

    % Adjust the axis limits
    padding = 0.02; % Adjust this value to get more or less space on edges of the plot
    xlim([min(F(:,1)) - padding, max(F(:,1)) + padding]);
    ylim([min(F(:,2)) - padding, max(F(:,2)) + padding]);
    
    % Set axis color to white for non visible axis
    ax = gca; 
    ax.XColor = [1 1 1]; % White color 
    ax.YColor = [1 1 1]; % White color

    % Optional: Remove the tick labels if not needed
    set(gca, 'XTick', []);
    set(gca, 'YTick', []);


end


