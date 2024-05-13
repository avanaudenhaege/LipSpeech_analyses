% (C) Copyright 2019 bidspm developers

%%% Script for statistical univariate analysis with bidspm()
%%% to run, needs bidspm() installed and path added and saved. 
clear;
clc;

%% initialize bidspm() for this matlab session
%addpath(/Applications/bidspm);
bidspm();

%% set up BIDS folders and path
this_dir = fileparts(mfilename('fullpath'));

root_dir = fullfile(this_dir, '..');
bids_dir = fullfile(root_dir);
output_dir = fullfile(root_dir, 'derivatives');
preproc_dir = fullfile(root_dir, 'derivatives', 'bidspm-preproc');

%% TO BE DEFINED subjects and tasks to analyze

subject_label = {'04', '05', '06', '07', '08', '09', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24', '26', '27'}; %can write several subj {'04', '05'}

task = 'Vis'; % Phono - Vis - TVSA

model_file = fullfile(this_dir, 'models', strcat('model-', task, 'Loc.json'));

%model_file = fullfile(this_dir, 'models', 'model-TVSALoc.json');
%model_file = fullfile(this_dir, 'models', 'model-VisLoc.json');
%model_file = fullfile(this_dir, 'models', 'model-PhonoLoc.json');
%model_file = fullfile(this_dir, 'models', 'model-MVPAAud.json');
%model_file = fullfile(this_dir, 'models', 'model-MVPAVis.json');


%% Define results saved as output

results.nodeName = 'subject_level';

if strcmp(task, 'Vis')
     results.name = {'face_gt_others', 'word_gt_others', 'house_gt_others'};
     results.montage.slices = -16:2:0;
elseif strcmp(task, 'Phono')
    results.name = {'SYL_gt_SCR', 'SCR_gt_SYL'};
    results.montage.slices = -10:2:10;
elseif strcmp(task, 'TVSA')
     results.name = {'VL_gt_NL', 'NL_gt_VL'};
     results.montage.slices = -10:2:10;
end 

results.png = true();
results.csv = true();
results.montage.do = true();
results.montage.orientation = 'axial';
results.montage.background = struct('suffix', 'T1w', ...
                                    'desc', 'preproc', ...
                                    'modality', 'anat');
%(add opt.results parameters to choose contrasts/corrections..)                                 
opt.results = results;


%% run bidspm stats
% % bidspm(bids_dir, output_dir, 'subject', ...
% %        'participant_label', subject_label, ...
% %        'action', 'stats', ...
% %        'preproc_dir', preproc_dir, ...
% %        'model_file', model_file, ...
% %        'space', {'individual', 'IXI549Space'}, ...
% %        'fwhm', 6, ...
% %        'options', opt);
       
%'action', 'contrasts', ... %to run only the contrasts


%% dataset level

opt.results = struct('nodeName',  'dataset_level', ...
                     'name', results.name, ...
                     'MC', 'none', ...
                     'p', 0.05, ...
                     'k', 10, ...
                     'nidm', true);
% % opt.results = struct('nodeName',  'dataset_level', ...
% %                      'name', {'face_gt_others', 'word_gt_others', 'house_gt_others'}, ...
% %                      'MC', 'none', ...
% %                      'p', 0.05, ...
% %                      'k', 10, ...
% %                      'nidm', true);

bidspm(bids_dir, output_dir, 'dataset', ...
       'participant_label', subject_label, ...
       'action', 'stats', ...
       'preproc_dir', preproc_dir, ...
       'model_file', model_file, ...
       'options', opt);
