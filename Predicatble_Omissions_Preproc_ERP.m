
datadir ='D:\Desktop\Predictable_omissions\EEG\';
savedir ='D:\Desktop\Predictable_omissions\Preprocessed\';
if ~exist(savedir,'dir'); mkdir(savedir); end

% add path for toolboxes and load in necessary files. 
addpath('\\trainorserv.mcmaster.ca\trainorlab\David_Prete\Infant_MMN_Omitt\Scripts\AB_NONGUI')
addpath('\\trainorserv.mcmaster.ca\trainorlab\David_Prete\Predicatble_Omissions\EEG\scripts\')
addpath('D:\David\MATLAB\fieldtrip-20170625\fieldtrip-20170625\');
%set fieldtrip defautls
ft_defaults

%load in the trigger labels
load('D:\Desktop\Predictable_omissions\Predictable_Omissions_stim\all_trigs.mat');
%load the output of the ICA analysis
load("D:\Desktop\Predictable_omissions\ICA\predticable_omissions_ica_mix_topo.mat");
%load cell array contating ICA comps to remove (found via visual inspection) 
load("D:\Desktop\Predictable_omissions\ICA\comps_to_remove2.mat");


%create Fieldtrip neighoubors based on Fieldtrip template of electrode caps
layoutFile       = 'D:\David\MATLAB\fieldtrip-20170625\fieldtrip-20170625\template\layout\biosemi64.lay';
cfg              = [];
cfg.method       = 'distance';
cfg.neighbordist = 4; 
cfg.layout       = layoutFile;
neighbours       = ft_prepare_neighbours(cfg);
labels           = {neighbours.label}';

%get the file name that will be read in and setup the output file name 
files  = dir([datadir,'*.bdf']);
output = [savedir, 'predicability_omissions_adults.mat'];



%%

%cell array all the preprocessing data will be stored 
allData_cleaned = cell(length(files),1);

for s = 1:length(files)

    disp('--------------------------------')
    disp(['Preprocessing: ', num2str(s)])
    disp('--------------------------------')
    output_comp = [savedir,files(s).name(1:end-4),'_ica_comps.mat'];
    
%Getting the file name of the EEG data
%loading the data into fieldtrip
disp('Loading data...')
    fileName = strcat(files(s).folder,'\',files(s).name);
    
    cfg         = [];
    cfg.dataset = fileName;
    data        = ft_preprocessing(cfg);
    
    %gets the timing of the trggers in the EEG data based on the STATUS channel
    status_chan  = diff(data.trial{1,1}(end,:));
    status_trigs = find(status_chan==1);
    
    %Exceptions due to technical errors
    if s==8
        all_trigs{s,1} = all_trigs{s,1}(17:end); 
    elseif s==14
        all_trigs{s,1} = all_trigs{s,1}(2:end);
        data.label([65:96]) = [];
        data.trial{1,1}(65:96,:)=[];
    elseif s==19
        all_trigs{s,1} = all_trigs{s,1}(2:end);
    end
    data.label(1:64) = labels;
    
    
    % define trials and store the start, end and offset in trl 
    cfg          = [];
    cfg.dataset  = fileName;
    cfg.fsample  = data.fsample;
    cfg.prestim  = -0.100; % 100 ms baseline
    cfg.poststim = 0.400;  % 400 ms post event onset  
    cfg.sample   = status_trigs;
    trl          = Predict_Omit_definetrial(cfg);
    
    %exception due to techincal errors 
    if s==8 || s==19
        trl(1,:) = [];
    elseif s ==21
        trl(1:9,:) =[];
    end


 % Highpass butterworth filter
    cfg                 = [];
    cfg.hpfilter   = 'yes';
    cfg.hpfreq     = 0.5;
    cfg.hpfiltord  = 4;
    cfg.channel    = {'eeg'};
    data           = ft_preprocessing(cfg, data);
    
 % Lowpass butter worth filter 
    cfg = [];
    cfg.lpfilter   = 'yes';
    cfg.lpfreq     = 25;
    cfg.lpfilttype = 'but';
    cfg.lpfiltord  = 4;
    data           = ft_preprocessing(cfg, data);
    
 % add the trial codes to the data 
    data.trialinfo = all_trigs{s,1};
    trl(:,4)       = all_trigs{s,1};
    
 % This conducts Artifcat blocking. It is a type of artifact correction for
 % high amplitude artifact such as blink, movements, jaw clenching. It is
 % based on the paper by Mourad, Reilly, Debruin and Hasey, 2007
    disp( 'Artifact Blocking...')
    Parameters            = [];
    Parameters.Approach   = 'Window'; 
    Parameters.Threshold  = 75; %voltage threshold 
    Parameters.Fs         = data.fsample;
    Parameters.WindowSize = 5; % unit in second
    Parameters.InData     = data.trial{1}; % may have to exclude the high-pass artifact before AB
    Parameters            = Run_AB(Parameters);
    data.trial{1}         = [Parameters.OutData];

    disp('Re-referencing to common average...')
    cfg            = []; 
    cfg.refchannel = {'all'};
    cfg.reref      = 'yes';
    data           = ft_preprocessing(cfg,data);
    
 % epoching dataing into trials 
    cfg       = [];
    cfg.trl   = trl;
    data      = ft_redefinetrial(cfg,data);

 % down-sampling    
    cfg            = [];
    cfg.resamplefs = 250;
    cfg.detrend    = 'yes';
    data           = ft_resampledata(cfg,data);

% Run ICA using unmixing and topographic matrices from previously
% calulated ICA 

    cfg            = [];
    cfg.method     = 'runica';
    cfg.unmixing   = all_ica{s,1}.unmixing;
    cfg.topolabel  = all_ica{s,1}.topolabel;
%     cfg.runica.pca = rank(data.trial{1,1})-1;
    comps          = ft_componentanalysis(cfg,data);
    
% Remove artifact components. Components were labeled as artifacts based on 
% visual inspection using ft_databrowser
    cfg           = [];
    cfg.component = comps_to_remove{1,s};
    data          = ft_rejectcomponent(cfg,comps,data);
%     clear comps
    
% find trials with range greater than 100 micro volts in any channel and
% removes it from the data. stores the trials labeled as artifacts
trials_to_reject = [];

for jj = 1:length(data.trial)
    trial = data.trial{jj};
    maximums = max(trial,[],2 );
    minimums = min(trial,[],2 );
    ranges   = abs(minimums-maximums);
    reject   = sum(ranges>=100);
    
    if reject>0
        trials_to_reject   = [trials_to_reject,jj];
    end
end

data.trial(trials_to_reject)     = [];
data.time(trials_to_reject)      = [];
data.trialinfo(trials_to_reject) = [];
data.artifacts                   = trials_to_reject;



% stores preprocessed data in cell array to be saved after preprocessing
% all participants data 
allData_cleaned{s,1} = data;

end 

% 
disp('Saving data...')
save(output,'allData_cleaned','-v7.3');
disp('Data Saved!')
% 
% 
%%
%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% ERP ANALYSIS %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%

savedir= 'D:\Desktop\Predictable_omissions\Preprocessed\';
if ~exist(savedir,'dir'); mkdir(savedir); end

load('D:\Desktop\Predictable_omissions\Preprocessed\predicability_omissions_adults2.mat')

allData_unpred_omit = cell(length(allData_cleaned),1);
allData_pred_omit   = cell(length(allData_cleaned),1);
allData_unpred_tone = cell(length(allData_cleaned),1);
allData_pred_tone   = cell(length(allData_cleaned),1);

for s = 1:length(allData_cleaned) 

    data = allData_cleaned{s,1};
    
    %linear detrending
    cfg            = []; 
    cfg.detrend    = 'yes';
    cfg.continuous = 'no';
    data           = ft_preprocessing(cfg,data);

    % separating the trials into the differnet stimuli types
    
    cfg              = [];
    cfg.keeptrials   = 'no';
    % unpredictable omissions
    cfg.trials       = data.trialinfo==4; 
    data_unpred_omit = ft_timelockanalysis(cfg,data);
    % tones in unpredictable block
    cfg.trials       = data.trialinfo==1;
    data_unpred_tone = ft_timelockanalysis(cfg,data);
    % predictable omissions
    cfg.trials       = data.trialinfo==8;
    data_pred_omit   = ft_timelockanalysis(cfg,data);
    % tones in predictable block
    cfg.trials       = data.trialinfo==2;
    data_pred_tone   = ft_timelockanalysis(cfg,data);
    
    
    %baseline correction 
    cfg              = [];
    cfg.baseline     = [-.1 0];
    cfg.parameter    = 'avg';
    data_unpred_omit = ft_timelockbaseline(cfg,data_unpred_omit);
    data_pred_omit   = ft_timelockbaseline(cfg,data_pred_omit);
    data_unpred_tone = ft_timelockbaseline(cfg,data_unpred_tone);
    data_pred_tone   = ft_timelockbaseline(cfg,data_pred_tone);

    % storing the ERPs for the different event types in different cell array
    allData_unpred_omit{s,1} = data_unpred_omit;
    allData_pred_omit{s,1}   = data_pred_omit;
    allData_unpred_tone{s,1} = data_unpred_omit;
    allData_pred_tone{s,1}   = data_pred_omit;


end



disp('Saving data...')
save([savedir,'predictable_omisions_omit_detrend_ERPs2.mat'],'allData_unpred_omit','allData_pred_omit','-v7.3');
save([savedir,'predictable_omisions_tone_detrend_ERPs2.mat'],'allData_unpred_tone','allData_pred_tone','-v7.3');
disp('Data Saved!')



 
