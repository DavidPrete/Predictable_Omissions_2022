%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% DIR AND DATA SETUP %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('\\trainorserv.mcmaster.ca\trainorlab\David_Prete\Infant_MMN_Omitt\Scripts\AB_NONGUI')
addpath('D:\Desktop\Predictable_omissions\EEG\scripts\')
addpath('D:\David\MATLAB\fieldtrip-20170625\fieldtrip-20170625\');
addpath('D:\David\MATLAB\raacampbell-shadedErrorBar-0dc4de5\')
ft_defaults

load('D:\David\MATLAB\fieldtrip-20170625\fieldtrip-20170625\template\neighbours\biosemi64_neighb.mat')
load('D:\Desktop\Predictable_omissions\Preprocessed\predictable_omisions_omit_detrend_ERPs2.mat')
savedir ='D:\Desktop\Predictable_omissions\Cluster_Stats\';
if ~exist(savedir,'dir'); mkdir(savedir); end

subjs  = length(allData_unpred_omit);
labels = allData_pred_omit{1,1}.label;
time   = allData_pred_omit{1,1}.time;

%%
%%%%%%%%%%%%%%%%%%%%%%%
%%%%% CHAN INFO %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%

FL = {'Fp1','AF7','AF3','F7','F5','F3','F1'};
FZ = {'Fpz','AFz','Fz','Fp1','FP2','AF3','AF4','F1','F2'};
FR = {'Fp2','AF8','AF4','F8','F6','F4','F2'};

CL = {'FC1','FC5','FC3','C1','C5','C3','CP5','CP3','CP1'};
CZ = {'FCz','Cz','CPz', 'FC1','C1','CP1','FC2','C2','CP2'};
CR = {'FC6','FC4','FC2','C6','C4','C2','CP6','CP4','CP2'};

PL = {'P7','P5','P3','PO7','PO3','O1', 'P1'};
PZ = {'Pz','POz','Oz','P1','P2', 'PO3', 'PO4', 'O1','O2'};
PR = {'P8','P6','P4','PO8','PO4','O2','P2'};

regions       = {FL;FZ;FR;CL;CZ;CR;PL;PZ;PR};
regions_names = {'Frontal Left';'Frontal Middle';'Frontal Right';...
                 'Central Left';'Central Middle';'Central Right';...
                 'Parietal Left';'Parietal Middle';'Parietal Right'};

[~,P1_start] = min(abs(time-0.05));
[~,P1_end]   = min(abs(time-0.12));

[~,P2_start] = min(abs(time-0.120));
[~,P2_end]   = min(abs(time-0.220));

[~,N1_start] = min(abs(time-0.09));
[~,N1_end]   = min(abs(time-0.150));

[~,N2_start] = min(abs(time-0.250));
[~,N2_end]   = min(abs(time-0.35));

[~,base_start] = min(abs(time-(-0.100)));
[~,base_end]   = min(abs(time-0));

time_indicies = [P1_start, P1_end;N1_start, N1_end;P2_start, P2_end;N2_start, N2_end;];


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% FINDING THE PEAK LATENCIES %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

allData = [allData_pred_omit; allData_unpred_omit];

cfg           = [];
cfg.parameter ='avg';
grand_avg     = ft_timelockgrandaverage(cfg,allData{:});
all_avg = zeros(length(allData),length(time),length(regions));


P1_unpred_locs = zeros(length(allData_pred_omit),length(regions));
N1_unpred_locs = zeros(length(allData_pred_omit),length(regions));
P2_unpred_locs = zeros(length(allData_pred_omit),length(regions));
N2_unpred_locs = zeros(length(allData_pred_omit),length(regions));

P1_pred_locs = zeros(length(allData_pred_omit),length(regions));
N1_pred_locs = zeros(length(allData_pred_omit),length(regions));
P2_pred_locs = zeros(length(allData_pred_omit),length(regions));
N2_pred_locs = zeros(length(allData_pred_omit),length(regions));

for s = 1:length(allData_unpred_omit)
    unpred_omit = allData_unpred_omit{s,1};
    pred_omit   = allData_pred_omit{s,1};
    
    for rr = 1:length(regions) 
        elecs       = matches(pred_omit.label,regions{rr});
        unpred_data = mean(unpred_omit.avg(elecs,:),1);
        pred_data   = mean(pred_omit.avg(elecs,:),1);

            
        [~,P1_unpred] = findpeaks(unpred_data(P1_start:P1_end),'NPeaks',1,'SortStr','descend');
        [~,N1_unpred] = findpeaks(-1*unpred_data(N1_start:N1_end),'NPeaks',1,'SortStr','descend');
        [~,P2_unpred] = findpeaks(unpred_data(P2_start:P2_end),'NPeaks',1,'SortStr','descend');
        [~,N2_unpred] = findpeaks(-1*unpred_data(N2_start:N2_end),'NPeaks',1,'SortStr','descend');
        
        [~,P1_pred] = findpeaks(pred_data(P1_start:P1_end),'NPeaks',1,'SortStr','descend');
        [~,N1_pred] = findpeaks(-1*pred_data(N1_start:N1_end),'NPeaks',1,'SortStr','descend');
        [~,P2_pred] = findpeaks(pred_data(P2_start:P2_end),'NPeaks',1,'SortStr','descend');
        [~,N2_pred] = findpeaks(-1*pred_data(N2_start:N2_end),'NPeaks',1,'SortStr','descend');
        
        
        if isempty(P1_unpred)
            P1_unpred = NaN;
        end
        if isempty(P1_pred)
            P1_pred = NaN;
        end
        if isempty(N1_unpred)
            N1_unpred = NaN;
        end
        if isempty(N1_pred)
            N1_pred = NaN;
        end
        if isempty(P2_unpred)
            P2_unpred = NaN;
        end
        if isempty(P2_pred)
            P2_pred = NaN;
        end
        if isempty(N2_unpred)
            N2_unpred = NaN;
        end
        if isempty(N2_pred)
            N2_pred = NaN;
        end
            
        
        P1_unpred_locs(s,rr) = P1_unpred;
        N1_unpred_locs(s,rr) = N1_unpred;
        P2_unpred_locs(s,rr) = P2_unpred;
        N2_unpred_locs(s,rr) = N2_unpred;

        P1_pred_locs(s,rr) = P1_pred;
        N1_pred_locs(s,rr) = N1_pred;
        P2_pred_locs(s,rr) = P2_pred;
        N2_pred_locs(s,rr) = N2_pred;

        
    end
end


P1_unpred_locs = replace_nan_with_avg(P1_unpred_locs, P1_start);
N1_unpred_locs = replace_nan_with_avg(N1_unpred_locs, N1_start);
P2_unpred_locs = replace_nan_with_avg(P2_unpred_locs, P2_start);
N2_unpred_locs = replace_nan_with_avg(N2_unpred_locs, N2_start);

P1_pred_locs = replace_nan_with_avg(P1_pred_locs, P1_start);
N1_pred_locs = replace_nan_with_avg(N1_pred_locs, N1_start);
P2_pred_locs = replace_nan_with_avg(P2_pred_locs, P2_start);
N2_pred_locs = replace_nan_with_avg(N2_pred_locs, N2_start);


base_pred_locs   = [];
base_unpred_locs = [];
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% CALCULATING MEAN AMPLITUDE %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

P1_unpred_amps = zeros(length(allData_pred_omit),length(regions));
N1_unpred_amps = zeros(length(allData_pred_omit),length(regions));
P2_unpred_amps = zeros(length(allData_pred_omit),length(regions));
N2_unpred_amps = zeros(length(allData_pred_omit),length(regions));

P1_pred_amps = zeros(length(allData_pred_omit),length(regions));
N1_pred_amps = zeros(length(allData_pred_omit),length(regions));
P2_pred_amps = zeros(length(allData_pred_omit),length(regions));
N2_pred_amps = zeros(length(allData_pred_omit),length(regions));

base_pred_amps   = zeros(length(allData_pred_omit),length(regions));
base_unpred_amps = zeros(length(allData_pred_omit),length(regions));


for s = 1:length(allData_unpred_omit)
    unpred_omit = allData_unpred_omit{s,1};
    pred_omit   = allData_pred_omit{s,1};
    
    for rr = 1:length(regions) 
        elecs       = matches(pred_omit.label,regions{rr});
        unpred_data = mean(unpred_omit.avg(elecs,:),1);
        pred_data   = mean(pred_omit.avg(elecs,:),1);

        P1_unpred = extract_means(unpred_data,P1_unpred_locs(s,rr));
        N1_unpred = extract_means(unpred_data,N1_unpred_locs(s,rr));
        P2_unpred = extract_means(unpred_data,P2_unpred_locs(s,rr));
        N2_unpred = extract_means(unpred_data,N2_unpred_locs(s,rr));
        
        P1_pred = extract_means(pred_data,P1_unpred_locs(s,rr));
        N1_pred = extract_means(pred_data,N1_unpred_locs(s,rr));
        P2_pred = extract_means(pred_data,P2_unpred_locs(s,rr));
        N2_pred = extract_means(pred_data,N2_unpred_locs(s,rr));
        
        
        P1_unpred_amps(s,rr) = P1_unpred;
        N1_unpred_amps(s,rr) = N1_unpred;
        P2_unpred_amps(s,rr) = P2_unpred;
        N2_unpred_amps(s,rr) = N2_unpred;

        P1_pred_amps(s,rr) = P1_pred;
        N1_pred_amps(s,rr) = N1_pred;
        P2_pred_amps(s,rr) = P2_pred;
        N2_pred_amps(s,rr) = N2_pred;
        
        base_pred_amps(s,rr)   = mean(pred_data(base_start:base_end));
        base_unpred_amps(s,rr) = mean(unpred_data(base_start:base_end));
        
    end
end
