%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% FOLDERS AND DATA SETUP %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('D:\David\MATLAB\raacampbell-shadedErrorBar-0dc4de5\')
addpath('D:\David\MATLAB\suplabel\')
load('D:\Desktop\Predictable_omissions\Preprocessed\predictable_omisions_omit_detrend_ERPs2.mat')
layoutFile = 'D:\David\MATLAB\fieldtrip-20170625\fieldtrip-20170625\template\layout\biosemi64.lay';

data = allData_pred_omit{1,1};


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% CHANNELS AND LAYOUT SETUP %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Frontal Left Electrode region
FL = {'Fp1','AF7','AF3','F7','F5','F3','F1'};
%Frontal Middle Electrode region
FZ = {'Fpz','AFz','Fz','Fp1','FP2','AF3','AF4','F1','F2'};
%Frontal Right Electrode region
FR = {'Fp2','AF8','AF4','F8','F6','F4','F2'};

%Central Left Electrode region
CL = {'FC1','FC5','FC3','C1','C5','C3','CP5','CP3','CP1'};
%Central Middle Electrode region
CZ = {'FCz','Cz','CPz', 'FC1','C1','CP1','FC2','C2','CP2'};
%Central Right Electrode region
CR = {'FC6','FC4','FC2','C6','C4','C2','CP6','CP4','CP2'};

%Parietal Left Electrode region
PL = {'P7','P5','P3','PO7','PO3','O1', 'P1'};
%Parietal Middle Electrode region
PZ = {'Pz','POz','Oz','P1','P2', 'PO3', 'PO4', 'O1','O2'};
%Parietal Right Electrode region
PR = {'P8','P6','P4','PO8','PO4','O2','P2'};


% storing the different regions in a cell array and as well as their 
regions       = {FL;FZ;FR;CL;CZ;CR;PL;PZ;PR};
region_labels = {'FL';'FM';'FR';...
                 'CL';'CM';'CR';...
                 'PL';'PM';'PR'};

% cfg        = [];
% cfg.layout = layoutFile;
% cfg.box    = 'no';
% ft_layoutplot(cfg,data);
% set(findall(gcf,'-property','FontSize'),'FontSize',13)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% MEANS & SEM  for all regions %%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Grand Average across participants
cfg               = [];
cfg.parameter     ='avg';
grand_unpred_omit = ft_timelockgrandaverage(cfg,allData_unpred_omit{:});
grand_pred_omit   = ft_timelockgrandaverage(cfg,allData_pred_omit{:});

%create matricies to store the average ERP per participant, per region for
%the unpredicatble omissions, predictable omissions and difference waveform
all_unpred_omit = zeros(length(allData_pred_omit),length(grand_pred_omit.time),...
                        length(regions));
all_pred_omit   = zeros(size(all_unpred_omit));
all_diff_omit   = zeros(size(all_unpred_omit));

for s = 1:length(allData_pred_omit)
    
    unpred_omit = allData_unpred_omit{s,1};
    pred_omit   = allData_pred_omit{s,1};
    
    cfg           = [];
    cfg.operation = 'subtract';
    cfg.parameter = 'avg';
    diff_omit     = ft_math(cfg,pred_omit,unpred_omit);
    
    for r = 1:length(regions) 

        elecs = matches(pred_omit.label,regions{r});
        
        sub_unpred = mean(unpred_omit.avg(elecs,:),1);
        sub_pred   = mean(pred_omit.avg(elecs,:),1);
        sub_diff   = mean(diff_omit.avg(elecs,:),1);
        
        all_unpred_omit(s,:,r) = sub_unpred;
        all_pred_omit(s,:,r)   = sub_pred;
        all_diff_omit(s,:,r)   = sub_diff;
    end 
end

%calculate the grand averages 
grand_avg_unpred = squeeze(mean(all_unpred_omit,1));
grand_avg_pred   = squeeze(mean(all_pred_omit,1));
grand_avg_diff   = squeeze(mean(all_diff_omit,1));

%calucalte the standard error of the mean
grand_sem_unpred = squeeze(std(all_unpred_omit,[],1)/sqrt(size(all_unpred_omit,1)));
grand_sem_pred   = squeeze(std(all_pred_omit,[],1)/sqrt(size(all_pred_omit,1)));
grand_sem_diff   = squeeze(std(all_diff_omit,[],1)/sqrt(size(all_diff_omit,1))); 

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% PLOTING ALL 9 GRAND AVERAGES %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f = figure();
f.WindowState = 'maximized';

for r = 1:length(regions) 
    
sph(r) = subplot(3,3,r);

shadedErrorBar(data.time,grand_avg_unpred(:,r),grand_sem_unpred(:,r),...
                'lineProps',{'-r','LineWidth',2.25}); hold on
shadedErrorBar(data.time,grand_avg_pred(:,r),grand_sem_pred(:,r),...
                'lineProps',{'--b','LineWidth',2.25}); hold on
            
%changes the axis to be on the origin instead of the bottom and left             
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
ax.XAxis.FontSize = 16;
ax.YAxis.FontSize = 16;

ylim([-1.7 1.7])
yticks([-1.5,-1,-0.5,0.5,1,1.5])
xlim([-0.1 0.4])
% xlabel('Time (s)', 'FontSize', 20)
% ylabel('Amplitude (\muV)','FontSize', 20)
title(region_labels{r,1},'FontSize', 18, 'VerticalAlignment','baseline')
set(sph(r),'position',get(sph(r),'position')+[-0.05 -0.05 0.05 0.05])

end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% PLOTTING TOPOGRAPHIES %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calucalting the grand averages 
cfg               = [];
cfg.parameter     ='avg';
grand_unpred_omit = ft_timelockgrandaverage(cfg,allData_unpred_omit{:});
grand_pred_omit   = ft_timelockgrandaverage(cfg,allData_pred_omit{:});

%Calucaltaing difference waves
cfg           = [];
cfg.operation = 'subtract';
cfg.parameter = 'avg';
grand_diff     = ft_math(cfg,grand_unpred_omit,grand_pred_omit);

% Peak latencies in 20 ms windows for the P1, N1, P2 and N2
lats = {[0.0779 0.0979],[0.1179 0.1379],[0.1659 0.1859],[0.2979 0.3179]};

%cell arrays to store the grand average for each latency
grands_unpred = cell(length(lats),1);
grands_pred   = cell(length(lats),1);

%store the grand averages in the latency windows above
for jj = 1:length(lats)

cfg                 = [];
cfg.parameter       ='avg';
cfg.latency         = lats{1,jj};
grands_unpred{jj,1} = ft_timelockgrandaverage(cfg,allData_unpred_omit{:});
grands_pred{jj,1}   = ft_timelockgrandaverage(cfg,allData_pred_omit{:});

end

%cell arrays to store the difference amplitude relative to the previous ERP
grands_unpred_diffs = cell(length(lats)-1,1);
grands_pred_diffs   = cell(length(lats)-1,1);
grands_diffs = cell(length(lats)-1,1);

%caluclates the difference amplitude relative to the previous ERP  
count = 1;
for jj = 2:length(lats)
 
    unpred_data1 = grands_unpred{jj};
    unpred_data2 =  grands_unpred{jj-1,1};
    unpred_data2.time = unpred_data1.time;
    
    pred_data1 = grands_pred{jj};
    pred_data2 =  grands_pred{jj-1,1};
    pred_data2.time = pred_data1.time;

    
    cfg           = [];
    cfg.operation = 'subtract';
    cfg.parameter = 'avg';
    grands_unpred_diffs{count,1} = ft_math(cfg, unpred_data1, unpred_data2);
    grands_pred_diffs{count,1}   = ft_math(cfg, pred_data1, pred_data2);
    grands_diffs{count,1} = ft_math(cfg,grands_unpred_diffs{count,1},grands_pred_diffs{count,1});
    count = count+1;
end



cfg           = [];
% cfg.xlim      = [0.0779 0.0979]; %time window of the P1 response
% cfg.zlim      = [-0.85 0.85];

% cfg.xlim      = [0.1179 0.1379]; %time window of the N1 response
% cfg.zlim      = [-0.65 0.65];
% 
% cfg.xlim      = [0.1659 0.1859]; %time window of the P2 response
% cfg.zlim      = [-1.05 1.05];

% cfg.xlim      = [0.2979 0.3179]; %time window of the N2 response
% cfg.zlim      = [-0.075 0.075];

cfg.layout    = layoutFile;
cfg.comment   = 'no'; 
cfg.parameter = 'avg'; % the default 'avg' is not present in the data
figure; ft_topoplotER(cfg,grands_diffs{3,1}); colorbar('FontSize',14)
% title(["Unpredictable - Predictable N2-P2"," "]) 

