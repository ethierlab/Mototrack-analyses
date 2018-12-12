

% select data table to analyse
% data_table = pre_initiation;
data_table = hold_time;

num_bins = size(data_table.spikes_bin{1},1);
bins_dur = 0.05; % 50ms
timeframe = (-(num_bins-1):0)*bins_dur;

trial_types = {'no reward','single','jackpot'};
trial_colors = {'k','b','r'};

sessions = unique(data_table.session_index);
num_sess = length(sessions);

% array to store averaged FR for each cluster
trial_data = cell(1,3); % {no reward, single, jackpot}
mean_data = cell(1,3); % {no reward, single, jackpot}

num_trials = zeros(num_sess,3);
clust_sidx = [];

for s = 1:num_sess

    data_s = data_table(data_table.session_index == sessions(s), :);
    num_clust = size(data_s.spikes_bin{1},2);
    
    % vector to identify which session each cluster comes from
    clust_sidx = [ clust_sidx sessions(s)*ones(1,num_clust)];
    
    %% extract FR data from different trial types

    for r = -1:1
        % in jackpot_bin, no_rew = -1, single=0, jackpot=1
        trial_data(sessions(s),r+2) = {data_s.spikes_bin(cell2mat(data_s.jackpot_bin)==r)};
    end
       
    %% average FR for all units independently, for each trial type
    mean_data_tmp = cell(1,3); % {no reward, single, jackpot}
    null_data = zeros(num_bins,num_clust);
    for r = 1:3
        mean_data_tmp{1,r} = null_data; % start with zeros
        num_trials(sessions(s),r) = size(trial_data{sessions(s),r},1);
        for t = 1:num_trials(sessions(s),r)
            mean_data_tmp{1,r} = mean_data_tmp{1,r} + trial_data{sessions(s),r}{t,:}./num_trials(sessions(s),r); %add FR/num_trials from each trials to get mean
        end
        mean_data{1,r} = [mean_data{1,r} mean_data_tmp{1,r}];
    end
    
end

%% plot average FR for each trial type

mean_FR = nan(num_bins,3);
sd_FR   = nan(num_bins,3);
se_FR   = nan(num_bins,3);
num_clust= size(mean_data{1,1},2);

for r = 1:3
    mean_FR(:,r) = mean(mean_data{1,r},2);
    sd_FR(:,r)   = std(mean_data{1,r},0,2);
    se_FR(:,r)   = sd_FR(:,r)./sqrt(num_clust);
end

fh = figure;
plotShadedSD(fh,timeframe,mean_FR,se_FR,trial_colors);
legend(trial_types); xlabel('Time (s)'); ylabel('mean FR accross units +- sem');
title('Mean FR for all clusters');
pretty_fig;

%% Calculate PCA space from averaged trials

% conc_mean_data is [M*3xN], M=number of time bins, N=Number of clusters(spikes) [num_bins*3 x numclust]
% it has the mean FR for all spikes in columns, time bins in rows, and vertically concatenated trial types [no_rew; single; jackpot];

%  [ mean(bin1_no_rew_clust1) mean(bin1_no_rew_clust2) ... mean(bin1_no_rew_clustN) ]
%  [ mean(bin2_no_rew_clust1) mean(bin2_no_rew_clust2) ... mean(bin2_no_rew_clustN) ]
%  [       ...                        ...                           ...             ]
%  [ mean(bin1_single_clust1) mean(bin1_single_clust2) ... mean(bin1_single_clustN) ]
%  [       ...                        ...                           ...             ]
%  [ mean(binM_jackpt_clust1) mean(binM_jackpt_clust2) ... mean(binM_jackpt_clustN) ]


conc_mean_data = cell2mat(mean_data');

% calculate PCA axes:
[coeff_m, score_m, latent, tsquared, explained ] = pca(conc_mean_data);

%% Calculate PCA space from individual trials
% TODO?


%% Plot mean trajectories in first 3 PCs

num_PCs = size(coeff_m,2);
mean_scorePC = nan(num_bins,num_PCs,3);
figure; hold on;
for r=1:3
    mean_scorePC(:,:,r) = mean_data{1,r}*coeff_m;
    plot3(mean_scorePC(:,1,r),mean_scorePC(:,2,r),mean_scorePC(:,3,r),trial_colors{r});
end
view(45,45); grid on;
xlabel('PC1');ylabel('PC2');zlabel('PC3');
title('Traj in first 3 PCs');
pretty_fig;
legend(trial_types);

for r=1:3
    % first point is a circle
    plot3(mean_scorePC(1,1,r),mean_scorePC(1,2,r),mean_scorePC(1,3,r),['o-' trial_colors{r}],'MarkerSize',14,'LineWidth',1.5);
    % last point is an x
    plot3(mean_scorePC(end,1,r),mean_scorePC(end,2,r),mean_scorePC(end,3,r),['x-' trial_colors{r}],'MarkerSize',14,'LineWidth',1.5);
end
    
%% Calculate PC projections for indiv trials:
trial_proj = cell(num_sess,3);

% calculate PC proj for every trials
for s=1:num_sess
    for r=1:3
        trial_proj{sessions(s),r} = cell(num_trials(sessions(s),r),1);
        for t=1:num_trials(sessions(s),r)
            trial_proj{sessions(s),r}{t} = trial_data{sessions(s),r}{t}*coeff_m(clust_sidx==sessions(s),:);
        end
    end
end




%% Mahalanobis distance between individual trials and distributions of trials of all types

MD   = cell(3,3); % [NR-NR, NR-S, NR-J; S-NR, S-S, S-J; J-NR, J-S, J-J]
mMD  = cell(3,3); % mean for each bin
mbMD = nan(3,3);  % overall mean across bins
sdMD = cell(3,3); % sd for each bin


for r1 = 1:3
    indiv_trials_proj = vertcat(trial_proj{:,r1});
    num_indiv_trials  = sum(num_trials(:,r1));
    
    for r2 = 1:3
        num_other_trials = sum(num_trials(:,r2));
        
        MD{r1,r2} = nan(num_bins,num_indiv_trials);

        for this_trial=1:num_indiv_trials
            
            if r1==r2
                % same trial_type: calculate MD between each trial and the distrib of the others
                other_trials = setdiff(1:num_indiv_trials,this_trial);
            else
                other_trials = 1:num_other_trials;
            end
            
            other_trials_proj = vertcat(trial_proj{:,r2});
            
            % concatenate other projections in a 3D array
            other_proj3 = cat(3,other_trials_proj{other_trials});
            
            % for each time bin, calculate PCproj of this trial to distrib of proj of other trials
            for b=1:num_bins    
                MD{r1,r2}(b,this_trial) = mahal(indiv_trials_proj{this_trial}(b,1:num_PCs), squeeze(other_proj3(b,1:num_PCs,:))');
            end
        end
        
        mMD{r1,r2}  = mean(MD{r1,r2},2);
        sdMD{r1,r2} = std(MD{r1,r2},0,2);
        mbMD(r1,r2) = mean(mMD{r1,r2});
    end
end

% plot overall mean MD in a confusion matrix
figure;
imagesc(mbMD);
colorbar;
title('confusion matrix of MD'); ylabel('single trials'); xlabel('distributions');
set(gca,'XTick',1:3,'XTickLabel',trial_types);
set(gca,'YTick',1:3,'YTickLabel',trial_types);

% plot MD for each bin
for r = 1:3
    mean_MD = cell2mat(mMD(r,:));
    se_MD   = cell2mat(sdMD(r,:))./sqrt(sum(num_trials(:,r)));
    plotShadedSD(timeframe,mean_MD,se_MD,trial_colors);
    legend(trial_types); xlabel('Time (s)'); ylabel('mean FR accross units');
    pretty_fig;
    title(sprintf('MD between indiv %s trials and other distrib',trial_types{r}));
    xlabel('Time (s)'); ylabel('MD^2 +- sem');
    legend(trial_types);
end






