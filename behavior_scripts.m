
%load('/Users/christianethier/Dropbox (EthierLab)/Michaël/mototrak_tables_pour_Christian/jados_80gr_800ms_noCNO.mat');
%load('/Users/christianethier/Dropbox (Personal)/z-scoring/jasm1_120_800_sansCNO.mat');
load('/Users/christianethier/Dropbox (Personal)/z-scoring/jasm1_80_800_sansCNO.mat');

%% all options and general variables
options_rat     = 'jasm1 80g 800ms sansCNO';
options_FTI     = {'0-1s'     , '300-600ms' , '500-800ms' ,'500-1500ms'};
FTI_values      = {'FTI0_1000', 'FTI300_600', 'FTI500_800','FTI500_1500'};

options_FTI     = {'500-800ms' ,'500-1500ms'};
FTI_values      = {'FTI500_800','FTI500_1500'};

SR_lag          = 10;

% 
% options_FTI     = { 'pre-init', '0-300ms', '300-600ms', '500-800ms' };
% FTI_values      = { 'FTIpre', 'FTI0_300', 'FTI300_600', 'FTI500_800'};
% 
% options_FTI     = { '500-1500ms'};
% FTI_values      = { 'FTI500_1500'};

options_success = {'Successful trials', 'Failed trials', 'All success'};
success_values  = {1,0,[0,1]};

options_reward  = {'No Reward', 'Single', 'Jackpot', 'All rewards'};
Rew_values      = {-1, 0, 1, [-1,0,1]};
rew_colors      = {'k','b','r', [.5 .5 .5]};

options_SR      = {'[0-33]%', '[33-66]%','[66-100]%','All'};
SR_values       = {[0 .33], [.33 .66], [.66 1], [0 1]};
SR_colorSat     = [cellfun(@mean, SR_values(1:end-1)) 1];

% SR_colorSat     = [cellfun(@(x) round(256*mean(x)), SR_values) 1];

pre_init_bin = temps_preinitiation/bin_size;

SR  = rolling_SR(table_debut,SR_lag);
FTI = cellfun(@sum,table_debut.force_bin);
num_trials = size(table_debut,1);


types_id     = [-1,0,1];
num_types    = numel(types_id);
types_labels = {'no reward', 'single', 'jackpot'};
type_colors  = {'k','b','r'};

success_id     = [0,1];
num_success_id = numel(success_id);
success_labels = {'successful','failed'};
success_colors = {'c','m'};
 

% % % 
% % % interval = 1:pre_init_bin-1; %pre-init
% % % FTIpre = cellfun(@(x) sum(x(interval)),table_debut.force_bin);
% % % 
% % % interval = pre_init_bin:(pre_init_bin+300/bin_size); %0-300ms
% % % FTI0_300 = cellfun(@(x) sum(x(interval)),table_debut.force_bin);
% % 
% % interval = pre_init_bin:(pre_init_bin+600/bin_size); %0-1000ms
% % FTI0_1000 = cellfun(@(x) sum(x(interval)),table_debut.force_bin);
% % 
% % interval = (pre_init_bin+300/bin_size):(pre_init_bin+600/bin_size); %300-600ms
% % FTI300_600 = cellfun(@(x) sum(x(interval)),table_debut.force_bin);
% % 
% % interval = (pre_init_bin+500/bin_size):(pre_init_bin+800/bin_size); %500-800ms
% % FTI500_800 = cellfun(@(x) sum(x(interval)),table_debut.force_bin);
% % 
% % interval = (pre_init_bin+500/bin_size):(pre_init_bin+1500/bin_size); %500-1500ms
% % FTI500_1500 = cellfun(@(x) sum(x(interval)),table_debut.force_bin);

%% Overall Success Rate and Force PSTH
I_j = select_trials_conditions(table_debut,'trial_type',1);
I_s = select_trials_conditions(table_debut,'trial_type',0);
I_nr= select_trials_conditions(table_debut,'trial_type',-1);

I_sj = select_trials_conditions(table_debut,'trial_type',1,'success',1);
I_ss = select_trials_conditions(table_debut,'trial_type',0,'success',1);
I_snr= select_trials_conditions(table_debut,'trial_type',-1,'success',1);

R_nr = sum(I_snr)/sum(I_nr);
R_s  = sum(I_ss)/sum(I_s);
R_j  = sum(I_sj)/sum(I_j);

R_values = [R_nr R_s R_j];

types_labels = {'no reward', 'single', 'jackpot'};
type_colors = {'k','b','r'};
figure; hold on;
for i=1:3
    bar(i,R_values(i));
    dat = get(gca,'Children');
    set(dat(1),'FaceColor',type_colors{i});
end

set(gca,'XTick',1:3);
set(gca,'XTickLabel',{'No Reward','Single','Jackpot'});
% set(gca,'XTickLabel',{['No Reward (' R_values(1) '%)'],['Single (' R_values(2) '%)'],['Jackpot  (' R_values(3) '%)']});

ylabel('Success Rate'); title(['Mean Overall Success Rate - ' options_rat]);
pretty_fig;


%% plot overall force PSTH
PSTH_t       = -temps_preinitiation:bin_size:(hold_time_duration+temps_postrecompense-bin_size);
PSTH_numbins = length(PSTH_t);

% big array of force signals for successful trials of each type:
F_j = cell2mat( cellfun(@(x)  x(1:PSTH_numbins)', table_debut.force_bin(I_sj),'UniformOutput',false));
F_s = cell2mat( cellfun(@(x)  x(1:PSTH_numbins)', table_debut.force_bin(I_ss),'UniformOutput',false));
F_nr = cell2mat( cellfun(@(x)  x(1:PSTH_numbins)', table_debut.force_bin(I_snr),'UniformOutput',false));

mean_Fjsnr = [mean(F_nr);mean(F_s);mean(F_j)]'; %mean force for every time bin. order: [nr s j]
IC_Fjsnr  = 1.96*[std(F_s)/sqrt(sum(I_snr)); std(F_s)/sqrt(sum(I_ss));std(F_j)/sqrt(sum(I_sj));]';  %1.96*SE for every time bin

plotShadedSD(PSTH_t,mean_Fjsnr,IC_Fjsnr,type_colors);
legend(types_labels); ylabel('Force (g)'); xlabel('Time (ms)');
ttl =  {['Mean Force - all reward - all success - ' options_rat]};
title(ttl);


%% plot Force PSTH according to ongoing SR
% for Succ_type = 1:length(success_values)

for Succ_type = 1:length(options_success)
    for Rew_type = 1:length(options_reward)
        
        ttl =  {['Mean Force ' options_reward{Rew_type} ' - ' options_success{Succ_type} ' - ' options_rat]};
        
        % table index for reward and success type: (e.g. successful jackpot)
        t_idx = select_trials_conditions(table_debut,'trial_type',Rew_values{Rew_type},'success',success_values{Succ_type});
        
        N = zeros(1,length(SR_values));
        F_mean = [];
        F_IC   = [];
        for SR_type = 1:length(SR_values)  %extract mean force for trials within given range of rolling success rate
            
            t_sidx     = t_idx & SR >= SR_values{SR_type}(1) & SR <= SR_values{SR_type}(2); % subset of index with SR in range
            N(SR_type) = sum(t_sidx);
            F_t        = cell2mat( cellfun(@(x)  x(1:PSTH_numbins)', table_debut.force_bin(t_sidx),'UniformOutput',false));
            F_mean  = [F_mean mean(F_t)'];
            F_IC    = [F_IC 1.96*std(F_t)'/sqrt(N(SR_type))];
        end
        
        % then plot
        figure;
        plot(PSTH_t,F_mean);
        lines = flipud(get(gca,'Children'));
        leg = {};
        for i = 1:length(lines)
            lines(i).Color = rew_colors{Rew_type};
            line_rgb = lines(i).Color;
            lines(i).Color = (1 - (1-line_rgb).*SR_colorSat(i)); % make low SR lighter color
            leg = [leg [options_SR{i} ', N=' num2str(N(i))]];
        end
        pretty_fig; title(ttl); xlabel('Time (ms)'); ylabel('Force (g)');
        legend(leg,'Location','SouthEast');
    end
end


%% Time to success
success_times(table_debut,hold_time_duration/1000);


%% FTI per reward type

for Force_type = 1:length(FTI_values)
    for Succ_type = 1:length(success_values)
        
        ttl = {['Mean FTI ' options_FTI{Force_type} ' - ' options_success{Succ_type} ' - ' options_rat]};
        FTIt = eval(FTI_values{Force_type});
        
        I_j = select_trials_conditions(table_debut,'trial_type',1,'success',success_values{Succ_type});
        I_s = select_trials_conditions(table_debut,'trial_type',0,'success',success_values{Succ_type});
        I_nr= select_trials_conditions(table_debut,'trial_type',-1,'success',success_values{Succ_type});
        
        n_j = sum(I_j);
        n_s = sum(I_s);
        n_nr= sum(I_nr);
        N = [n_nr n_s n_j];
        
        FTI_j = FTIt(I_j);
        FTI_s = FTIt(I_s);
        FTI_nr= FTIt(I_nr);
        
        FTIm = [mean(FTI_nr) mean(FTI_s) mean(FTI_j)];
        FTIe = 1.96*[std(FTI_nr)/sqrt(n_nr) std(FTI_s)/sqrt(n_s) std(FTI_j)/sqrt(n_j)];
        
        % bar plot:
        figure;
        for i=1:num_types
            barwitherr(FTIe(i),i,FTIm(i));
            hold on;
            dat = get(gca,'Children');
            bars = dat(2);
            set(bars,'FaceColor',type_colors{i});
        end
        
        %labels and formating:
        set(gca,'XTick',1:num_types);
        xtl = get(gca,'XTickLabel');
        for ttype = 1:num_types
            xtl{ttype} = sprintf([types_labels{ttype} ' (n=%d)'],N(ttype));
        end
        set(gca,'XTickLabel',xtl);
        ylabel('FTI (±1.96*SE)'); title(ttl);
        pretty_fig;
        
    end
end

%% FTI per Success Type (regardless of reward types)

for Force_type = 1:length(FTI_values)
    
    ttl = {['Mean FTI - ' options_FTI{Force_type} ' - ' options_rat ]};
    FTIt = eval(FTI_values{Force_type});
    
    I_su = select_trials_conditions(table_debut,'success',1);
    I_fa = select_trials_conditions(table_debut,'success',0);
    
    n_su = sum(I_su);
    n_fa = sum(I_fa);
    N= [n_su n_fa];
    
    FTI_su = FTIt(I_su);
    FTI_fa = FTIt(I_fa);
    
    FTIm = [mean(FTI_su) mean(FTI_fa)];
    FTIe = 1.96*[std(FTI_su)/sqrt(n_su) std(FTI_fa)/sqrt(n_fa)];
    
    % bar plot:
    figure;
    for i=1:num_success_id
        barwitherr(FTIe(i),i,FTIm(i));
        hold on;
        dat = get(gca,'Children');
        bars = dat(2);
        set(bars,'FaceColor',success_colors{i});
    end
    
    %labels and formating:
    set(gca,'XTick',1:num_success_id);
    xtl = get(gca,'XTickLabel');
    for ttype = 1:num_success_id
        xtl{ttype} = sprintf([success_labels{ttype} ' (n=%d)'],N(ttype));
    end
    set(gca,'XTickLabel',xtl);
    ylabel('FTI (±1.96*SE)');
    title(ttl);
    pretty_fig;
    
end

%% FTI vs SR

for Force_type = 1:length(FTI_values)
    for Rew_type = 1:length(options_reward)
        
        FTIt = eval(FTI_values{Force_type});
        
        %smooth FTIt (moving average filter function:)
        windowSize = 10; 
        b = (1/windowSize)*ones(1,windowSize);
        a = 1;
        FTIts = filtfilt(b,a,FTIt);
        
        I = select_trials_conditions(table_debut,'trial_type',Rew_values{Rew_type});
                 N = sum(I);
       
        figure;
        subplot(2,1,1);
        plot(1:N,FTIt(I),'r');
        hold on;
        plot(1:N,FTIts(I),'k');
        pretty_fig;
        ylabel('FTI (au)');
        title([options_FTI{Force_type} ' - ' options_reward{Rew_type} ' - ' options_rat]);
        
        SRs = filtfilt(b,a,SR);
        
        subplot(2,1,2);
        plot(1:N,SR(I),'b'); hold on;
        plot(1:N,SRs(I),'k'); pretty_fig;
        title(['Rolling SR for 5 past trials (all trials, linear weights) - at times of '  options_reward{Rew_type} ' trials - ' options_rat]);
        xlabel([options_reward{Rew_type} 'Trials']);
        
        I = I & ~isnan(SRs);
        figure;
        plot(SRs(I),FTIts(I),'.'); hold on;
        [Fo, Go] = fit(SRs(I),FTIts(I),'poly1'); 
        plot(Fo);
        ylabel(['FTI ' options_FTI{Force_type}]);xlabel('Ongoing SR');
        title([options_FTI{Force_type} ' for ' options_reward{Rew_type} ' trials vs SR']);
        legend(['LinFit, R2 = ' num2str(Go.rsquare)]); pretty_fig;
        
    end
    
end



