
load('/Users/christianethier/Dropbox (EthierLab)/Michaël/mototrak_tables_pour_Christian/jados_80gr_800ms_noCNO.mat');

SR  = rolling_SR(table_debut,5);
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

interval = 30:50; %first second only
FTI1 = cellfun(@(x) sum(x(interval)),datatable.force_bin);

interval = 40:60;
FTI15 = cellfun(@(x) sum(x(interval)),datatable.force_bin);

interval = 50;
FTIend = cellfun(@(x) sum(x(interval:end)),datatable.force_bin);

%% FTI per trial type
%
FTIt = FTIend; % <------ Choose!!! (FTI, FTI1, FTI15...) --------
% % %
% I_j = select_trials_conditions(table_debut,'trial_type',1);
% I_s = select_trials_conditions(table_debut,'trial_type',0);
% I_nr= select_trials_conditions(table_debut,'trial_type',-1);

I_j = select_trials_conditions(table_debut,'trial_type',1,'success',1);
I_s = select_trials_conditions(table_debut,'trial_type',0,'success',1);
I_nr= select_trials_conditions(table_debut,'trial_type',-1,'success',1);

ttl = {'Mean FTI after 1st sec - Failed trials- Jados 80g,800ms'};

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

%% FTI per Success

FTIs = FTI15; % <------ Choose!!! (FTI, FTI1, FTI15...) --------
ttl = {'Mean FTI - 500-1500ms - Jados 80g,800ms,noCNO'};

I_su = select_trials_conditions(table_debut,'success',1);
I_fa = select_trials_conditions(table_debut,'success',0);
n_su = sum(I_su);
n_fa = sum(I_fa);
N= [n_su n_fa];

FTI_su = FTIs(I_su);
FTI_fa = FTIs(I_fa);

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

%% Overall Success Rate
I_j = select_trials_conditions(table_debut,'trial_type',1);
I_s = select_trials_conditions(table_debut,'trial_type',0);
I_nr= select_trials_conditions(table_debut,'trial_type',-1);

I_sj = select_trials_conditions(table_debut,'trial_type',1,'success',1);
I_ss = select_trials_conditions(table_debut,'trial_type',0,'success',1);
I_snr= select_trials_conditions(table_debut,'trial_type',-1,'success',1);

R_nr = sum(I_snr)/sum(I_nr);
R_s  = sum(I_ss)/sum(I_s);
R_j  = sum(I_sj)/sum(I_j);

values = [R_nr R_s R_j];

types_labels = {'no reward', 'single', 'jackpot'};
type_colors = {'k','b','r'};
figure; hold on;
for i=1:3
    bar(i,values(i));
    dat = get(gca,'Children');
    set(dat(1),'FaceColor',type_colors{i});
end

set(gca,'XTick',1:3);
set(gca,'XTickLabel',{'No Reward','Single','Jackpot'});

ylabel('Success Rate'); title('Mean Overall Success Rate - Jados 80g,800ms');
pretty_fig;

%% FTI vs SR

figure;
subplot(3,1,1);
plot(1:num_trials,FTI,'r'); pretty_fig;
ylabel('FTI (au)');
title('Force-Time Integral for whole trial - Jados 80, 800ms, no CNO');

subplot(3,1,2);
plot(1:num_trials,FTI15,'m'); pretty_fig;
ylabel('FTI (au)');
title('Force-Time Integral 500-1500ms - Jados 80, 800ms, no CNO');

subplot(3,1,3);
plot(1:num_trials,SR,'b'); pretty_fig;
title('Rolling Success Rate, 5 past trials, linear weights - Jados 80, 800ms, no CNO');
xlabel('Trials');

idx = ~isnan(SR);
figure;
plot(SR(idx),FTI1(idx),'.');pretty_fig;
ylabel('FTI1');xlabel('Ongoing SR'); title('FTI1 vs SR');

figure;
plot(SR(idx),FTI15(idx),'.');pretty_fig;
ylabel('FTI15');xlabel('Ongoing SR'); title('FTI15 vs SR');





