function [sr_all, sr_sess, N_sess] = success_rates(bindata, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% returns success rats for each trial type ("no_reward", "single", or "jackpot")
% usage: sr = success_rates(bindata,[plot_flag])
%
% input args:
%       bindata         : binned data extracted with Michaël's mysterious function
%       plot_flag       : [optional] default = 1. Set to 0 to avoid also producing histograms and time-performance plots
%
% output args:
%       sr_all          : 1x3 array of success rates for each trial type, in this order: {no_reward, single, jackpot}
%       sr_sess         : cell array of success rates for each trial type, for each session, in this order: {no_reward, single, jackpot}
%
%%%%%% EthierLab - 03/2019 - CE %%%%%%%%%%%%%%%%%


plot_flag = 1;
if nargin>1
    plot_flag = varargin{1};
end

types_id     = [-1,0,1];
num_types    = numel(types_id);
types_labels = {'no reward', 'single', 'jackpot'};
type_colors  = {'k','b','r'};

folder_idx = strfind(bindata.data_folder,'\');
folder_str = bindata.data_folder(folder_idx(end)+1:end);
data_info = [bindata.animal_name ' (HT=' num2str(bindata.hold_time*1000) 'ms, FT=' num2str(bindata.threshold) 'g, ' folder_str ')'];
data_info = strrep(data_info,'_','\_');

sess = unique(bindata.bin_initial_tout.session_index);
num_sess = numel(sess);

% to store success rate in that order: {no_reward, single, jackpot}
sr_all  = nan(1,num_types);
sr_sess = cell(1,num_types);
N       = nan(1,num_types);
N_sess  = cell(1,num_types);

% get overall sr for each trial type
for ttype = 1:num_types
    trial_select = bindata.bin_initial_tout.jackpot_bin == types_id(ttype);
    num_t = sum(trial_select);
    N(ttype) = num_t;
    num_suc = sum(bindata.bin_initial_tout.successful_bin(trial_select));
    sr_all(1,ttype) = num_suc/num_t;
end

% get sr for each session
for s = 1:num_sess
    % get sr for each trial type
    for ttype = 1:num_types
        trial_select = bindata.bin_initial_tout.jackpot_bin == types_id(ttype) & bindata.bin_initial_tout.session_index == sess(s);
        num_t = sum(trial_select);
        N_sess{1,ttype} = [N_sess{1,ttype}; num_t];
        num_suc = sum(bindata.bin_initial_tout.successful_bin(trial_select));
        sr_sess{1,ttype} = [sr_sess{1,ttype}; num_suc/num_t];
    end
end

% plot bar plot
figure; hold on;
% error bar: 1.96*SE accross sessions
Err = cellfun(@(x) 1.96*std(x(~isnan(x)))/sqrt(length(x(~isnan(x)))),sr_sess);
for ttype = 1:num_types
    barwitherr(Err(ttype), ttype, sr_all(ttype));
    hold on;
    dat = get(gca,'Children');
    bars = dat(2);
    set(bars,'FaceColor',type_colors{ttype});
end
set(gca,'XTick',1:num_types);
xtl = get(gca,'XTickLabel');
for ttype = 1:num_types
    xtl{ttype} = sprintf([types_labels{ttype} ' (n=%d)'],N(ttype));
end
set(gca,'XTickLabel',xtl);

ylabel('Success rate');
title(sprintf('Success rate across %d sessions (1.96*SE)\n%s',num_sess,data_info));
pretty_fig;