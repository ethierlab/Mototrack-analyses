function t2s = success_times(bindata,hold_time, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% returns trial times for successful trials of each type ("no_reward", "single", or "jackpot")
% usage: t2s = success_times(bindata,[plot_flag])
%
% input args:
%       bindata         : binned data extracted with Michaël's mysterious function
%       hold_time       : required hold time, in seconds
%       plot_flag       : [optional] default = 1. Set to 0 to avoid also producing histograms and edcf plots
%
%
% output args:
%       t2s             : cell array of time to success for each trial type, in this order: {no_reward, single, jackpot}
%
%%%%%% EthierLab - 03/2019 - CE %%%%%%%%%%%%%%%%%

plot_flag = 1;
if nargin>2
    plot_flag = varargin{1};
end

types_id     = [-1,0,1];
num_types    = numel(types_id);
types_labels = {'no reward', 'single', 'jackpot'};
type_colors = {'k','b','r'};

% folder_idx = strfind(bindata.data_folder,'\');
% folder_str = bindata.data_folder(folder_idx(end)+1:end);
% data_info = [bindata.animal_name ' (HT=' num2str(bindata.hold_time*1000) 'ms, FT=' num2str(bindata.threshold) 'g, ' folder_str ')'];
% data_info = strrep(data_info,'_','\_');

% to store time 2 success in that order: {no_reward, single, jackpot}
t2s = cell(1,num_types);

% get t2s for each trial type
for ttype = 1:num_types
    trial_select = bindata.jackpot_bin == types_id(ttype) & bindata.successful_bin;
    t2s{ttype} = bindata.time2success(trial_select);
    
    %data in array is actually time before a successful hold time, so have to add hold time
    t2s{ttype} = t2s{ttype} + hold_time;
end

% number of successes:
N = cellfun(@length,t2s);

% plot if wanted
if plot_flag

    edges = 0:0.1:6;
         
%     % plot distribution histograms on top of each other
%     figure;
%     hold on;
%     for i=1:num_types
%         histogram(t2s{i},edges,'Normalization', 'probability','FaceColor',type_colors{i},'FaceAlpha',0.3);
%     end
%     xlabel('Reward time after trial onset (s)'); ylabel('proportion of successful trials');
%     title(sprintf('Time to success per trial type\n%s',data_info))
%     legend(sprintf('no reward (n=%d)',Ns(1)),sprintf('single (n=%d)',Ns(2)),sprintf('jackpot (n=%d)',Ns(3)));
%     pretty_fig;
    
    % plot distribution histograms in 3 different subplots
    figure;
    for i=1:num_types
        
        subplot(3,1,i);
        
        histogram(t2s{i},edges,'Normalization', 'probability','FaceColor',type_colors{i});
        ylabel('proportion of successful trials');
        legend(sprintf([types_labels{i} ' (n=%d)'],N(i)));
        if i==1
            title('Time to success per trial type');
        end
        if i==3
            xlabel('Reward time after trial onset (s)');
        end
        pretty_fig;
    end
    
    pretty_fig;
    
    
    
    % plot cummulative distributions
    figure; 
    hold on;
    lh = [];
    
    for i = 1:num_types
       [F,X,FLO,FUP] = ecdf(t2s{i});
       
       %plot confidence intervals:
       ytop = FUP(2:end-1);
       ybot = FLO(2:end-1);
       yarea = [ytop; ybot(end:-1:1)];
       numpts = length(X);
       xIdx  = [2:numpts-1 numpts-1:-1:2];
       xarea = X(xIdx);
       patch(xarea,yarea,type_colors{i},'LineStyle','none','FaceAlpha',0.1);
       
       % plot lines
       lh = [lh, plot(X,F,type_colors{i})];
    end
    xlim([0 6]);
    pretty_fig;
    legend(lh,sprintf('no reward (n=%d)',N(1)),sprintf('single (n=%d)',N(2)),sprintf('jackpot (n=%d)',N(3)),'Location','NorthWest');
    title('Cummulative distribution of time to success (with CI)');
    ylabel('Cummulative fraction of successful trials'); xlabel('Reward time after trial onset (s)');
    
    
    % plot bar plot
    
    figure;
    % error bar: 1.96*SE
    Err = cellfun(@(x) 1.96*std(x(~isnan(x)))/sqrt(length(x(~isnan(x)))),t2s);
    
    for ttype = 1:num_types
         barwitherr(Err(ttype), ttype, mean(t2s{ttype}));
%         barwitherr(Err(ttype),mean(t2s{ttype}));
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
    
    ylabel('Time to success');
    title('Time to success for all trials (1.96*SE)');
    pretty_fig;
    
end