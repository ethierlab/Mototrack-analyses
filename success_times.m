function t2s = success_times(bindata,varargin)
% returns trial times for successful trials of each type ("no_reward", "single", or "jackpot")

plot_flag = 0;
types_id = [-1,0,1];
types_labels = {'no reward', 'single', 'jackpot'};
type_colors = {'k','b','r'};

if nargin>1
    plot_flag = varargin{1};
    %     if nargin >2
    %         types_id = varargin{2}; %pass types_id as an input arg? not for now...
    %     end
end
num_types = numel(types_id);

% to store time 2 success in that order: {no_reward, single, jackpot}
t2s = cell(1,num_types);

% get t2s for each trial type
for ttype = 1:num_types
    trial_select = bindata.bin_initial_tout.jackpot_bin == types_id(ttype) & bindata.bin_initial_tout.successful_bin;
    t2s{ttype} = bindata.bin_initial_tout.time2success(trial_select);
    
    %data in array is actually time before a successful hold time, so have to add hold time
    t2s{ttype} = t2s{ttype} + bindata.hold_time;
end

% number of successes:
Ns = cellfun(@length,t2s);

% plot if wanted
if plot_flag
    
    % plot distribution histograms on top of each other
    figure;
    hold on;
    edges = [0:0.1:6];
    for i=1:num_types
        histogram(t2s{i},edges,'Normalization', 'probability','FaceColor',type_colors{i},'FaceAlpha',0.3);
    end
    xlabel('Reward time after trial onset (s)'); ylabel('proportion of successful trials');
    folder_idx = strfind(bindata.data_folder,'\');
    folder_str = bindata.data_folder(folder_idx(end)+1:end);
    data_info = [bindata.animal_name ' (HT=' num2str(bindata.hold_time*1000) 'ms, FT=' num2str(bindata.threshold) 'g, ' folder_str ')'];
    data_info = strrep(data_info,'_','\_');
    title(sprintf('Time to success per trial type\n%s',data_info))
    legend(sprintf('no reward (n=%d)',Ns(1)),sprintf('single (n=%d)',Ns(2)),sprintf('jackpot (n=%d)',Ns(3)));
    pretty_fig;
    
    % plot distribution histograms in 3 different subplots
    figure;
    for i=1:num_types;
        
        subplot(3,1,i);
        
        histogram(t2s{i},edges,'Normalization', 'probability','FaceColor',type_colors{i});
        ylabel('proportion of successful trials');
        legend(sprintf([types_labels{i} ' (n=%d)'],Ns(i)));
        if i==1
            title(sprintf('Time to success per trial type\n%s',data_info));
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
    
    for i = 1:num_types;
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
    legend(lh,sprintf('no reward (n=%d)',Ns(1)),sprintf('single (n=%d)',Ns(2)),sprintf('jackpot (n=%d)',Ns(3)),'Location','NorthWest');
    title(sprintf('Cummulative distribution of time to success (with CI)\n%s',data_info));
    ylabel('ECDF'); xlabel('Reward time after trial onset (s)');
    
end