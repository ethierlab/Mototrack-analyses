function SR = rolling_SR(datatable, n_trials)
%%%%%
% usage: [SR, I] = rolling_SR(datatable, varargin)
%
% this function calculates a vector of ongoing sucess rates for the mototrak lever-pulling task
%   the success rate is 
%
%   inputs:
%       datatable: mielbaz's format (jackpot_bin, force_bin, spike_bin, etc.)
%       n_trials: number of previous trials to average success rate
%   
%   outputs:
%       SR : success rate moving average, using linear weights.
%
%   important considerations: 
%       no-rewards are not considered successful trials by default. 
%       the first n_trials-1 of every sessions are assigned NaN as SR
%       the SR associated with a trial is the value reflecting past n_trials only
%           (it does not consider success of current trial)
%   
%    Example :
%       successful_bin : :older trial-> [0 0 0 1 1 1 1 0 1 1 0 0 0] <-newer trial
%       using n_trials = 3
%           ( weights = [ 0.1667    0.3333    0.5000] )
%
%       SR = [NaN NaN NaN 0 0.5 0.83 1 1 0.5 0.6667 0.8333 0.5 0.1667]
%
%
%%%% CE - June 2020 %%%%%%

consider_NR = false; %todo: input argument

SR = nan(size(datatable,1),1);

sessions = unique(datatable.session_index);

w = linspace(0,2/(n_trials+1),n_trials+1);
w = w(2:end); % linear weigths for > importance of recent trials. sum of weigths = 1;

for s = 1:length(sessions)
   
    
    trials_idx = find(datatable.session_index == sessions(s));
    
    if length(trials_idx) < n_trials+1
        continue    %skip session if not enough trials
    end
    
    for t = trials_idx(n_trials+1):trials_idx(end)
        if ~consider_NR
            w_temp = w.*(datatable.jackpot_bin(t-n_trials:t-1)>=0)'; % weight of 0 for NR trials, regardless of outcome
        else
            w_temp = w;
        end
        SR(t) = w_temp * datatable.successful_bin(t-n_trials:t-1);
    end
end
    