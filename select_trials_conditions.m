function I = select_trials_conditions(datatable,varargin)
% usage: I = select_trials_conditions(datatable,varargin)
% this function finds which trials (rows) of datatable meets
% the specified conditions.
%
%   outputs: 
%       I : 1xNtrials logical vector indicating which trials meet
%           the conditions
%       
%   inputs:
%       datatable: mielbaz's mototrak data format
%                   (colomns are "jacpot_bin, force_bin, spikes_bin,
%                   sucessful_bin, time2success and session_index")
%
%       varargin: pairs of optional "condition_names" and "condition_values"
%                 (see default conditions below)
%
%   Example:
%      I = select_trials_conditions(table_debut,'trial_type',[0,1],...
%                                                'success',1,'sessions',[1:10])
%
%      returns logical index indicating which rows of 'table_debut'
%      correspond to successful single and jacpot trials of sessions 1 to
%      10.


% default conditions
conditions = struct(...
    'trial_type'   , [-1,0,1],... %i.e. [no-reward, single, jackpot]
    'success'       , [0,1],... %[failure, success]
    'sessions'      , 0);   % [all]

%update missing params with default values
conditions = parse_input_params(conditions,varargin);

I_type      = ismember(datatable.jackpot_bin,conditions.trial_type);
I_success   = ismember(datatable.successful_bin,conditions.success);
if conditions.sessions
    I_sess  = ismember(datatable.session_index,conditions.sessions);
else
    I_sess = ones(size(datatable,1));
end

I = all([I_type I_success I_sess],2);

