function FTI = FTI_trials(datatable)
% This function returns the force-time integral for each trial
% taking into input argument mielbaz's mototrak data tables.

    FTI = cellfun(@sum,datatable.force_bin();