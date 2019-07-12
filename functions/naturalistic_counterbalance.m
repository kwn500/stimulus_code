function order = naturalistic_counterbalance(participant)

    % specify the fmri run conditions. 
    conditions = {'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'};

    % seed the random number generator based on the current time to ensure 
    % randomisation is different for each participant.
    rng('shuffle');  
    order = conditions(randperm(length(conditions))); % randomly shuffle the order of conditions. 
    
    % need to add in cd to the relevant participant directory- make their
    % directory.
    cd('C:\Users\Ryan Maloney\Documents\Kirstie\naturalistic\experiment\data');
    
    pp_dir = sprintf('R%d', participant);
    [~,~] = mkdir(pp_dir); cd(pp_dir);
    
    % write the randomised order of conditions to a participant-specific
    % .csv file. 
    dlmwrite(sprintf('R%d_naturalistic_order.csv', participant), order);  
    
    % save the randomised order of conditions to a participant-specific
    % .mat file. 
    save(sprintf('R%d_naturalistic_order.mat', participant), 'order'); 
end
%% --------------------------------------------------------------------- %% 