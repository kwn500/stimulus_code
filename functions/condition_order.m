function order = condition_order(toggle, pp, run)

% calculate the order of runs for minimum motion, contrast detection and
% chromatic threshold staircases and writes these orders to .mat files. 

rng('shuffle')
start_dir = pwd;

if isequal(toggle, 'minimum_motion')
    conditions = {'RG','BY'};
    
    % randomise order of conditions and write to .mat file. 
    conditions = conditions(randperm(length(conditions)));
    order{1} = conditions{1}; order{2} = conditions{2};
    
    if pp < 10 
        cd('C:\Users\PS-RMLab\Documents\Kirstie\colour\psychophysics_data\minimum_motion\'); 
        pp_dir = sprintf('pp_0%d',pp);
        mkdir(pp_dir); cd(pp_dir);
        save(strcat(sprintf('pp_0%d', pp),'_minimum_motion_order'), 'order');
    else
        cd('C:\Users\PS-RMLab\Documents\Kirstie\colour\psychophysics_data\minimum_motion\');
        pp_dir = sprintf('pp_%d',pp);
        mkdir(pp_dir); cd(pp_dir);
        save(strcat(sprintf('pp_%d', pp),'_minimum_motion_order'), 'order');
    end
    cd(start_dir);
    
elseif isequal(toggle, 'contrast_detection')
    conditions = {'Luminance', 'RG', 'BY'};
    
    % repeat the same process for the remaining conditions. 
    conditions = conditions(randperm(length(conditions)));
    order{1} = conditions{1}; order{2} = conditions{2}; order{3} = conditions{3};
    
    if pp < 10 
        cd('C:\Users\PS-RMLab\Documents\Kirstie\colour\psychophysics_data\contrast_detection\');
        pp_dir = sprintf('pp_0%d',pp);
        mkdir(pp_dir); cd(pp_dir);
        save(strcat(sprintf('pp_%0d', pp), sprintf('_run_0%d', run),'_contrast_detection_order'), 'order');
    else
        cd('C:\Users\PS-RMLab\Documents\Kirstie\colour\psychophysics_data\contrast_detection\');
        pp_dir = sprintf('pp_%d',pp);
        mkdir(pp_dir); cd(pp_dir);
        cd(sprintf('run_0%d', run));
        save(strcat(sprintf('pp_%d', pp),sprintf('_run_0%d', run),'_contrast_detection_order'), 'order');
    end
    cd(start_dir);
    
elseif isequal(toggle, 'feature_thresholds')
    
    conditions = {'Orientation Luminance', 'Contrast Luminance', 'Shape Luminance',...
                    'Orientation RG', 'Contrast RG', 'Shape RG',...
                    'Orientation BY', 'Contrast BY', 'Shape BY',};
                
    conditions = conditions(randperm(length(conditions)));
    order{1} = conditions{1}; order{2} = conditions{2}; order{3} = conditions{3};
    order{4} = conditions{4}; order{5} = conditions{5}; order{6} = conditions{6};
    order{7} = conditions{7}; order{8} = conditions{8}; order{9} = conditions{9};

    if pp < 10 
        cd('C:\Users\PS-RMLab\Documents\Kirstie\colour\psychophysics_data\staircase_thresholds\');
        pp_dir = sprintf('pp_0%d',pp);
        mkdir(pp_dir); cd(pp_dir);
        run_dir = sprintf('run_0%d', run);
        save(strcat(sprintf('pp_%0d', pp), sprintf('_run_0%d', run),'_staircase_thresholds_order'), 'order');
    else
        cd('C:\Users\PS-RMLab\Documents\Kirstie\colour\psychophysics_data\staircase_thresholds\');
        pp_dir = sprintf('pp_%d',pp);
        mkdir(pp_dir); cd(pp_dir);
        save(strcat(sprintf('pp_%d', pp),sprintf('_run_0%d', run),'_staircase_thresholds_order'), 'order');
    end
    cd(start_dir);
end
%% --------------------------------------------------------------------- %%  