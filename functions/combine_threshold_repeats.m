function combine_threshold_repeats(participant,toggle)

% toggle is either 'features', or 'contrast'. 

% ----------- CALCULATE AVERAGED FEATURE-SPECIFIC THRESHOLDS ------------ %

% loads in participant, feature, colour and scan-specifc .mat threshold
% files and create averaged thresholds across scan repetitions. Writes out
% data in participant, feature, colour and scan-specific .mat files for
% loading into fmri experiment, and saves all data to a .csv file.

% KWN (12\02\2018)

% ----------------------------------------------------------------------- %
%                             INITIAL SETUP                               %
% ----------------------------------------------------------------------- %
% clear; clc; % general housekeeping.

start_directory = pwd;
% participant = 1;
wpid = sprintf('%d', participant); % specify participant to process.

% format participant number with leading 0.
if participant < 10
    wpid = sprintf('0%d', participant);
end

if isequal(toggle,'features')
    feature_tasks = {'orientation', 'contrast', 'shape'}; % specify feature task.
    
    % specify and move into data storage directory.
    threshold_directory = 'C:\Users\PS-RMLab\Documents\Kirstie\colour\psychophysics_data\staircase_thresholds\'; cd(threshold_directory);
    
    % specify and move into the participant-specific data storage directory.
    cd(sprintf('pp_%s', wpid));
    
    % find the number of run repetitions and store the names of the directories.
    run_numbers = dir('run*'); run_numbers = {run_numbers.name};
    
    % ----------------------------------------------------------------------- %
    %                  EXTRACT FEATURE X COLOUR THRESHOLDS                    %
    % ----------------------------------------------------------------------- %
    
    for i = 1:length(run_numbers) % for each of the runs in turn:
        
        cd(run_numbers{i}); % move into the run-specific directory.
        
        for x = 1:length(feature_tasks) % for each of the feature tasks in turn:
            
            cd(feature_tasks{x}); % move into the feature-specific directory.
            
            % load in the blue-yellow feature and scan-specific .mat threshold
            % directory.
            load(sprintf('by_%s_feature_thresholds_pp_%s_run_%s.mat', feature_tasks{x},...
                wpid, run_numbers{i}(5:6)));
            
            if x == 1 % if the feature_task is orientation:
                
                % store the colour-specific orientation sd and mean in a
                % scan-specific position.
                byo_sd{i} = orientationSD; byo_mean{i} = orientationT;
                
                % perform a similar process for contrast and shape feature tasks.
            elseif x == 2
                byc_sd{i} = contrastSD; byc_mean{i} = contrastT;
            elseif x == 3
                bys_sd{i} = shapeSD; bys_mean{i} = shapeT;
            end
            
            % clear variables to ensure they're refreshed for each colour condition.
            clear orientationSD orientationT contrastSD contrastT RAMSD RAMT
            
            % repeat a similar process for red-green and luminance threshold values.
            load(sprintf('rg_%s_feature_thresholds_pp_%s_run_%s.mat', feature_tasks{x},...
                wpid, run_numbers{i}(5:6)));
            
            if x == 1
                rgo_sd{i} = orientationSD; rgo_mean{i} = orientationT;
            elseif x == 2
                rgc_sd{i} = contrastSD; rgc_mean{i} = contrastT;
            elseif x == 3
                rgs_sd{i} = shapeSD; rgs_mean{i} = shapeT;
            end
            
            clear orientationSD orientationT contrastSD contrastT RAMSD RAMT
            
            load(sprintf('luminance_%s_feature_thresholds_pp_%s_run_%s.mat', feature_tasks{x},...
                wpid, run_numbers{i}(5:6)));
            
            if x == 1
                lumo_sd{i} = orientationSD; lumo_mean{i} = orientationT;
            elseif x == 2
                lumc_sd{i} = contrastSD; lumc_mean{i} = contrastT;
            elseif x == 3
                lums_sd{i} = shapeSD; lums_mean{i} = shapeT;
            end
            
            clear orientationSD orientationT contrastSD contrastT RAMSD RAMT
            
            % move into the original scan number directory.
            cd(threshold_directory); cd(sprintf('pp_%s', wpid)); cd(run_numbers{i});
        end
        % move into the original participant directory.
        cd(threshold_directory); cd(sprintf('pp_%s', wpid));
    end
    
    % calculate the mean of the thresholds and standard deviations for each
    % feature task and colour condition.
    rgo_tmean = mean(cell2mat(rgo_mean)); rgo_sdmean = mean(cell2mat(rgo_sd));
    rgc_tmean = mean(cell2mat(rgc_mean)); rgc_sdmean = mean(cell2mat(rgc_sd));
    rgs_tmean = mean(cell2mat(rgs_mean)); rgs_sdmean = mean(cell2mat(rgs_sd));
    
    byo_tmean = mean(cell2mat(byo_mean)); byo_sdmean = mean(cell2mat(byo_sd));
    byc_tmean = mean(cell2mat(byc_mean)); byc_sdmean = mean(cell2mat(byc_sd));
    bys_tmean = mean(cell2mat(bys_mean)); bys_sdmean = mean(cell2mat(bys_sd));
    
    lumo_tmean = mean(cell2mat(lumo_mean)); lumo_sdmean = mean(cell2mat(lumo_sd));
    lumc_tmean = mean(cell2mat(lumc_mean)); lumc_sdmean = mean(cell2mat(lumc_sd));
    lums_tmean = mean(cell2mat(lums_mean)); lums_sdmean = mean(cell2mat(lums_sd));
    
    % create and move into an 'averaged thresholds' data storage directory.
    output_dir = 'averaged_thresholds'; [~,msg] = mkdir(output_dir); cd(output_dir);
    
    % save the mean feature and colour-specific thresholds and standard deviations to
    % associated .mat files for loading into later fmri experiments.
    save(sprintf('rg_averaged_feature_thresholds_pp_%s.mat',wpid), 'rgo_tmean', 'rgo_sdmean', 'rgc_tmean', 'rgc_sdmean',...
        'rgs_tmean', 'rgs_sdmean');
    save(sprintf('by_averaged_feature_thresholds_pp_%s.mat',wpid), 'byo_tmean', 'byo_sdmean', 'byc_tmean', 'byc_sdmean',...
        'bys_tmean', 'bys_sdmean');
    save(sprintf('lum_averaged_feature_thresholds_pp_%s.mat',wpid), 'lumo_tmean', 'lumo_sdmean', 'lumc_tmean', 'lumc_sdmean',...
        'lums_tmean', 'lums_sdmean');
    
    % combine all data into a singular cell array.
    all_data = {rgo_tmean, rgo_sdmean, rgc_tmean, rgc_sdmean, rgs_tmean, rgs_sdmean, byo_tmean, byo_sdmean, byc_tmean,...
        byc_sdmean, bys_tmean, bys_sdmean, lumo_tmean, lumo_sdmean, lumc_tmean, lumc_sdmean, lums_tmean, lums_sdmean};
    
    % convert all data to a cell with associated headings.
    data_table = cell2table(all_data, 'VariableNames', {'rgo_tmean' 'rgo_sdmean', 'rgc_tmean' 'rgc_sdmean' 'rgs_tmean'...
        'rgs_sdmean' 'byo_tmean' 'byo_sdmean', 'byc_tmean' 'byc_sdmean' 'bys_tmean' 'bys_sdmean' 'lumo_tmean' 'lumo_sdmean'...
        'lumc_tmean' 'lumc_sdmean' 'lums_tmean' 'lums_sdmean'});
    
    % write the combined averaged threshold data to a .csv file.
    writetable(data_table, sprintf('averaged_feature_thresholds_pp_%s.csv', wpid));
    % ----------------------------------------------------------------------- %
    
elseif isequal(toggle, 'contrast')
    % specify and move into data storage directory.
    detection_directory = 'C:\Users\PS-RMLab\Documents\Kirstie\colour\psychophysics_data\contrast_detection\'; cd(detection_directory);
    
    % specify and move into the participant-specific data storage directory.
    cd(sprintf('pp_%s', wpid));
    
    % find the number of run repetitions and store the names of the directories.
    run_numbers = dir('run*'); run_numbers = {run_numbers.name};
    
    % ----------------------------------------------------------------------- %
    %                  EXTRACT FEATURE X COLOUR THRESHOLDS                    %
    % ----------------------------------------------------------------------- %
    
    colour_conditions = {'rg', 'by', 'luminance'};
    
    for i = 1:length(run_numbers) % for each of the runs in turn:
        
        cd(run_numbers{i}); % move into the run-specific directory.
        
        for x = 1:length(colour_conditions)
            
            % load in the blue-yellow feature and scan-specific .mat threshold
            % directory.
            load(sprintf('%s_contrast_detection_thresholds_pp_%s_run_%s.mat', colour_conditions{x},...
                wpid, run_numbers{i}(5:6)));
            
            if x == 1 % if the feature_task is orientation:
                
                % store the colour-specific orientation sd and mean in a
                % scan-specific position.
                rgcontrast_sd{i} = contrastdetectionSD;
                rgcontrast_mean{i} = contrastdetectionT;
                
                % perform a similar process for contrast and shape feature tasks.
            elseif x == 2
                bycontrast_sd{i} = contrastdetectionSD;
                bycontrast_mean{i} = contrastdetectionT;
            elseif x == 3
                lumcontrast_sd{i} = contrastdetectionSD;
                lumcontrast_mean{i} = contrastdetectionT;
            end
            
            % move into the original scan number directory.
            cd(detection_directory); cd(sprintf('pp_%s', wpid)); cd(run_numbers{i});
        end
        % move into the original participant directory.
        cd(detection_directory); cd(sprintf('pp_%s', wpid));
    end
    
    % calculate the mean of the thresholds and standard deviations for each
    % feature task and colour condition.
    
    rg_min_contrast_mean = mean(cell2mat(rgcontrast_mean));
    rg_min_contrast_sd = mean(cell2mat(rgcontrast_sd));
    
    by_min_contrast_mean = mean(cell2mat(bycontrast_mean));
    by_min_contrast_sd = mean(cell2mat(bycontrast_sd));
    
    lum_min_contrast_mean = mean(cell2mat(lumcontrast_mean));
    lum_min_contrast_sd = mean(cell2mat(lumcontrast_sd));
    
    % create and move into an 'averaged thresholds' data storage directory.
    output_dir = 'averaged_thresholds'; [~,msg] = mkdir(output_dir); cd(output_dir);
    
    % save the mean feature and colour-specific thresholds and standard deviations to
    % associated .mat files for loading into later fmri experiments.
    save(sprintf('rg_averaged_contrast_detection_threshold_pp_%s.mat', wpid),...
        'rg_min_contrast_mean', 'rg_min_contrast_sd');
    save(sprintf('by_averaged_contrast_detection_threshold_pp_%s.mat', wpid),...
        'by_min_contrast_mean', 'by_min_contrast_sd');
    save(sprintf('lum_averaged_contrast_detection_threshold_pp_%s.mat', wpid),...
        'lum_min_contrast_mean', 'lum_min_contrast_sd');
    
    % combine all data into a singular cell array.
    all_data = {rg_min_contrast_mean, rg_min_contrast_sd, by_min_contrast_mean,...
        by_min_contrast_sd, lum_min_contrast_mean, lum_min_contrast_sd};
    
    % convert all data to a cell with associated headings.
    data_table = cell2table(all_data, 'VariableNames', {'rg_tmean', 'rg_sdmean',...
        'by_tmean', 'by_sdmean', 'lum_tmean', 'lum_sdmean'});
    
    % write the combined averaged threshold data to a .csv file.
    writetable(data_table, sprintf('averaged_contrast_detection_thresholds_pp_%s.csv', wpid));
    % ----------------------------------------------------------------------- %
    cd(start_directory)
end

%% --------------------------------------------------------------------- %%