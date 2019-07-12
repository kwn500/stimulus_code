%% --------STORES INDIVIDUAL PARTICIPANT THRESHOLDS AS .CSV FILE-------- %%
pid = 100;    %enter participants unique id number.
condition_task = 'temporal'; %temporal or spatial task toggle.
threshold_repeats = 3;

% Establish counter variables.
orientation_mean_counter = 0; contrast_mean_counter = 0; RAM_mean_counter = 0;
orientation_SD_counter = 0; contrast_SD_counter = 0; RAM_SD_counter = 0;

% cd into data storage directory.
base_directory = cd('\Users\PS-RMLab\Documents\Kirstie\updated\data_output\experimental_thresholds\');

% cd into individual participants threshold directory.
if pid < 10
    threshold_dir = sprintf('0%d_experimental_thresholds',pid);
else
    threshold_dir = sprintf('%d_experimental_thresholds',pid);
end
cd(threshold_dir);

for scan_number = 1:threshold_repeats
    % load matlab file in this location containing thresholds from pilot scans.
    if pid < 10
        load(sprintf('thresholds_%s_scan_0%d_pp_0%d.mat',condition_task, scan_number, pid));
    else
        load(sprintf('thresholds_%s_scan_0%d_pp_%d.mat',condition_task, scan_number, pid));
    end
    
    % Store a running total of thresholds for later average calculations.
    orientation_mean_counter = orientation_mean_counter + orientationT;
    orientation_SD_counter = orientation_SD_counter + orientationSD;
    contrast_mean_counter = contrast_mean_counter + contrastT;
    contrast_SD_counter = contrast_SD_counter + contrastSD;
    RAM_mean_counter = RAM_mean_counter + RAMT;
    RAM_SD_counter = RAM_SD_counter + RAMSD;
    
    % cd into individual participants threshold directory.
    all_data = {orientationT orientationSD contrastT contrastSD RAMT RAMSD};
    data_table = cell2table(all_data, 'VariableNames', {'Orientation_Threshold',...
        'Orientation_SD', 'Contrast_Threshold',...
        'Contrast_SD', 'RAM_Threshold', 'RAM_SD'});
    
    % save thresholds as .csv file with participant identifier.
    if pid < 10
        writetable(data_table, sprintf('experimental_thresholds_pp_0%d_scan_0%d.csv', pid, scan_number));
    else
        writetable(data_table, sprintf('experimental_thresholds_pp_%d_scan_0%d.csv', pid, scan_number));
    end
end

orientation_mean_counter = orientation_mean_counter / threshold_repeats;
orientation_SD_counter = orientation_SD_counter / threshold_repeats;
contrast_mean_counter = contrast_mean_counter / threshold_repeats;
contrast_SD_counter = contrast_SD_counter / threshold_repeats;
RAM_mean_counter = RAM_mean_counter / threshold_repeats;
RAM_SD_counter = RAM_SD_counter / threshold_repeats;

all_average_data = {orientation_mean_counter; orientation_SD_counter; contrast_mean_counter; contrast_SD_counter;...
    RAM_mean_counter; RAM_SD_counter}';
average_data_table = cell2table(all_average_data, 'VariableNames',...
    {'Orientation_Threshold', 'Orientation_SD', 'Contrast_Threshold',...
    'Contrast_SD', 'RAM_Threshold', 'RAM_SD'});
if pid < 10
    writetable(average_data_table, sprintf('mean_experimental_thresholds_pp_0%d.csv', pid));
else
    writetable(average_data_table, sprintf('mean_experimental_thresholds_pp_%d.csv', pid));
end
%% --------------------------------------------------------------------- %%