function [marked_response] = d_prime_classification(response, type_changes, block_condition)

% if the participant made a response, indicating a change in a particular 
% visual feature:
if isequal(response,'change')
    
    % if this was an within an attend-orientation block:
    if block_condition == 1
        
        % if the target stimulus for this trial included an orientation
        % change:
        if strfind(type_changes,'ORIENTATION') > 0
            
            % mark the participant response as a hit (i.e. they responded
            % to a change in the attended-feature). 
            marked_response = 'Hit';                          
            
            % else, if changes in any of the other features did not occur,
            % and yet a change response was made:                                                 
        elseif isempty(double(strfind(type_changes,'CONTRAST'))) ||...
                isempty(double(strfind(type_changes,'SHAPE')))
            
            % mark this response as a false alarm. 
            marked_response = 'False Alarm';      
        end
    
    % repeat the same process for attend-contrast and attend-shape trials. 
    elseif block_condition == 2
        if strfind(type_changes,'CONTRAST') > 0
            marked_response = 'Hit'; 
        elseif isempty(double(strfind(type_changes,'ORIENTATION'))) ||...
                isempty(double(strfind(type_changes,'SHAPE'))) 
            marked_response = 'False Alarm';
        end
        
    elseif block_condition == 3 
        if strfind(type_changes,'SHAPE') > 0
            marked_response = 'Hit'; 
        elseif isempty(double(strfind(type_changes,'CONTRAST'))) ||...
                isempty(double(strfind(type_changes,'ORIENTATION'))) 
            marked_response = 'False Alarm';
        end
    end
    
% else, if the participant made a response indicating no change in the 
% attended feature:
elseif isequal(response,'no_change')
    
    % if this was within an attend-orientation block:
    if block_condition == 1
        
        % if the orientation of the stimulus did change. 
        if strfind(type_changes, 'ORIENTATION') >= 1
            
            % mark this lack of response as a miss. 
            marked_response = 'Miss';                    
        
        % else, if there was not a change in the orientation of the
        % stimulus:
        elseif isempty(strfind(type_changes, 'ORIENTATION'))
            
            % mark this lack of response as a correct rejection. 
            marked_response = 'Correct Rejection';           
        end
    
    % repeat the same process for attend-contrast and attend-shape blocks. 
    elseif block_condition == 2   
        if strfind(type_changes, 'CONTRAST') >= 1 
            marked_response = 'Miss';
        elseif isempty(strfind(type_changes, 'CONTRAST'));
            marked_response = 'Correct Rejection'; 
        end
    elseif block_condition == 3  
        if strfind(type_changes, 'SHAPE') >= 1
            marked_response = 'Miss';
        elseif isempty(strfind(type_changes, 'SHAPE'));
            marked_response = 'Correct Rejection'; 
        end
    end
end

% if this trial was part of a passive-attention block, or the participant 
% made no response:
if block_condition == 4 || isequal(response, 'no')
    marked_response = 'N/A';  % set the marked response as non-applicable. 
end
end

