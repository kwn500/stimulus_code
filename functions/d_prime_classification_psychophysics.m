function [marked_response] = d_prime_classification_psychophysics(response, type_changes, specific_condition, block_condition)

% takes in attention condition, stimulus changes and participant response
% and classifies each response as hit, miss, false alarm or correct
% rejection. 

% if response was made in the non-selective condition: 
if block_condition == 1  
    
    % and if the participant indicated a change: 
    if isequal(response,'change')   
        
            %if the stimulus included a change in any of the features:
            if isempty(strfind(type_changes,'ORIENTATION')) == 0  ||...
                    isempty(strfind(type_changes,'CONTRAST')) == 0||...
                    isempty(strfind(type_changes,'SHAPE')) == 0
                marked_response = 'Hit';   %mark this response as a hit. 
                
            % else if the participant made a response, but no feature 
            % changed, mark as false alarm.     
            elseif isempty(strfind(type_changes,'ORIENTATION')) &&...
                    isempty(strfind(type_changes,'CONTRAST')) && isempty(strfind(type_changes, 'SHAPE')); 
                marked_response = 'False Alarm'; 
            end
            
            % if the participant did not indicate a change:
    elseif isequal(response, 'no_change')  
            if isempty(strfind(type_changes,'ORIENTATION')) == 0  ||...
                    isempty(strfind(type_changes,'CONTRAST')) == 0||...
                    isempty(strfind(type_changes,'SHAPE')) == 0
                
                % if there was a change in at least one stimulus feature, 
                % mark response as a miss. 
                marked_response = 'Miss';   
                
                % if there was no change in any stimulus feature, mark
                % response as a correct rejection.
            elseif isempty(strfind(type_changes,'ORIENTATION')) &&...
                    isempty(strfind(type_changes,'CONTRAST')) &&...
                    isempty(strfind(type_changes, 'SHAPE')); 
                marked_response = 'Correct Rejection';    
            end
    end
end

% if within the selective condition;
if block_condition == 2 
    
     % if the block had an orientation attentional focus:
    if isequal(specific_condition,'ORIENTATION') 
        
        % and if the participant detected a change:
        if isequal(response,'change')     
            
            % if the stimulus change in this trial included an orientation change:
            if isempty(strfind(type_changes,'ORIENTATION')) == 0 
                
                % mark this response as a hit (they responded to a change 
                % in the feature of attentional focus).
                marked_response = 'Hit';    
                
            % if the stimulus did not include an orientation change:    
            elseif isempty(strfind(type_changes,'ORIENTATION'))   
                
                 % mark this response as a false alarm. 
                marked_response = 'False Alarm';                
            end
            
            % if the participant did not detect a change:
        elseif isequal(response, 'no_change')    
            
             % if there was a change in orientation:
             if isempty(strfind(type_changes,'ORIENTATION')) == 0 
                 
                % mark response as a miss. 
                marked_response = 'Miss'; 
                
                % if there was not a change in orientation:
            elseif isempty(strfind(type_changes,'ORIENTATION'))
                
                % mark this response as a correct rejection.
                marked_response = 'Correct Rejection';                
             end
        end
        
     % these two loops follow the same format as the above loop. 
     elseif isequal(specific_condition,'CONTRAST')     
         if isequal(response,'change')
            if isempty(strfind(type_changes,'CONTRAST')) == 0 
                marked_response = 'Hit';  
            elseif isempty(strfind(type_changes,'CONTRAST')) 
                marked_response = 'False Alarm';
            end
        elseif isequal(response, 'no_change')  
             if isempty(strfind(type_changes,'CONTRAST')) == 0
                marked_response = 'Miss';
            elseif isempty(strfind(type_changes,'CONTRAST')) 
                marked_response = 'Correct Rejection';
             end
         end
    elseif isequal(specific_condition,'SHAPE')
        if isequal(response,'change')     
            if isempty(strfind(type_changes,'SHAPE')) == 0
                marked_response = 'Hit';  
            elseif isempty(strfind(type_changes,'SHAPE')); 
                marked_response = 'False Alarm';
            end
        elseif isequal(response, 'no_change') 
             if isempty(strfind(type_changes,'SHAPE')) == 0
                marked_response = 'Miss';
            elseif isempty(strfind(type_changes,'SHAPE'))
                marked_response = 'Correct Rejection';
             end
        end
    end
end

if isequal(response, 'no_response')
    marked_response = 'N/A';
end
end

%% --------------------------------------------------------------------- %%