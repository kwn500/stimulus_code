function [response_value, response, myBeep, answer] = response_classify(keyLabel, ii, feature_task)
    
if keyLabel == 'u';   %if they pressed the 'l' key; 
    if isequal(feature_task, 'Orientation'); 
        response = 'clockwise'; response_value = (ii<=.5); 
    elseif isequal(feature_task, 'Contrast');
        response = 'higher contrast'; response_value = (ii>.5); 
    elseif isequal(feature_task, 'RAM');
        response = 'spikier'; response_value = (ii>.5); 
    end
elseif keyLabel == 'n';    %if they pressed the 'a' key
    if isequal(feature_task, 'Orientation'); 
        response = 'anticlockwise'; response_value = (ii>.5); 
    elseif isequal(feature_task, 'Contrast');
        response = 'lower contrast'; response_value = (ii<=.5); 
    elseif isequal(feature_task, 'RAM');
        response = 'smoother'; response_value = (ii<=.5);
    end
end
        
if response_value == 1  %if participant was correct, play high pitched beep. 
   myBeep=MakeBeep(500, 0.2, 10000); answer = 'correct';
elseif response_value == 0
    myBeep = MakeBeep(500, 0.1, 100000); answer = 'incorrect'; 
end
end
%% --------------------------------------------------------------------- %%

