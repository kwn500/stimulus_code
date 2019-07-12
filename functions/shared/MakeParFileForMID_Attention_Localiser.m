
function MakeParFileForMID_Attention_Localiser

% Generate a .par file describing the timing in the Localiser for the attention to Motion In Depth experiment.
% This is in tab-delimited text format, as expected by mrVista when setting up GLM parameters.
% This includes the labelling of all fixation/baseline blocks as 'Fix'
% Only 1 file should be needed, since the localiser (should) be the same for all subjects.
% NOTE: there are no 'dummy' TRs in the localiser, the first 4 TRs being fixation/baseline.
% At first it was thought this wouldn't matter, because the first few noisy TRs would be lumped in with the fixation/baseline condition.
% However, the inhomogeneities of these first 4 dummies could stuff up the analysis, so we will treat the first 4 TRs as dummies, as if they were added 
% manually to the start of a timing sequence as they would be if it were generated using Optseq2.
% This means we will need to adjust the timing so that TIME 0 is actually the start of the 5th TR (as we will discard the first 4 TRs).
% R Maloney 7 Feb 2017

%The file name to be saved:
parFileName = 'MID_attn_localiser_design_no_dummy_vista.par';

%Most of below comes directly from Run_EoO_localiser where the timing is determined
%Set a bunch of variables important in determining the timing of stimuli:
DummyTRs = 4; % dummy volumes at the beginning;
BlockLengthTRs = 4; %number of TRs in an block
TR = 3; %length of volume (TR), in sec
BlockLengthSec = BlockLengthTRs * TR;

%Set up the block conditions:
%These provide the given condition for each event of the scan, and are 'read in' when displaying the stimuli
Inner = 1;
Outer = 2;
Fix = 0; %NOTE: In MrVista conventions, 'fixation' periods are always labelled as 'Fix' and given code 0

%Work out a single cycle:
BlockOrder = [Fix; Inner; Fix; Outer];

%Replicate for 7 cycles, and place one final Fixation block at the end:
BlockOrder = [repmat(BlockOrder,7,1); Fix];

%Insert the durations for each Block:
BlockOrder(:,2) = ones(length(BlockOrder),1) * BlockLengthSec;

% NOW: remove the 4 dummy TRs from the timing sequence (ie one fixation block =  4 TRs = 12 sec):
BlockOrder = BlockOrder(2:end,:);

% Determine the ONSET TIMES for each BLOCK. 
% Ie the first is set as time 0: this is the same as with the event-related scans as determined using Optseq, where the dummy scans are also inserted manually (ie not by Optseq)
% We want the onset times for each block only, so subtract the length of the dummy scans, and one block, so the final value is the ONSET time of the FINAL block.
% It is as if the dummies never happened. 
BlockOnsets = 0:BlockLengthSec: (348-DummyTRs*TR) -BlockLengthSec; %(348 is total length, IN SEC of scan).

%Now make a cell array containing the condition names as strings.
%This seems optional, but we will include it because the mrVista files seem to have it:
for ii = 1: length(BlockOrder)
    
    switch BlockOrder(ii,1)
        case Fix
            BlockName{ii} = 'Fix'
        case Inner
            BlockName{ii} = 'Inner'
        case Outer
            BlockName{ii} = 'Outer'
    end
end


%Now open and save the appropriate text file:
fileID = fopen(parFileName,'wt'); %open the file to be saved, 'w' for writing, 't' for text

%Set up the formatting for the file, this is the order of each item in each row of the file.
% 3.2f means fixed point notation, so up to 3 integer parts and 2 decimal places (for onset time)
% d means a signed integer  (for block code)
% s means a string          (for block/condition name)
% \t means separate by a horizontal tab 
% \n means go to a new line.
formatSpec = '%3.2f\t %d\t %s\t\n'; 

%loop across each line and add it to the file:
for ii = 1:length(BlockOrder)
    
    fprintf(fileID, formatSpec, BlockOnsets(ii), BlockOrder(ii), BlockName{ii});
    
end

fclose(fileID); %close the file when done. It should be saved in the pwd


    






