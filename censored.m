function [dwellc,sampling,model] = loadDwelltimes( filename, options )
% loadDwelltimes  Combined list of dwell-times for each state.
%
%   [DWELLS,SAMPLING,MODEL] = loadDwelltimes(FILENAME) creates the cell array
%   DWELLS with a list of all dwell times (in seconds) from each state. 
%   SAMPLING is the time resolution in ms and MODEL is the FRET model.
%
%   [...] = loadDwelltimes(...,'removeBlinks') will remove dwells in the 
%   zero-FRET state (state 1) associated with blinking.
%
%   See also: removeBlinks, dwellhist, lifetime_exp.

%   Copyright 2007-2015 Cornell University All Rights Reserved.



%% Prompt user for file names if not given.
if nargin<1,
    filename = getFile('*.dwt','Choose a dwell-time file');
    if isempty(filename), return; end  %user hit cancel.
end


% If .traces files are given, quietly look for the corresponding .qub.dwt
[p,f,e] = fileparts(filename);
    
if ~strcmp(e,'.dwt'),
    filename = fullfile(p,[f '.qub.dwt']);
    if ~exist(filename, 'file'),
        filename = fullfile(p,[f '.dwt']);
    end
    if ~exist(filename,'file'),
        error('Input file is not a .dwt and no associated .dwt file could be found');
    end
end



%% Load the dwell times and remove dark-state dwells.
[dwells,sampling,offsets,model] = loadDWT(filename);
assert( numel(dwells)>0, 'Empty or invalid dwell-time file' );

nStates = numel(model)/2;

if nargin>=2 && strcmpi(options,'removeBlinks'),
    dwells = removeBlinks(dwells,offsets);
end
dwellc=cell(size(dwells,2),nStates);
% Merge dwells from all traces
for k=1:size(dwells,2)
size(dwells,2);% trace number
a=dwells{1,k};


states = a(:,1);
times = a(:,2);

% Build list of dwell times in each class
ac = cell(nStates, 1);


for i=1:nStates,
    ac{i} = sum(times(states==i))*sampling/1000;
end

dwellc(k,:)=ac';

end

response_time_file=uigetfile('*.xlsx');
response_time=xlsread(response_time_file);
response_time
min(response_time(:,1))
censored=cell2mat(dwellc(:,2))-min(response_time(:,1));
fid = fopen('censored_lifetime.csv','w'); 
dlmwrite('censored_lifetime.csv',censored,'-append') 
fclose(fid);
