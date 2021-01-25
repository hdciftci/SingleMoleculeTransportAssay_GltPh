function varargout = lifetimes_perstate_pertrace_double( filename, options )
% meantime nStates gives NAN instead of 0!!!!!
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
%setting parameters
    params.logX = false;
    params.dx = 0.25;  %log-scale bin width (0.1=25%, 0.2=60%, 0.5=3-fold, 1=10-fold)
    params.removeBlinks = true;
    params.normalize = 'state';
    params.fitSingle=true;
    
%% Prompt user for file names if not given.
[varargout{1:nargout}] = deal([]);
if nargin<1,
    filename = getFiles('*.dwt','Choose a dwell-time file');
   if isempty(filename), return; end  %user hit cancel.
end


% If .traces files are given, quietly look for the corresponding .qub.dwt
for i=1:numel(filename)
lastslash_pos = find(filename{i} == '/', 1, 'last');
fname1=sprintf('lowFRET_%s.csv', filename{i}(lastslash_pos+1:end-8));
fname2=sprintf('highFRET_%s.csv', filename{i}(lastslash_pos+1:end-8));
fid1 = fopen(fname1,'w');
fid2 = fopen(fname2,'w');

%% Load the dwell times and remove dark-state dwells.
[dwells,sampling,offsets,model] = loadDWT(filename{i});% dwells is going to be a cell array 1Xtrace number dimensions. if you look at dwells{trace number} youll  see a 2d matrix ( visiting state vs. frame length it stayed there) 
assert( numel(dwells)>0, 'Empty or invalid dwell-time file' );

nStates = numel(model)/2;% numel gives you total element number so when you dived it to 2 you find the number of states

if nargin>=2 && strcmpi(options,'removeBlinks'), % if the number of input arguments is more than or equal to 2 it means you passed some options, and one of those is remove blinks, this removes blinks
    dwells = removeBlinks(dwells,offsets);
end

dwells_trace=cell(numel(dwells),1); % creates a cell array with number of traces (numel(dwells))x1 dimensions
dwells_state=cell(nStates,1);% creates a cell array with nstates x 1 dimensions each nstate matrix is 1x1
for trace=1:numel(dwells) % strating from the first trace we'll do the same thing for all traces
    times=dwells{trace}(:,2); % number of frames the trace stays at each visit 1 column vector
    states=dwells{trace}(:,1); % corresponding states for those visiting times 1 column vector
    for state=1:nStates
        dwells_state{state}=times(states==state)*sampling/1000;% sampling is in ms timescale by dividing it you turn it into seconds and if you multiply it by the frame number you'll get the visit durations in seconds
    end
    dwells_trace{trace}=dwells_state;
end


% Get dwell time limits for setting axes limits later.
maxTime=0;
totalTime=zeros(numel(dwells_trace),nStates); % creating 2d matrices with trace number of rows each column is for a different state
meanTime=zeros(numel(dwells_trace),nStates);
for trace=1:numel(dwells_trace)
   dwellc=dwells_trace{trace};
   maxTime = max( maxTime, max(vertcat(dwellc{:})) );
    %totalTime(trace,:) = cellfun(@sum, dwellc)'; % total time has rows eeach row has the information of total time that that trace spent in each of the states, thus it has state number of columns
    %meanTime(trace,:)  = cellfun(@mean, dwellc)'; % mean time has rows eeach row has the information of mean time that that trace spent in each of the states, thus it has state number of columns
end

sampling=sampling/1000; % CONEVERSION TO SECONDS FROM THIS POINT ON
%% Calculate dwell time bins (EDGES)
if ~params.logX,
    % Linear X-axis in seconds.
    dwellaxis = 0:sampling:maxTime;
else
   
    % Create a log time axis with a fixed number of bins.
    % histcounts uses bin edges of E(k) <= X(i) < E(k+1).
    dwellaxis = log10(sampling):params.dx:log10(maxTime*3);
    
    % Force the bins edges to be exact intervals of the time resolution.
    % The histogram will better sample the discrete nature of the data.
    maxFrames = ceil(maxTime*3/sampling);
    fullaxis = log10( (1:maxFrames)*sampling )';
    dwellaxis = unique( nearestBin(dwellaxis, fullaxis) );
%     dwellaxis = unique( floor(fullaxis/sampling)*sampling );
    
    % Normalization factor to account for varying-sized bins.
    dlx = dwellaxis(2:end) - dwellaxis(1:end-1);
    dlx = [dlx dlx(end)];
end
 
%% Calculate histograms
histograms = cell(numel(dwells),nStates);


for trace=1:numel(dwells),
      
    for state=2:nStates, 
  %keep int states separate !!!
        if state>=3 && state<nStates
            continue;
        end
        if state==2 % since ifs/ofs and ofs/ofs will have ofs dwell time information in it you should analyze them together
        
            dwellc = dwells_trace{trace}{state}; % here we are fetching the visits of one trace to a particular state in seconds, we converted this remember!!! 
            for i=2:size (dwells{trace},1)
                if dwells{trace}(i,1)>2 && dwells{trace}(i,1)<nStates
                    if dwells{trace}(i-1,1)==nStates
                       t=dwells{trace}(i,2)*sampling
                       dwellc=[dwellc;t];
                    end
                end
            end
        end
        if state==nStates % since ifs/ofs and ofs/ofs will have ofs dwell time information in it you should analyze them together
            dwellc = dwells_trace{trace}{state}; % here we are fetching the visits of one trace to a particular state in seconds, we converted this remember!!! 
            for i=2:size (dwells{trace},1)
                if dwells{trace}(i,1)>2 && dwells{trace}(i,1)<nStates
                    if dwells{trace}(i-1,1)==2 || dwells{trace}(i-1,1)==2
                       t=dwells{trace}(i,2)*sampling
                       dwellc=[dwellc;t];
                    end
                end
            end
        end
        
        meanTime(trace,state)=nanmean(dwellc);
        totalTime(trace,state)=sum(dwellc);
        if isempty(dwellc)==1; continue; end
        if nnz(dwellc)<4; continue; end
        % Make linear-scale survival plot.
        if ~params.logX,
            counts = histc( dwellc, dwellaxis );
            histdata = sum(counts) - cumsum(counts); % cumulative counts
            histdata = histdata/histdata(1); % normalization to 1
        else
            counts = histc( log10(dwellc)', dwellaxis );
      
            if size(counts,2)==0
                continue;
            end
            histdata = counts./dlx ; %normalize by log-space bin size
            histdata = histdata/sum(histdata);  %normalize to 1
            
            switch params.normalize
                case {'none','off'}  %raw dwell counts
                    histdata = histdata*ndwells(state);
                case 'state'  %fraction of counts in each bin for this state
                    histdata = 100*histdata;
                case 'file'  %fraction of counts in each bin across entire file
                    histdata = 100*histdata *ndwells(state)/sum(ndwells);
                case 'time'  %fraction of dwells per total observation time
                    histdata = histdata*ndwells(state)/sum(totalTime(file,:));
            end
        end
        
        histograms{trace,state} = to_col(histdata); % the outward facing states and inward facing states are state 2 and state 5 right, at least in the model it is but we will only get survival plots of outward facing state and inward facing state 
        
        
    end

end

%%

if params.logX,
    dwellaxis = 10.^dwellaxis;
end

output = {dwellaxis,histograms,meanTime,sampling};
meanTime(isnan(meanTime))=0;
[varargout{1:nargout}] = output{1:nargout};
[lifetimes]=lifetime_pertrace(dwellaxis,histograms,meanTime,sampling,maxTime);
nonzero_high=sum(lifetimes.high~=0,2)>=0;
nonzero_low=sum(lifetimes.low~=0,2)>=0;
cHeader = {'fit_type' 'trace ID' 'meanTime' 'Rsquared' 'Tau_one' 'Tau_two' 'amp_of_Tau_one' 'Rsquared_of_rejected_fit'};
commaHeader = [cHeader;repmat({','},1,numel(cHeader))]; %insert commaas
commaHeader = commaHeader(:)';
textHeader = cell2mat(commaHeader); %cHeader in text with commas
fprintf(fid1,'%s\n',textHeader);
fprintf(fid2,'%s\n',textHeader);
dlmwrite(fname1,lifetimes.low(nonzero_low,:),'-append');
dlmwrite(fname2,lifetimes.high(nonzero_high,:),'-append');
fclose(fid1);
fclose(fid2);


end
end
%% Fit survival plots to exponential functions
function varargout= lifetime_pertrace(dwellaxis,histograms,meanTime,sampling,maxTime)
params.removeBlinks = true;  % merge blinks into previous dwell

params.logX = false;  %no Sine/Sigworth transform.
params.fitSingle = false;  %otherwise, double exponential fitting.

nStates = size(histograms,2); % here we are going to be extracting the dwell times of inward facing and outward facing states we designed the loaddwelltimes function to concatanate the dwell times of asymetric states with either inward facing or outward facing state
nTraces=size(histograms,1);  % first dimension of survival cell array gives the number of traces

fit_type=1;
lifetimes.high= zeros(nTraces,8);% average lifetimes (fit)
lifetimes.low= zeros(nTraces,8);
fits = cell(nTraces,nStates);    % fit structures for plotting

for trace=1:nTraces,
    for state=2:nStates,
       
      
        if state>=3 && state<nStates
            continue;
        end
         if state==2
            fieldname='low';
         else
            fieldname='high';
         end
        y = histograms{trace,state};
       
        if isempty(y)==1 % this is when there is histogram created for that state and trace; could be due to two reasons either there are no visits to that state or the number of visits is less than 4 so a meaningful histogram cannot be created.
            lifetimes.(fieldname)(trace,1)=1;
            lifetimes.(fieldname)(trace,2)=trace;
            lifetimes.(fieldname)(trace,3)=meanTime(trace,state);
            disp('fit function state is empty');
            continue;
        end
        if nnz(y)==1
            disp('one item')
            lifetimes.(fieldname)(trace,1)=1;
            lifetimes.(fieldname)(trace,2)=trace;
            lifetimes.(fieldname)(trace,3)=meanTime(trace,state);
            continue;
        end 
        plen = find(y>0.01,1,'last');
        if dwellaxis(plen)<maxTime/3 
        x = dwellaxis(1:plen*3)';
        y = y(1:plen*3);
        else
        x = dwellaxis(1:plen)';
        y = y(1:plen);
        end
        if nnz(y)==sum(y(1:nnz(y)),1) % if there is only one dwell time 
            disp('one item after cut')
            lifetimes.(fieldname)(trace,1)=1;
            lifetimes.(fieldname)(trace,2)=trace;
            lifetimes.(fieldname)(trace,3)=meanTime(trace,state);
            continue;
        end
        if size(y,1)<=4
            disp('less than three data points')
            lifetimes.(fieldname)(trace,1)=1;
            lifetimes.(fieldname)(trace,2)=trace;
            lifetimes.(fieldname)(trace,3)=meanTime(trace,state);
            continue;
        end
        % Fit the survival plot to a double exponetial first
           disp('fitting double')
           ft2 = fittype('(a-p)*exp(b*x)*S+ (a-p)*exp(c*x)*(1-S)+ p','independent','x','dependent','y');
           opts2 = fitoptions( 'Method', 'NonlinearLeastSquares' );
           opts2.Lower = [0 1 -Inf -Inf 0];
           opts2.StartPoint = [0.5 1 meanTime(trace,state) meanTime(trace,state) 0];
           opts2.Upper = [1 1 0 0 0];
           [fits{trace,state},gof2] = fit( x, y, ft2,opts2);
           %plot(fits{trace,state},x,y);
           coeffs2 = coeffvalues( fits{trace,state} );
           coeffnames2=coeffnames(fits{trace,state}); % S a b c p
           % if the tau's are less than 3 fold different fit the lifetimes to single exponential
       % if coeffs(3)/coeffs(4)<3 && coeffs(3)/coeffs(4)>=1 || coeffs(3)/coeffs(4)>0.3 && coeffs(3)/coeffs(4)<=1
           disp('fitting single')
           ft = fittype('(a-p)*exp(b*x)+p','independent','x','dependent','y');
           opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
           opts.Lower = [1 -Inf 0];
           opts.StartPoint = [1 meanTime(trace,state) 0];
           opts.Upper = [1 0 0];
           [fits{trace,state},gof1] = fit( x, y, ft,opts);
           coeffs1 = coeffvalues( fits{trace,state} );
           coeffnames1=coeffnames(fits{trace,state}); % a b p
        %end
        
        if abs(gof1.adjrsquare-gof2.adjrsquare)/max(gof1.adjrsquare,gof2.adjrsquare)>0.1
            rsquare=max(gof1.adjrsquare,gof2.adjrsquare);
            if rsquare==gof1.adjrsquare
                fit_type=1;
            else 
                fit_type=2;
            end
        else 
            rsquare=gof1.adjrsquare;
            fit_type=1;
        end

        lifetimes.(fieldname)(trace,1)=fit_type;
        lifetimes.(fieldname)(trace,2)=trace;
        lifetimes.(fieldname)(trace,4)=rsquare;
       
        if isnan(lifetimes.(fieldname)(trace,4))==1
           lifetimes.(fieldname)(trace,3)=meanTime(trace,state) ;
        end
        
        if isnan(lifetimes.(fieldname)(trace,4))==0 
            if lifetimes.(fieldname)(trace,4) >=0.15 && fit_type==2
               lifetimes.(fieldname)(trace,5) = -1/coeffs2(1,3);
               lifetimes.(fieldname)(trace,6) = -1/coeffs2(1,4);
               lifetimes.(fieldname)(trace,7) =coeffs2(1,1);
              if isnan(gof1.adjrsquare)==0
               lifetimes.(fieldname)(trace,8) =gof1.adjrsquare;
              end
            end
            if lifetimes.(fieldname)(trace,4) >=0.15 && fit_type==1
               lifetimes.(fieldname)(trace,3) = -1/coeffs1(1,2);
            end
            if lifetimes.(fieldname)(trace,4)<0.15
                lifetimes.(fieldname)(trace,3)=meanTime(trace,state);
            end
            
        end
        
    end

end  % for each sample
%% Write the fitted lifetimes and their sigworth distribution to csv files
output = {lifetimes};
[varargout{1:nargout}] = output{1:nargout};

end

 
function [newVal,idx] = nearestBin( values, bins )
% For each VALUE, find the BIN with the closest value.

newVal = zeros( size(values) );
idx = zeros( size(values) );

for i=1:numel(values),
    [~,idx(i)] = min( abs(bins-values(i)) );
    newVal(i) = bins(idx(i));
end

end



