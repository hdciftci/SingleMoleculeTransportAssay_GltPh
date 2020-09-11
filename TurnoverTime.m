function varargout = TurnoverTime(varargin)

% Default parameter values
params.truncateLen = 1000;  %frames to calculate over
params.min_fret = 0.5;  % minimum fret value, below which we assume there is no FRET.
createplot=true;% should a plot be created of each trajectory with its fit?
params.dx=0.25;% bin width for when you want to create a distribution
liposome_size_check=true;% should liposomes smaller than a certain size be excluded right now the cutoff is at 25 nm radius
limx=100;% limit of x axis of the plots, if you wish to create them
%% Process input arguments
narginchk(0,3);
nargoutchk(0,1);
[varargout{1:nargout}] = deal([]);
[cax,args] = axescheck(varargin{:});

switch numel(args)
    case 0
        files = getFiles();
    case 1
        files = args{1};
    case 2
        [files,inputParams] = args{:};
        params = mergestruct(params, inputParams);
end

if ~iscell(files), files={files}; end
nFiles = numel(files);
if nFiles==0,  return;  end

for i=1:numel(files)% this file should only have traces that respond 
    % Load FRET data and truncate to target length
    data = loadTraces( files{1} );
    data.nFrames = params.truncateLen;
  % N is the number of frames you want to use to fit the data, it helps to
  % get a better fits especially when the Thc (half cycle time) is short
  % and the trajectory itself is very long.
 
  coefficients_2=zeros(size(data.fret,1),21);% whether a trajectory passes the criteria to be fitted or not all the information of trajectories in datasets are collected here
  coefficients_2(1,19)=params.truncateLen;
  coefficients_2(1,20)=params.min_fret;
  r2_ktot=zeros(size(data.fret,1),20);% fit information of trajectories that have fits with an r square value above a certaain cutoff
  r2_ktot_cutoff=0.2;% r square value cutoff
  r2_ktot_not=zeros(size(data.fret,1),20);% information of trajectories that either did not pass the criteria to be fitted or could not be fitted with an r-square above cutoff
  r2_ktot_i=1;% counter for r2_ktot
  r2_ktot_j=1;%counter for r2_ktot_not
  filename=uigetfile('*.xlsx');% user must input the response time excel sheet, with row numbers corresponding to trajectories numbers in the traces file selected
  response_time=xlsread(filename);
  kd=5.6;% affinity of the sensor used in the experiment
  
  w = waitbar(0,'Please wait...');
    for trace=1:size(data.fret,1)
        %CRITERIA: if there is no corresponding response time for the trajectory continue to the next one. 
        if isnan(response_time(trace,1))==1 
           coefficients_2(trace,13)=1;
           continue;
        end 
        %fetching fret values bigger than params.min_fret (you can change this value in the beginning of the function) in each trace
        nonzero= data.fret(trace,:) >= params.min_fret; 
        %forming time(sec) FRET matrix named fret_trace 
        fret_trace=horzcat(transpose(data.time(:,nonzero)/1000),transpose(data.fret(trace,nonzero)));
        fret_trace=fret_trace(1:end-1,:);
        %calculating the temporal resolution of the trace
        resolution=fret_trace(2,1)-fret_trace(1,1);
        %x is the time component of the trace
        x=fret_trace(:,1);
        %response times are given as seconds in the excel sheet. convert them to frame number
        
        response_timee=max(response_time(trace,1));
        if response_timee==1000
            continue;
        end
        response_frame=find(abs(x-response_timee)<0.05);
        %y is the fret component of the trace
        y=fret_trace(:,2);
        %trim the traces to align time 0 with response time
           x_trimmed_jump=fret_trace(response_frame+1:length(x),1)-(response_time(trace,1)+resolution);
           y_trimmed_jump=fret_trace(response_frame+1:length(y),2);
        
        %CRITERIA: if the trajectory does not have any frames after the
        %jump continue to the next trajectory
        if isempty( y_trimmed_jump) || size(y_trimmed_jump,1)<7
            coefficients_2(trace,14)=1;
            continue;
        end
        %calculate initial FRET value; the first 5 frames right before the
        %response frame are averaged to find this value, if there is less
        %than 5 frames before response time all those are averaged, however
        %many they may be.
        if response_frame>10
           init_f=mean(fret_trace(response_frame-10:response_frame,2));
        else
           init_f=mean(fret_trace(1:response_frame,2));
        end
        %maximum FRET efficiency change between saturated (closed) and
        %unbound(open) sensor
        diff=0.22;
        %CRITERIA: if trajectory has an initial FRET higher than 0.65 or
        %lower than 0.55 that is not normal continue to the next one, this might happen because of
        %noise or misfolding or there is already some L-asp in the liposome
       
        if init_f>0.67 || init_f<0.55
            coefficients_2(trace,15)=1;
            continue;
        end
        %the following condition is put due to trajectories that start out
        %very noisy, when this happens it messes up the intial fret but im
        %not sure if this should be in there 
        %if init_f<0.6
        %   init_f=0.6;
        %end
        % saturated sensor FRET efficiency (max_F) is calculated. 
        max_F=init_f+diff;
        %first increase in FRET(jump_F) occurs when the first asp comes into the
        %liposome. depending on the size of the liposome first step FRET
        %changes.it is assumed that the first step FRET value is the avg of
        %next three frames after the response frame.
        jump_F=mean(y_trimmed_jump(1:3,1));
        % CRTERIA: if the trajectory is very noisy calculated initial FRET is sometimes higher than jump FRET value,for such cases jump percentage becomes negative/meaningless 
        
        if init_f>jump_F 
            coefficients_2(trace,16)=1;
            continue;
        end 
        jump_size=jump_F-init_f;
        jump_percentage=jump_size/diff;
        C=10^30/(6.02*10^23*4/3*pi);
        jp= jump_percentage;
        r=(C/( (jp*kd)/(jp^2 - 2*jp + 1)))^(1/3);
        %CRITERIA: of the jump_percentage is more than 80% or if the jump_F
        %is more than 0.76, it is either because there were already some
        %L-asp molecules inside liposome or the liposome is too small, in
        %either case these trajectories cannot be used to extract further
        %information.
        if liposome_size_check
        if jump_F> 0.76 || jump_percentage>0.8
            coefficients_2(trace,17)=1;
            continue;
        end
        end
        %calculation of 1 molecule of sensor inside liposomes. this depends on r (radius of liposome) 
        etotal=((1/(6.02*10^23))/(4/3*pi*r^3*10^(-24)))*10^6;
       
        idx=max(find(y_trimmed_jump>=max(0.85,max_F),15));
        N=min(size(y_trimmed_jump,1),idx+20);
        if N<5
            continue;
        end
        if isempty(N)
            N=size(y_trimmed_jump,1);
        end 
        coefficients_2(trace,18)=N-size(y_trimmed_jump,1);
        %N=size(y_trimmed_jump,1);
        %fit of individual traces, the fits only optimizes ktot(1/turnover
        %time) and max_F, other parameters are constrained.  
        ft = fittype( '[((etotal+[x*ktot+etotal]+kd)-((etotal+[x*ktot+etotal]+kd)^2-(4*[x*ktot+etotal]*etotal))^(1/2))/ 2]/etotal*(MAX_F-init_f)+init_f', 'independent', 'x', 'dependent', 'y' );
        opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
        opts.Display = 'Off';
        opts.Lower = [max_F-0.03 etotal init_f kd  0 ];
        opts.StartPoint = [max_F etotal init_f kd  1 ];
        opts.Upper = [max_F+0.07  etotal init_f kd Inf ];
        [fitresult_rate,gof] = fit( x_trimmed_jump(1:N), y_trimmed_jump(1:N), ft, opts );
       coeff=coeffvalues(fitresult_rate);
       coefficients_2(trace,3:7)=coeff(1:5);
       coefficients_2(trace,1)=gof.rsquare;
       coefficients_2(trace,2)=r;
        %conversion of turnover time unit to # of molecules/second: based on how much of a conc. one molecule creates and ktot (transport rate constant in umolar/sec units) value.
        coefficients_2(trace,8)=etotal/coefficients_2(trace,7);
        coefficients_2(trace,9)=jump_F;
        coefficients_2(trace,10)=response_time(trace,1)-min(response_time(:,1));
        coefficients_2(trace,21)=response_time(trace,1);
        coefficients_2(trace,11)=diff;
        coefficients_2(trace,12)=jump_size;
      if coefficients_2(trace,1)<r2_ktot_cutoff
        r2_ktot_not(r2_ktot_j,1)=trace;
        r2_ktot_not(r2_ktot_j,2:11)=coefficients_2(trace,1:10);
        r2_ktot_not(r2_ktot_j,12)=jump_size;
        r2_ktot_j=r2_ktot_j+1;
      end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
      if coefficients_2(trace,1)>r2_ktot_cutoff
        r2_ktot(r2_ktot_i,1)=trace;
        r2_ktot(r2_ktot_i,2:11)=coefficients_2(trace,1:10);
        r2_ktot(r2_ktot_i,12)=jump_size;
        r2_ktot(r2_ktot_i,13)=coefficients_2(trace,21);
        r2_ktot_i=r2_ktot_i+1;
      
       if createplot 
        %create plot
        figure1=figure( 'Name', ' fits' );
        % Create axes
        axes1 = axes('Parent',figure1);
        hold(axes1,'on');

        plot(x(1:length(x),1),y(1:length(y),1),'DisplayName','avgFRET vs. time(sec)','LineWidth',0.5,...
        'Color',[0 0 1]);
        plot((x_trimmed_jump(1:length(x_trimmed_jump),1)+response_time(trace,1)+resolution),fitresult_rate(x_trimmed_jump(1:length(x_trimmed_jump),1)),'DisplayName','ktot(full turnover) fit','LineWidth',1.5,...
        'Color',[1 0 0]);
        title({'trace',trace},'FontSize',24);
       % Create xlabel
        xlabel('Time(sec)');
        % Create ylabel
        ylabel('FRET','Interpreter','none');
        xlim(axes1,[0 limx]);
        ylim(axes1,[0.4 1]);
        box(axes1,'on');
        % Set the remaining axes properties
        set(axes1,'FontSize',24);
        % put waiting time and radius of liposome to plot
       ylim_text=get(gca,'ylim');
       xlim_text=get(gca,'xlim');
       txt=sprintf('waiting time = %.1f sec \nliposome size = %.1f nm\ntransport time = %.1f sec\nraquare=%.1f', response_time(trace,1),r,  coefficients_2(trace,8), coefficients_2(trace,1));
       text(xlim_text(2)-490,ylim_text(2)-0.08,txt,'FontSize',14);
       fname='fits';
       filename=sprintf('trace %d',trace);
       saveas(figure1,fullfile(fname, filename),'jpeg')
       
       end
      end
     
     waitbar(trace/size(data.fret,1),w)  
    end
    delete(w)
     
end

cHeader = {'rsquare' 'radius of liposome' 'MAX_F' 'etotal' 'init_f' 'kd' 'ktot' 'sec/asp' 'jump_f' 'response time' 'diff' 'jump_size' 'No response time detected' 'Trace is too short' 'Contamination of misfolding' 'Noisy trace' 'Liposome is smaller than 25 nm' '#frames fitted-#total frame' '# frames provided as input' 'photobleach FRET cutoff' };
commaHeader = [cHeader;repmat({','},1,numel(cHeader))]; %insert commaas
commaHeader = commaHeader(:)';
textHeader = cell2mat(commaHeader); %cHeader in text with commas

fid = fopen('coefficients.csv','w'); 
fprintf(fid,'%s\n',textHeader);
fclose(fid);

dlmwrite('coefficients.csv',coefficients_2,'-append') 

cHeader = {'trace #' 'rsquare' 'radius of liposome' 'MAX_F' 'etotal' 'init_f' 'kd' 'ktot' 'sec/asp' 'jump_f' 'response time' };
commaHeader = [cHeader;repmat({','},1,numel(cHeader))]; %insert commaas
commaHeader = commaHeader(:)';
textHeader = cell2mat(commaHeader); %cHeader in text with commas


%{
maxTime=0;
dwellc=r2_ktot(:,9);
maxTime = max( maxTime, max(dwellc));

dwellaxis = log(resolution):params.dx:log(maxTime*3);
    
    % Force the bins edges to be exact intervals of the time resolution.
    % The histogram will better sample the discrete nature of the data.
    maxFrames = ceil(maxTime*3/resolution);
    fullaxis = log( (1:maxFrames)*resolution )';
    dwellaxis = unique( nearestBin(dwellaxis, fullaxis) );
    transpose(dwellaxis)

    
    % Normalization factor to account for varying-sized bins.
    dlx = dwellaxis(2:end) - dwellaxis(1:end-1);
    dlx = [dlx dlx(end)];

    counts = histc( log(dwellc)', dwellaxis );
    histdata = counts./dlx;  %normalize by log-space bin size
    histdata = histdata/sum(histdata);  %normalize to 1
    histdata=histdata.*100;
    antilog_dwellaxis=exp(dwellaxis);
    r2_ktot(1:length(histdata),15)=transpose(histdata) ;  
    r2_ktot(1:length(dwellaxis),13)=transpose(dwellaxis);
    r2_ktot(1:length(dwellaxis),14)=transpose(antilog_dwellaxis);
%}
fid = fopen('rate_constants.csv','w'); 
fprintf(fid,'%s\n',textHeader);
dlmwrite('rate_constants.csv',r2_ktot,'-append') 
fclose(fid);
fid = fopen('rate_constants_not.csv','w'); 
fprintf(fid,'%s\n',textHeader);
dlmwrite('rate_constants_not.csv',r2_ktot_not,'-append') 
fclose(fid);

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







