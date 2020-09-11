filename=uigetfile('*.xlsx');
filename1=uigetfile('*.xlsx');
response_time_w_turnover=xlsread(filename);
response_time=xlsread(filename1);
output=setdiff(response_time(:,1), response_time_w_turnover(:,1));
output=output-min(response_time(:,1));
filename=[filename(1:end-5),'_wo.csv'];
fid = fopen(filename,'w'); 
dlmwrite(filename,output,'-append'); 
fclose(fid);


