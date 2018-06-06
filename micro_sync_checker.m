function [data,duration] = sync_checker(whole_data)
%
% This program extracts the part of image acquisition from the whole behavioral data. 
% The existance of frame numbers is examined and if not present, new
% numbers are calculated from the TTL signal and assigned.
%
% INPUT ARGUMENT
% whole_data: whole data matrix read from the behavior data file
%
% OUTPUT ARGUMENT
% starttime: absolute time point when imaging started  
% endtime: absolute time point when imaging ended
% duration: duration of image acquisition
% startpoint: absolute time point when imaging started
% endpoint: absolute time point when imaging ended
% data: the extracted part of data matrix which corresponds to image acquisition
%
% reviewed 11-15-2014

data1=whole_data(1,:);% time
data4=whole_data(4,:);% frame number
data8=whole_data(8,:);% microscope TTL

startframe=1;
endframe=8989;
sync_sum=sum(data4);
sync_last=data4(end);

%%
if sync_sum~=0 & endframe==sync_last
    
    startpoint=find(data4==startframe,1,'first');
    endpoint=find(data4==endframe-1,1,'last')+3;
    datapoints=endpoint-startpoint+1;
    datapoints_per_frame=datapoints/endframe;
    
    starttime=data1(startpoint);
    endtime=data1(endpoint);
    duration=endtime-starttime;

elseif  sync_sum==0 | endframe~=sync_last
    
    startpoint=find(data8<1,1,'first');% find the first point or low TTL    
    diffsig=diff(data8);
    endpoint=find(diffsig > 1.5,1,'last');% find the last point or high TTL

    datapoints=endpoint-startpoint+1;
    datapoints_per_frame=datapoints/endframe; 
    
    starttime=data1(startpoint);
    endtime=data1(endpoint);
    duration=endtime-starttime;
          
    if sync_sum==0
            disp('WARNING --- Sync signal is missing')
    elseif endframe~=sync_last
            disp('WARNING --- Frame number mismatch')
    end
    
end

%% time normalization
data=whole_data(:,startpoint:endpoint);%% crop the data
data(1,:)=data(1,:)-data(1,1);%% time normalization

%% assign new frame numbers when sync signals are missing or end of signals is not the last frame number 
if  sync_sum==0 | endframe~=sync_last
    data(4,:)=ceil((1:datapoints)/datapoints_per_frame);  
    disp('New frame numbers assigned')
    disp('   ')
end


end

