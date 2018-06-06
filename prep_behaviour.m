

function [bh] = prep_behaviour(para,pathData,loadData)
  
  disp('preparing behavioural data...')
  
  pathSave = pathcat(pathData,'behaviour_data.mat');
  if exist(pathSave,'file') && loadData
    disp(sprintf('Behaviour file already exists. Loading from %s...',pathSave))
    load(pathSave)
    return
  end
  
  lr_min = 30;                     %% number of minimal frames considered to be a "long run"-event
  
  bh(para.nSes) = struct;
  
  for s=1:para.nSes
    disp(sprintf('Processing session %d',s))
    
    tic
    %%% load data
    pathSession = pathcat(pathData,sprintf('Session%02d',s));
    bhfile=dir(pathcat(pathSession,'*.txt'));
    fid = fopen(pathcat(pathSession,bhfile.name),'r');
    fgets(fid);
    whole_data=fscanf(fid, '%f', [8, inf]);
    fclose(fid);
    
    %% checking sync signal
    [crop_bhdata, bh(s).duration] = micro_sync_checker(whole_data);
    
    %% defining run/rest epochs
    speed     = crop_bhdata(2,:);
    frame_num = crop_bhdata(4,:);
    position  = crop_bhdata(6,:);
    
    [bh(s).runrest bh(s).all_speed] = runrest_identifier(speed, frame_num);
    
    bh(s).location = zeros(para.nframe,1);
    for i=1:para.nframe
      bh(s).location(i)=mean(position(find(frame_num==i)))+para.pos_offset;
    end
    
    bh(s).longrunperiod = bwareaopen(bh(s).runrest,lr_min);
    longwalkind = find(bh(s).longrunperiod);   %% find all connected areas > 30
      
    %% create dwell time histogram
    bh(s).dwelltime=zeros(1,para.nbin);
    binnum=floor((bh(s).location(bh(s).longrunperiod))/para.binwidth)+1;
    binnum = min(binnum,para.nbin);
    for i=1:length(binnum)
      bh(s).dwelltime(binnum(i)) = bh(s).dwelltime(binnum(i)) + 1/para.f;
    end
    bh(s).norm_dwelltime = bh(s).dwelltime/sum(bh(s).dwelltime);
  end
%    
  save(pathSave,'bh','-v7.3')
  disp(sprintf('Behavioural data saved under %s',pathSave))
%    
  plt = false;
  if plt
%      
    for s=1:para.nSes
      time_arr = linspace(1/para.f,bh(s).duration,para.nframe);
      
      figure('position',[100 100 900 400])
      %% plot position
%        ax1 = subplot(4,1,1);
      
      ymax = ceil(max(bh(s).all_speed));
      
      yyaxis left
      hold on
      bar(time_arr,bh(s).longrunperiod*ymax,1,'FaceColor',[0.8 0.8 0.8])
      plot(time_arr,bh(s).location*ymax/para.totallength,'b')
      plot([0,time_arr(end)],[60,60]*ymax/para.nbin,'r--','LineWidth',2)
      plot([0,time_arr(end)],[20,20]*ymax/para.nbin,'g--','LineWidth',2)
      hold off
      ylim([-ymax ymax])
      yticks([0,ymax])
      yticklabels([0,1])
      ylabel('Position')
      
      yyaxis right
      %% plot speed
%        ax2 = subplot(4,1,2);
%        plot(time_arr,bh(s).speed,'k')
      plot(time_arr,-bh(s).all_speed,'r')
      ylim([-ymax ymax])
      yticks([-ymax,0])
      yticklabels([-ymax,0])
      ylabel('velocity [cm/s]')
      
      xlabel('time [s]')
      %% plot longrunperiod
%        ax3 = subplot(4,1,3);
%        hold on
%        bar(time_arr,bh(s).longrunperiod,1,'r')
%        bar(time_arr,bh(s).runrest*0.5,1,'b')
%        
%        linkaxes([ax1,ax2,ax3],'x')
      
%        subplot(4,2,7)
%        barh(bh(s).norm_dwelltime)
%        title('normalized dwelltime')
      
      waitforbuttonpress;
    end
  end
    
end



function [runrest speed_rs] = runrest_identifier(speed, frame_num)
%
% This program identifies running and resting periods from mouse's runnind
% speed data. First performs moving average with a window of 100ms and
% resamples data at 15Hz by averaging. Then, defines running periods by
% thresholding at 0.5cm/s and pads gaps between running (< 0.33s, i.e. 5 
% consecutive data points). 
%
% INPUT ARGUMENT
% speed: vector of raw data of mouse's running speed (acquired at 50Hz)
% microsync: vector of data of image frame numbers (usually 1-8989, acquired at 50Hz)
%
% OUTPUT ARGUMENT
% runrest: a vector for run/rest status, 1 is running, 0 is resting.
% speed_rs: filtered and resampled speed data
%
% reviewed 11-15-2014

  nframe=8989;
  datapoints=length(speed);
  
%    sm=5;
%    mafilter_sp=ones(1,sm);% moving average window=200ms (samples obtained @ 50Hz (?))
  smooth_speed = imgaussfilt(speed,10);
%    smooth_speed=filter(mafilter_sp,sm,speed);    %% normalized by sm
  
  speed_rs = zeros(1,nframe);
  for i=1:nframe % resampling by averaging
    speed_rs(i)=mean(smooth_speed(find(frame_num==i)));
  end
  
  %% defining run epochs
  runthres=0.5;%(cm/sec)
  runrest=zeros(1,nframe);
  runrest(find(speed_rs >= runthres))=1;
  
  sm = ones(1,8+1);
  runrest = imerode(imdilate(runrest,sm),sm);   %% filling holes of size up to 5 (=1/3 sec)
end