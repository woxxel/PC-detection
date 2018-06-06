

function [PC_fields] = placecell(pathData,clusters)
  
  if nargin < 2 || isempty(clusters)
    load(pathcat(pathData,'clusters.mat'))
  end
  
  nSes = size(clusters,2);
  para = struct;
  para.f = 15;                 %% Hz frame rate
  nframe = 8989;
  totallength = 1600;     %% length of the linear track
  pos_offset = totallength;   %% offset of position values
  para.nbin=80;                %% number of divisions of the linear track
  para.binwidth = totallength/para.nbin;
  
  lr_min = 30;                     %% number of minimal frames considered to be a "long run"-event
  prc = 20;               %% percentile for computing noise level of S
  nsd = 3;                %% number of std above noise level to be considered "activity"
  
  repnum=500;            %% number of randomized traces to compare MI to
  
  PC_fields(length(clusters)) = struct;
  
  for c = 1:length(clusters)
    PC_fields(c).firingmap = zeros(nSes,para.nbin);
    PC_fields(c).firingmap(:) = NaN;
    
    PC_fields(c).dwellmap = zeros(nSes,para.nbin);
    PC_fields(c).dwellmap(:) = NaN;
    
    PC_fields(c).status = false(nSes,1);
    PC_fields(c).mean_fr = zeros(nSes,1);
    PC_fields(c).mean_fr(:) = NaN;
    
    PC_fields(c).max_fr = zeros(nSes,1);
    PC_fields(c).max_fr(:) = NaN;
    
    PC_fields(c).max_pos = zeros(nSes,1);
    PC_fields(c).max_pos(:) = NaN;
    
    PC_fields(c).MI = zeros(nSes,1);
    PC_fields(c).MI(:) = NaN;
    
    PC_fields(c).MI_frac = zeros(nSes,1);
    PC_fields(c).MI_frac(:) = NaN;
  end
  
  
  for s=1:nSes
    disp(sprintf('Processing session %d',s))
    
    tic
    %%% load data
    pathSession = pathcat(pathData,sprintf('Session%02d',s));
    bhfile=dir(pathcat(pathSession,'*m.txt'));
    
    fid = fopen(pathcat(pathSession,bhfile.name),'r');
    fgets(fid);
    whole_data=fscanf(fid, '%f', [8, inf]);
    fclose(fid);
    
    %% checking sync signal
    [starttime, endtime, duration, startpoint, endpoint, crop_bhdata] = micro_sync_checker(whole_data);
    
    %%% processing behavioural data
    bh = struct;
    
    %% defining run/rest epochs
    speed     = crop_bhdata(2,:);
    frame_num = crop_bhdata(4,:);
    position  = crop_bhdata(6,:);
    
    [runrest all_speed] = runrest_identifier(speed, frame_num);
    
    bh.location = zeros(nframe,1);
    for i=1:nframe
      bh.location(i)=mean(position(find(frame_num==i)))+pos_offset;
    end
    
    bh.longrunperiod = bwareaopen(runrest,lr_min);
    longwalkind = find(bh.longrunperiod);   %% find all connected areas > 30
      
    %% create dwell time histogram
    bh.dwelltime=zeros(1,para.nbin);
    for i=1:length(longwalkind)
        binnum=floor((bh.location(longwalkind(i)))/para.binwidth)+1;
        if binnum <= para.nbin & ~isnan(binnum)
          bh.dwelltime(binnum) = bh.dwelltime(binnum) + 1;
        end
    end
    bh.norm_dwelltime = bh.dwelltime/sum(bh.dwelltime);
    
    %% set parameters
    IDs = [clusters(:,s).ROI_ID];
    detected = ~isnan(IDs);
    idx_clusters = find(detected);
    IDs = IDs(detected);
    ncell = sum(detected);
    
    endframe=nframe;
    
    %%% get firing rates
    for c = idx_clusters
      modeS = prctile(clusters(c,s).S(clusters(c,s).S>0),prc);
      clusters(c,s).activity = floor(clusters(c,s).S/(modeS*nsd));
      clusters(c,s).spikes = sum(clusters(c,s).activity);
      PC_fields(c).mean_fr(s) = clusters(c,s).spikes/duration;
    end
    
    %% place cell detection from mutual information
    MI=zeros(1,ncell);
    status_tmp = zeros(ncell,1);
    for c=idx_clusters
        
        PC_fields(c).dwellmap(s,:) = bh.dwelltime;
        [MI, PC_fields(c).firingmap(s,:)] = calc_MI(bh,clusters(c,s).activity,para);
        
        %% calculate MI from random data 
        MI_rand=zeros(1,repnum);
        rand_ind_vec=randperm(nframe);
        
        for L=1:repnum
            rand_ind=rand_ind_vec(L);
            activity_rand = [clusters(c,s).activity(rand_ind:nframe), clusters(c,s).activity(1:rand_ind-1)];% random permutation --- no ?!?!
            
            [MI_rand(L), ~] = calc_MI(bh,activity_rand,para);
        end
        
        status_tmp(c) = MI > prctile(MI_rand,95);    %% comparison with 95 percentile
        PC_fields(c).status(s) = status_tmp(c);
        
        PC_fields(c).MI(s) = MI;
        PC_fields(c).MI_frac(s) = MI / prctile(MI_rand,95);
        
%          if PC_fields(c).status(s)
        [PC_fields(c).max_fr(s),PC_fields(c).max_pos(s)] = max(PC_fields(c).firingmap(s,:));
%          else
%            PC_fields(c).max_pos(s) = NaN;
%          end
    end % c
    
    toc
    idx_PC = find(status_tmp);
    pc_ct = length(idx_PC);
    disp(sprintf('found %d place cells from %d neurons',pc_ct,ncell))
  end
  
  savePath = pathcat(pathData,'place_fields.mat');
  save(savePath,'PC_fields','-v7.3')
  disp(sprintf('data saved under %s',savePath))

end


function [MI, norm_firingmap] = calc_MI(bh,activity,para)
  
  MI = 0;
  firingmap=zeros(1,para.nbin);
  
  longrunspikes = activity & bh.longrunperiod;
  spikes = activity(longrunspikes);
  binnum = floor(bh.location(longrunspikes)/para.binwidth)+1;
  binnum = min(binnum,80);
  %% check not really needed
  if any(binnum>para.nbin | isnan(binnum))
    disp('somethings wrong')
    binnum
    bh.location(longrunspikes)
    binnum(binnum>para.nbin | isnan(binnum)) = [];
  end
  
  for i=1:length(binnum)
    firingmap(binnum(i)) = firingmap(binnum(i))+spikes(i);
  end
  norm_firingmap = firingmap./(bh.dwelltime/para.f);
  mean_firingmap = mean(norm_firingmap);
  
  for k=1:para.nbin
    if norm_firingmap(k)>0
      MI = MI + bh.norm_dwelltime(k) * (norm_firingmap(k)/mean_firingmap)...
              * log2(norm_firingmap(k)/mean_firingmap);
    end
  end
end
