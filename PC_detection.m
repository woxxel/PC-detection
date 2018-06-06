

function [PC_fields] = PC_detection(clusters,bh,para,repnum,pathData,suffix,loadData)
  
  disp('detecting place fields...')
  
  pathSave = pathcat(pathData,sprintf('place_fields%s_nrep=%d.mat',suffix,repnum));
  if exist(pathSave,'file') && loadData
    disp(sprintf('PC field file %s already exists. Loading data...',pathSave))
    load(pathSave)
    return
  end
  
  nSes = para.nSes;
  prc = 20;               %% percentile for computing noise level of S
  nsd = 3;                %% number of std above noise level to be considered "activity"
  
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
    
    PC_fields(c).MI_rand = zeros(nSes,1);
    PC_fields(c).MI_frac(:) = NaN;
  end
  
%    bh(nSes) = struct;
  
  for s=1:nSes
    disp(sprintf('Processing session %d',s))
    
    tic
    %%% load data
    pathSession = pathcat(pathData,sprintf('Session%02d',s));
    
    %% set parameters
    IDs = [clusters(:,s).ROI_ID];
    detected = ~isnan(IDs);
    idx_clusters = find(detected);
    IDs = IDs(detected);
    ncell = sum(detected);
    
    location = bh(s).location(bh(s).longrunperiod);
    binloc = floor(location/para.binwidth)+1;
    binloc = min(binloc,para.nbin);
    %% check not really needed
    if any(binloc>para.nbin | isnan(binloc))
      disp('somethings wrong')
      binloc(binloc>para.nbin | isnan(binloc)) = [];
    end
    
    %% place cell detection from mutual information
    MI=zeros(1,ncell);
    status_tmp = zeros(ncell,1);
    for c=idx_clusters
        
        modeS = prctile(clusters(c,s).S(clusters(c,s).S>0),prc);                          %% get mode from overall activity
        activity = floor(sqrt(clusters(c,s).S(bh(s).longrunperiod)/(modeS*nsd)));         %% only activity from actual times
        %%% store, where longrunperiods start and shuffle only within those
        
        spike_times = find(activity);
        spikes = activity(spike_times);
        ISI = diff(spike_times);
        T = length(activity);
        
        [MI, PC_fields(c).firingmap(s,:)] = calc_MI(bh(s),binloc,spike_times,spikes,T,ISI,para);
        %% calculate MI from random data 
        MI_rand=zeros(1,repnum);
        for L=1:repnum
        
            [MI_rand(L), ~] = calc_MI(bh(s),binloc,spike_times,spikes,T,ISI,para,'dithershift');
%              [MI_rand(L), ~, ~] = calc_MI(bh(s),activity_rand,para,true);
        end
        
        status_tmp(c) = MI > prctile(MI_rand,95);    %% comparison with 95 percentile
        PC_fields(c).status(s) = status_tmp(c);
        
        PC_fields(c).MI(s) = MI;
        PC_fields(c).MI_frac(s) = MI / prctile(MI_rand,95);
        PC_fields(c).MI_rand = MI_rand;
        
        [PC_fields(c).max_fr(s),PC_fields(c).max_pos(s)] = max(PC_fields(c).firingmap(s,:));
        
%          figure
%          hold on
%          histogram(MI_rand)
%          plot(MI,0,'ro','MarkerSize',10)
%          hold off
%          figure
%          subplot(2,1,1)
%          bar(bh(s).dwelltime)
%          
%          subplot(2,1,2)
%          bar(PC_fields(c).firingmap(s,:))
%  %          
%          waitforbuttonpress;
        
        
    end % c
    
    toc
    idx_PC = find(status_tmp);
    pc_ct = length(idx_PC);
    disp(sprintf('found %d place cells from %d neurons',pc_ct,ncell))
  end
  
  save(pathSave,'PC_fields','-v7.3')
  disp(sprintf('Place field data saved under %s',pathSave))

end
