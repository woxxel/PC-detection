

function test_MI(PC_fields,clusters,bh,para)
  
  prc = 20;               %% percentile for computing noise level of S
  nsd = 3;                %% number of std above noise level to be considered "activity"
  
  stat_tmp = [PC_fields.status];
  
  %% get random trace (and output index)
  [s_ind, c_ind] = find(~isnan(stat_tmp));
  idx_tmp = randi(length(s_ind));
  
  s = s_ind(idx_tmp)
  c = c_ind(idx_tmp)
  
  %% get base MI
  modeS = prctile(clusters(c,s).S(clusters(c,s).S>0),prc);
  activity = floor(sqrt(clusters(c,s).S/(modeS*nsd)));
  
  [MI_base, ~] = calc_MI(bh(s),activity,para);
  
%    longrunspikes = activity & bh(s).longrunperiod;
%    spikes = activity(longrunspikes);
  binnum = floor(bh(s).location(bh(s).longrunperiod)/para.binwidth)+1;
  binnum = min(80,binnum);
  run_pos = zeros(1,para.nbin);
  for i=1:length(binnum)
    run_pos(binnum(i)) = run_pos(binnum(i)) + 1;
  end
  
  %% do all possible shifts and store MI(shift)
  nframe = length(clusters(c,s).S);
  MI_shift = zeros(1,nframe);
  
  for shift = 1:nframe
    activity_rand = [activity(shift:nframe),activity(1:shift-1)];
    
    [MI_shift(shift), ~] = calc_MI(bh(s),activity_rand,para);
  end
  
  
  %% calculate and display histogram of MI with some "nearby shifts" not included
  
  figure('position',[100 100 1200 900])
  
  subplot(2,2,3)
%    bar(PC_fields(c).dwellmap(s,:),'FaceColor','r')
  bar(bh(s).norm_dwelltime,'FaceColor','r')
  
%    subplot(3,2,3)
%    bar(run_pos,'FaceColor','y')
  
  subplot(2,2,1)
  bar(PC_fields(c).firingmap(s,:),'FaceColor','b')
  
  
  half_idx = floor(nframe/2);
  
  MI_shift = [MI_shift(half_idx:end),MI_shift(1:half_idx-1)];
  
  subplot(1,2,2)
  hold on
  plot(-half_idx:half_idx,MI_shift,'k.')
  plot(0,MI_base,'ro')
  plot(0,prctile(MI_shift,95),'rx')
  hold off
  xlabel('shift')
  ylabel('MI(shift)')

end