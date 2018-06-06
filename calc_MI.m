%% should maintain positions, that are covered during longrunperiods
%% spike train with similar statistics should be applied to this one



function [MI, norm_firingmap] = calc_MI(bh,binloc,spike_times,spikes,T,ISI,para,shuffle)
  
  MI = 0;
  firingmap=zeros(1,para.nbin);
  
  if (nargin==8 && ~isempty(shuffle) && ~isempty(spike_times))
%      disp('shuffle')
    %% keep in mind: discontinuities in between long run periods
    %%% maybe rather feed single compartments in, one by one?
%      shuffled_spike_train = spike_shuffling('shift',true,spike_times,spikes,T,ISI,5);
%      shuffled_spike_train = spike_shuffling('dither',true,spike_times,spikes,T,ISI,5);
    shuffled_spike_train = spike_shuffling(shuffle,true,spike_times,spikes,T,ISI,2*para.f);
    
    spike_times = find(shuffled_spike_train);
    spikes = shuffled_spike_train(spike_times);
  end
  
  binloc = binloc(spike_times);
  for i=1:length(binloc)
    firingmap(binloc(i)) = firingmap(binloc(i))+spikes(i);
  end
  norm_firingmap = firingmap./bh.dwelltime;
  mean_firingmap = nanmean(norm_firingmap);        %% this shouldnt be the "cleared" firingmap, but the overall mean value during running events
  
  MI_arr = zeros(1,para.nbin);
  for k=1:para.nbin
    if bh.dwelltime(k) < 1
      MI_arr(k) = NaN;
    elseif norm_firingmap(k)>0
      MI_arr(k) = bh.norm_dwelltime(k) * (norm_firingmap(k)/mean_firingmap)...
            * log2(norm_firingmap(k)/mean_firingmap);
      MI = MI + bh.norm_dwelltime(k) * (norm_firingmap(k)/mean_firingmap)...
            * log2(norm_firingmap(k)/mean_firingmap);
    end
  end
  
%  %    if ~shuffle
%      figure
%      
%      subplot(3,1,1)
%      hold on
%      bar(firingmap/mean(firingmap),'r')
%      bar(-bh.dwelltime,'FaceColor','b')
%      ylabel('dwelltime / firing map')
%      title(sprintf('mean firing rate: %5.3g',mean_firingmap))
%      ylim([-5,20])
%      
%      subplot(3,1,2)
%      hold on
%      bar(norm_firingmap/mean_firingmap,'r')
%      bar(-bh.dwelltime,'FaceColor','b')
%      ylabel('dwelltime / firing map')
%      title(sprintf('mean firing rate: %5.3g',mean_firingmap))
%      ylim([-5,20])
%      
%      subplot(3,1,3)
%      plot(MI_arr)
%      ylabel('MI')
%      title(sprintf('overall MI: %5.3g vs %5.3g',MI, nansum(MI_arr)))
%      ylim([-0.2,2])
%      
%  %      subplot(3,1,3)
%  %      bar(bh.dwelltime,'FaceColor','b')
%  %      ylabel('dwelltime')
%      
%      waitforbuttonpress;
%      
%  %    end
%    
end