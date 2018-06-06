

function post_PC_fields()
  
  
  for c = 1:nC
    for s = 1:nSes
      
      if PC_fields(c).status 
        %% smooth firing rate
        sigma = 2;
        smoothed_map = imgaussfilt(PC_fields(c).firingmap(s,:),sigma);
        fr_mean = mean(smoothed_map);
        
        fields = smoothed_map > fr_mean;
        fields = bwareaopen(fields,3);
        % apply
        
        %% find regions above average firing rate (for a few bins)
        %% and mark those as coding area
      
      
    end
  end
  
end