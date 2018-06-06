%%% function to test the method of obtaining mutual information from Skaggs formula and comparing it to MI from shuffled datasets

%%% Inputs:
%%%       - location:   location dataset from experiment
%%%       - nSpikes:    numbers of spikes to be distributed (rather rate?)
%%%       - place:      number of bin for which the cell should code. 0 if non-coding (random)
%%%       - nShuffle:   number of shuffles to obtain MI-histogram

function test_Skaggs(dwelltime,para,peakRate,place,baseRate,nShuffle)
  
  %%% construct artificial data
  %%% how to get locations???
  spikes = zeros(1,length(dwelltime));
  nSpikes = peakRate*sum(dwelltime(place));
  nSpikes
  for s=1:nSpikes
    p = randi(length(place));
    pos = round(normrnd(place(p),1));     %% place field should be ~4-5cm wide?
    spikes(pos) = spikes(pos) + 1;
  end
  
  for s=1:baseRate*sum(dwelltime)
    pos = randi(para.nbin);               %% uniform noise
    spikes(pos) = spikes(pos) + 1;
  end
  
  figure
  subplot(2,1,1)
  bar(spikes,'FaceColor','r')
  
  subplot(2,1,2)
  bar(dwelltime,'FaceColor','b')
  title(sprintf('total time measured: %5.3g',sum(dwelltime)))
  
  
  %%% compute MI from this
  
  
end