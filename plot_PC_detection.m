

function plot_PC_detection(PC_fields,clusters,c)
  
  nSes = length(PC_fields(c).status);
  
  figure('position',[200 200 1200 900])
  for s = 1:nSes
    subplot(6,3,s)
    hold on
    bar(PC_fields(c).firingmap(s,:))
    plot([PC_fields(c).max_pos(s) PC_fields(c).max_pos(s)],[0,10],'r--')
    title(sprintf('status: %d, MI: %5.3g',PC_fields(c).status(s),PC_fields(c).MI(s)))
  end
  
%    subplot(6,2,11)
  %behaviour
%    bar(bh(s).dwelltime)
  
  subplot(6,2,12)
  s_arr = linspace(1,15,15);
  mask = find(PC_fields(c).status);
  plot(PC_fields(c).max_pos(mask),s_arr(mask),'kx')
  xlim([0,80])
  ylim([0,16])
  xlabel('Location')
  ylabel('session')
  suptitle('Firing fields')
  
  
  figure('position',[500 500 800 600])
  hold on
  for s = 1:nSes
    if ~isnan(clusters(c,s).ROI_ID)
      col = 1-[s/nSes,s/nSes,s/nSes];
      contour(reshape(clusters(c,s).A,512,512),[0.1 0.1]*max(clusters(c,s).A),'LineColor',col)
      centr = [clusters(c,s).centroid(1),clusters(c,s).centroid(2)];
    end
  end
  hold off
  
  xlim([centr(2)-15,centr(2)+15])
  ylim([centr(1)-15,centr(1)+15])
  