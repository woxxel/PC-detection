%%% plotting all kind of results

%% what about recruitment? check, where from and to receptive fields are changed

function plot_PC_fields(PC_fields,clusters,bh,para,nSes)
  
  nC = size(clusters,1);
  
  scrsz = get(0,'ScreenSize');
  
  fr_tmp = [PC_fields.max_fr];
  MI_tmp = [PC_fields.MI];
  stat_tmp = [PC_fields.status];
  pos_tmp = [PC_fields.max_pos];
  
  for s=1:nSes
    IDs = [clusters(:,s).ROI_ID];
    detected = ~isnan(IDs);
    ncell(s) = sum(detected);    
  end
  PC_occurence = sum(stat_tmp,1);
  ncell_PC = sum(stat_tmp,2);
  
  figure('position',[100 100 1500 600])
  subplot(2,2,1)
  tight_axes(gca)
  histogram(PC_occurence(PC_occurence>0),linspace(0,16,17))
%    set(gca, 'YScale', 'log')
  xlabel('# sessions')
  ylabel('# of cells')
  title(sprintf('total number of PCs: %d',nnz(sum(stat_tmp,1))))
  
  
  %%% display # and ratio of PC per session
  subplot(2,2,2)
  tight_axes(gca)
  hold on
  plot(para.t_s,ncell,'DisplayName','Total number of cells')
  plot(para.t_s,ncell_PC,'DisplayName','Number of place cells')
  hold off
  xlabel('time [h]')
  ylim([0 1000])
  legend('Location','NorthWest')
  
  ax = axes('position',[0.85,0.8,0.12,0.15]);
  tight_axes(ax)
  plot(para.t_s,ncell_PC'./ncell)
  xlabel('time [h]')
  ylabel('ratio PCs')
  
  
  %%% add recall etc here
  
  %%% common # of PC per adjacent sessions
  recall = zeros(nSes,para.t_s(end));
  recall_right = zeros(nSes,para.t_s(end));
  
  stable_bin = zeros(80,1);
  tolerance = 10/1.25;
  n_total = zeros(para.t_s(end),1);
  
%    for d = 1:nSes-1
%      n_total(d) = sum(ncell_PC(1:(nSes-d)));
%    end
  for s = 1:nSes-1
    d = 1;
    while s+d <= nSes
      dt = para.t_s(s+d)-para.t_s(s);
      n_total(dt) = n_total(dt) + ncell_PC(s);
      d=d+1;
    end
  end
  
  for c=1:nC
    for s = 1:nSes-1
      d = 1;
      while s+d <= nSes
        dt = para.t_s(s+d)-para.t_s(s);
        if stat_tmp(s,c) && stat_tmp(s+d,c)
          recall(s,dt) = recall(s,dt) + 1;
          
          if abs(PC_fields(c).max_pos(s) - PC_fields(c).max_pos(s+d)) < tolerance
            recall_right(s,dt) = recall_right(s,dt) + 1;
            stable_bin(PC_fields(c).max_pos(s)) = stable_bin(PC_fields(c).max_pos(s)) + 1;
          end
        end
        d = d+1;
      end
    end
  end
  
  hours = 1:para.t_s(end);
  
  tot_recall = sum(recall,1);
  tot_recall_right = sum(recall_right,1);
  
%    figure('position',[100 100 1500 600])
  subplot(2,3,4)
  tight_axes(gca)
  hold on
  plot(hours(tot_recall>0),tot_recall(tot_recall>0),'k','DisplayName','anywhere')
  plot(hours(tot_recall_right>0),tot_recall_right(tot_recall_right>0),'r','DisplayName','same location \pm 10cm)')
  hold off
  legend()
  xlabel('time [h]')
  ylabel('# Recalled PCs')
  plt_rec = tot_recall./n_total';
  plt_rec_right = tot_recall_right./n_total';
  
  %%% ratio retained PC location vs session difference
  subplot(2,3,5)
  tight_axes(gca)
  hold on
  plot(hours(tot_recall>0),plt_rec(tot_recall>0),'k','DisplayName','anywhere')
  plot(hours(tot_recall_right>0),plt_rec_right(tot_recall_right>0),'r','DisplayName','same location \pm 10cm')
  hold off
  xlabel('time [h]')
  ylabel('# Recalled PCs')
  legend()
  
  %%% PC stability vs location
  subplot(2,3,6)
  tight_axes(gca)
  bar(stable_bin)
  xlabel('Location')
  title('stable positions')
  
  
  
  
  %%% display MI vs FR of PC & nPC
  m_PFR = zeros(nSes,1);
  m_MI = zeros(nSes,1);
  
  figure('position',[100 100 1200 400])
  subplot(1,3,3)
  hold on
  for s = 1:nSes
    PC_ind = find(stat_tmp(s,:));
    scatter(fr_tmp(s,:),MI_tmp(s,:))
    scatter(fr_tmp(s,PC_ind),MI_tmp(s,PC_ind),'filled')
  end
  hold off
  xlim([0,20])
  xlabel('peak firing rate')
  ylabel('Mutual Information')
  
%    figure('position',[100 100 800 600])
  
  ax_fr = subplot(1,3,1);
  ax_MI = subplot(1,3,2);
  hold(ax_fr,'on')
  hold(ax_MI,'on')
  
  for s=1:nSes
    PC_ind = find(stat_tmp(s,:));
    m_PFR(s) = mean(fr_tmp(s,PC_ind));
    m_MI(s) = mean(MI_tmp(s,PC_ind));
    plot(ax_fr,para.t_s,fr_tmp(:,PC_ind),'kx')
    plot(ax_MI,para.t_s,MI_tmp(:,PC_ind),'rx')
  end
  plot(ax_fr,para.t_s,m_PFR,'k-','LineWidth',3)
  plot(ax_MI,para.t_s,m_MI,'r-','LineWidth',3)
  
  hold(ax_fr,'off')
  hold(ax_MI,'off')
  
  xlabel(ax_fr,'time [h]')
  xlabel(ax_MI,'time [h]')
  
  ylabel(ax_fr,'peak firing rate')
  ylabel(ax_MI,'Mutual Information')
  
  
  %%% coverage of place over all sessions
  figure('position',[100 100 scrsz(3)/2 scrsz(4)])
  for s=1:nSes
    [~, IX] = sort(pos_tmp(s,:));
    
    fm_plt = zeros(ncell_PC(s),80);
    i=0;
    for c=IX
      if PC_fields(c).status(s)
        i=i+1;
        fm_plt(i,:) = PC_fields(c).firingmap(s,:) / max(PC_fields(c).firingmap(s,:));
      end
    end
    
    
%      plot_pos = s+floor((s-1)/5)*5;% + floor((s-1)/5)*5
%      plot_pos2 = s+floor((s-1)/5+1)*5;% + floor((s-1)/5)*5
    
    ax = subplot(5,3,s);
    tight_axes(ax)
    hold on
    imagesc(fm_plt)
    bar(-bh(s).dwelltime*0.2*ncell_PC(s)/max(bh(s).dwelltime))
    plot([60 60],[-0.2*ncell_PC(s) ncell_PC(s)],'r--','LineWidth',2)
    plot([20 20],[-0.2*ncell_PC(s) ncell_PC(s)],'g--','LineWidth',2)
    hold off
    title(sprintf('Session %d, # PCs: %d',s,ncell_PC(s)))
    set(gca,'YTick',[])
    set(gca,'YDir','reverse')
    ylim([-0.2*ncell_PC(s) ncell_PC(s)])
    
%      plot_pos3 = s+floor((s-1)/5+2)*5 + floor((s-1)/5)*5;
%      subplot(9,5,plot_pos3)
%      bar(bh(s).dwelltime)
%      ylim([0,300])
  end

  
  
  
  
  %%% color coded location per session
  
%    stat_tmp(stat_tmp==0) = nan;
  color_pos = stat_tmp.*pos_tmp;
  color_pos(~stat_tmp) = nan;
  color_pos(end+1,:) = 1:size(color_pos,2);
%    color_pos
  
  ordered_fields = sort_PC_fields(color_pos);
%    size(ordered_fields)
%    ordered_fields
  
  new_fields = zeros(para.nSes,nC);
  for i = 1:nC
    new_fields(:,i) = color_pos(1:end-1,ordered_fields(i));
  end
    
  %%% ordered by...
  %% first occurence
  %% decision tree to order stuff
  %% iterative call of function
  
  figure
  imagesc(new_fields)
  colorbar
  
  
  %%% not so telling...
  remap_hist = zeros(para.nbin+1);
  
  for s = 1:nSes-1
    for c = 1:nC
      if isnan(color_pos(s,c))
        idx_y = para.nbin + 1;
      else
        idx_y = color_pos(s,c);
      end
      if isnan(color_pos(s+1,c))
        idx_x = para.nbin + 1;
      else
        idx_x = color_pos(s+1,c);
      end
%        if ~(isnan(color_pos(s,c)) || isnan(color_pos(s+1,c)))
      remap_hist(idx_y,idx_x) = remap_hist(idx_y,idx_x) + 1;
%        end
    end
  end
  
  figure
  bar3(remap_hist)
  xlabel('s+1')
  ylabel('s')
  
  
  
  
  %%% take into account 2nd peak as well!!
  
  idx_others = ones(1,para.nbin+1);
  idx_others(15:25) = 0;
  idx_others(55:65) = 0;
  idx_others(end) = 0;
  
  idx_others = find(idx_others);
  
  gt_any = sum(remap_hist(15:25,1:end-1));
  rw_any = sum(remap_hist(55:65,1:end-1));
  others_any = sum(remap_hist(idx_others,1:end-1));
  nPC_any = remap_hist(end,1:end-1);
  
  
  any_gt = sum(remap_hist(1:end-1,15:25),2);
  any_rw = sum(remap_hist(1:end-1,55:65),2);
  any_others = sum(remap_hist(1:end-1,idx_others),2);
  any_nPC = remap_hist(1:end-1,end);
  
  
  
  figure
  subplot(4,2,1)
  bar(gt_any)
  title('gate -> any')
  subplot(4,2,2)
  bar(any_gt)
  title('any -> gate')
  
  subplot(4,2,3)
  bar(rw_any)
  title('reward -> any')
  subplot(4,2,4)
  bar(any_rw)
  title('any -> reward')
  
  subplot(4,2,5)
  bar(others_any)
  title('others -> any')
  subplot(4,2,6)
  bar(any_others)
  title('any -> others')
  
  subplot(4,2,7)
  bar(nPC_any)
  title('nPC -> any')
  xlabel('bins')
  
  subplot(4,2,8)
  bar(any_nPC)
  title('any -> nPC')
  xlabel('bins')
  
  
  
  figure('position',[300 300 1500 1200])
  subplot(1,2,1)
  bar([gt_any;rw_any;others_any]','stacked')
  title('... -> any')
  
  subplot(1,2,2)
  bar([any_gt';any_rw';any_others']','stacked')
  title('any -> ...')
%    subplot(4,2,7)
%    bar(sum(remap_hist(1:end-1,end)))
%    title('any -> nPC')
  
%    PC_stab = sum(stat_tmp,1);
%    c_idx = find(PC_stab>6);
%    c = c_idx(2);
%    nSes = length(PC_fields(c).status);
%    
%    figure('position',[200 200 1200 900])
%    for s = 1:nSes
%      subplot(6,3,s)
%      hold on
%      bar(PC_fields(c).firingmap(s,:))
%      plot([PC_fields(c).max_pos(s) PC_fields(c).max_pos(s)],[0,10],'r--')
%      title(sprintf('status: %d, MI: %5.3g',PC_fields(c).status(s),PC_fields(c).MI(s)))
%    end
%    
%  %    subplot(6,2,11)
%    %behaviour
%  %    bar(bh(s).dwelltime)
%    
%    subplot(6,2,12)
%    s_arr = linspace(1,15,15);
%    mask = find(PC_fields(c).status);
%    plot(PC_fields(c).max_pos(mask),s_arr(mask),'kx')
%    xlim([0,80])
%    ylim([0,16])
%    xlabel('Location')
%    ylabel('session')
%    suptitle('Firing fields')
%    
%    
%    figure('position',[500 500 800 600])
%    hold on
%    for s = 1:nSes
%      if ~isnan(clusters(c,s).ROI_ID)
%        col = 1-[s/nSes,s/nSes,s/nSes];
%        contour(reshape(clusters(c,s).A,512,512),[0.1 0.1]*max(clusters(c,s).A),'LineColor',col)
%        centr = [clusters(c,s).centroid(1),clusters(c,s).centroid(2)];
%      end
%    end
%    hold off
%    
%    xlim([centr(2)-15,centr(2)+15])
%    ylim([centr(1)-15,centr(1)+15])
end



function [ordered_entries] = sort_PC_fields(field_array)
  
  first_entries_idx = find(field_array(1,:)>0);
  first_entries = field_array(end,first_entries_idx);
  
  last_entries_idx = find(isnan(field_array(1,:)));
  last_entries = field_array(end,last_entries_idx);
  
  if size(field_array,1) > 2
    if length(first_entries_idx)
      [ordered_entries_first] = sort_PC_fields(field_array(2:end,first_entries_idx));
    else
      ordered_entries_first = [];
    end
    if length(last_entries_idx)
      [ordered_entries_last] = sort_PC_fields(field_array(2:end,last_entries_idx));
    else
      ordered_entries_last = [];
    end
    ordered_entries = [ordered_entries_first ordered_entries_last];
%      first_entries = first_entries_lower & first_entries;
  else
    ordered_entries = [first_entries,last_entries];
  end
  
end





%  %% display firing maps of place cells sorted by peaks of firing maps
%  
%  for L=1:good_pcnum
%    sorted_place_map(L,:)=firingmaps(place_cell_order(L),:);
%    norm_sorted_place_map(L,:)=sorted_place_map(L,:)/max(sorted_place_map(L,:));  
%  end
%  
%  figure
%  imagesc(norm_sorted_place_map);xlabel('position'); ylabel('sorted place cells');
%  hold on
%      line([20.5 20.5],[1 good_pcnum],'Color','r','lineStyle',':','LineWidth',1)% reward site
%      line([60.5 60.5],[1 good_pcnum],'Color','r','lineStyle',':','LineWidth',1)% gate
%  hold off
%  
%  %% display raster plot in which place cells are sorted by their place field peaks
%  non_pc=setdiff(goodFRcells, good_sig_pc_ind);
%  reorg_binaryU=[binaryU(:,pfp_sorted_good_sig_pc_ind) binaryU(:,non_pc)];
%  
%  figure('Position',[scrsz(3)/4 0 scrsz(3)/2 scrsz(4)*0.8])%
%  subplot(8,1,1);plot(bh.location);axis ([1, nframe, -1700, 100])
%      hold on 
%      line([1 nframe],[-400 -400],'Color','r','lineStyle',':')% reward site
%      line([1 nframe],[-1200 -1200],'Color','r','lineStyle',':')% gate
%      ylabel('position');xlabel('time');title('position-rs')
%      hold off
%  subplot(8,1,[2 7]);imagesc(reorg_binaryU');colormap gray;
%  title('sorted stop cells and non-stop cells')
%  hold on
%  line([1 nframe],[good_pcnum good_pcnum],'Color','w')%
%  hold off
%  ylabel('cell #');xlabel('frame #')
%  
%  
%  %% histogram of MI distribution
%  xvec=linspace(0,10,41);
%  MI_elem=histc(MI(goodFRcells),xvec);
%  good_sig_MI_elem=histc(MI(good_sig_pc_ind),xvec);
%  
%  histmat=[(MI_elem-good_sig_MI_elem); good_sig_MI_elem];
%  figure
%  bar(xvec, histmat', 'stacked');xlim([-0.5 10.5]);xlabel('MI')
%  
%  %% histogram of peak firing rate distribution
%  PFR_elem=histc(PFR(goodFRcells),xvec);
%  good_sig_PFR_elem=histc(PFR(good_sig_pc_ind),xvec);
%  
%  histmat2=[(PFR_elem-good_sig_PFR_elem) good_sig_PFR_elem];
%  figure
%  bar(xvec, histmat2, 'stacked');xlim([-0.5 10.5]);xlabel('PFR')
%  
%  %%  histogram of place field distrobution
%  good_sig_pc_pfpeak=pfpeak(good_sig_pc_ind);
%  
%  figure
%  pfdist_xvec=0.5:1:nbin+0.5;
%  pfdist_elem=histc(good_sig_pc_pfpeak,pfdist_xvec);
%  bar(pfdist_xvec,pfdist_elem);xlim([min(pfdist_xvec) max(pfdist_xvec)])


%% saving data

%  dirname=savedata(datafolder, ['pc_rep1000_' num2str(sd)]);

%  dlmwrite([datafolder '_pc_rep1000_' num2str(sd) '_index.txt'], good_sig_pc_ind,'\t');
%  dlmwrite([datafolder '_pc_rep1000_' num2str(sd) '_peak.txt'], pfpeak, '\t');
%  dlmwrite([datafolder '_pc_rep1000_' num2str(sd) '_firingmap.txt'], norm_firingmap, '\t');

%  diary off
%  movefile(logfilename,dirname);

%  toc
%  
%  cd ..
%  











%      sd=0.1
    
    
    
%      binaryU=U;
%      binaryU(find(binaryU<sd))=0;
%      binaryU(find(binaryU>=sd))=1;
    
  %    bh.longrunperiod=zeros(1,nframe);
  %    bh.longrunperiod(longwalkind)=1;
    
      %% calcurate firing rate    
%        firingrate=sum(binaryU,1)/600;
%        mFR=mean(firingrate);
%        stdFR=std(firingrate);
      
%        highFRcells=find(firingrate > mFR+3*stdFR);% high FR cells > mean + 3SD
%        lowFRcells=find(firingrate < 0.1);%  low FR cells < 0.1 Hz
      
%        nhighFRcells=length(highFRcells)
%        nlowFRcells=length(lowFRcells)
%        
%        allcells=1:ncell;
%        abFRcells=union(highFRcells, lowFRcells);
%        goodFRcells=setdiff(allcells, abFRcells);
%        
%        goodFRCvec=zeros(1,ncell);
%        goodFRCvec(goodFRcells)=1;
%        
      
    
    