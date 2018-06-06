

function [PC_fields,clusters,bh,para] = detect_PC(basePath,mouse,p_thr,thr,repnum,clusters)
  
  loadData = true;
  
  pathData = sprintf('%s%d',basePath,mouse)
  loadPath = pathcat(pathData,sprintf('ROI_final_p=%4.2g.mat',p_thr))
  load(loadPath)
  
  suffix = sprintf('_s=%d_p=%4.2g',thr(1),thr(2));
  
  if nargin < 6
    clusters = prep_clusters(ROI_cluster_final,ROI_data,thr,pathData,suffix,loadData);
  end
  
  para = set_paras(mouse);
  
  bh = prep_behaviour(para,pathData,loadData);
  PC_fields = PC_detection(clusters,bh,para,repnum,pathData,suffix,loadData);
  
%    plot_PC_fields(PC_fields,clusters,para.nSes);
  
end