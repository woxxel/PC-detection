
function [clusters] = prep_cluster(clusters_old,data,thr,pathData,suffix,loadData)
  
  disp('preparing clusters for PC field detection...')
  
  pathSave = pathcat(pathData,sprintf('clusters%s.mat',suffix));
  if exist(pathSave,'file') && loadData
    disp(sprintf('Clusters file %s already exists. Loading data...',pathSave))
    load(pathSave)
    return
  end
  
  tic
  counts = [clusters_old.ct];
  scores = [clusters_old.score];
  
  keep = find(counts >= thr(1) & scores >= thr(2));
  clusters_old = clusters_old(keep);
  
  nC = length(clusters_old);
  nSes = length(clusters_old(1).list);
  
  clusters(nC,nSes) = struct;
  
  for s = 1:nSes
    disp(sprintf('loading / reading session %d',s))
    pathSession = pathcat(pathData,sprintf('Session%02d',s));
    load(pathcat(pathSession,'resultsCNMF_MF1_LK1.mat'),'C2','S2');
    
    %%% load imaging data
    for c=1:nC
      n = clusters_old(c).list(s);
      if n
        clusters(c,s).ROI_ID = n;
        clusters(c,s).CaTrace = single(C2(n,:));
        clusters(c,s).S = single(S2(n,:));
        clusters(c,s).A = data(s).A(:,n);
        clusters(c,s).centroid = data(s).centroid(n,:);
      else
        clusters(c,s).ROI_ID = NaN;
      end
    end
    
  end
  toc
  
  save(pathSave,'clusters','-v7.3')
  disp(sprintf('data saved under %s',pathSave))
  
end