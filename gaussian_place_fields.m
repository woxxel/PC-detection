function gaussian_place_fields(datafolder)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


cd(datafolder)
close all

datafolder

pcfile=dir('*_pc_rep1000*0.1_index.txt');
firingmapfile=dir('*_pc_rep1000*0.1_firingmap.txt');
pfpeakfile=dir('*_pc_rep1000*0.1_peak.txt');
nbin=80;
scrsz = get(0,'ScreenSize');

pc_ind=importdata(pcfile.name);
pc_num=length(pc_ind)
pfpeak=importdata(pfpeakfile.name);
pc_pfpeak=pfpeak(pc_ind);
firingmap=importdata(firingmapfile.name);

pc_firingmap=firingmap(pc_ind,:);
[peakheight peakheight_bin]=max(pc_firingmap,[],2);
[sorted_peakheight_bin sorted_peakheight_bin_ind]=sort(peakheight_bin);

smfactor=5;
gauss_filter=gausswin(smfactor);

for i=1:pc_num
norm_pc_firingmap(i,:)=pc_firingmap(i,:)./peakheight(i);
tri_norm_pc_firingmap=[norm_pc_firingmap(i,:) norm_pc_firingmap(i,:) norm_pc_firingmap(i,:)];
tri_norm_pc_firingmap=filtfilt(gauss_filter,smfactor,tri_norm_pc_firingmap);
gauss_norm_pc_firingmap(i,:)=tri_norm_pc_firingmap(81:160);
gauss_norm_pc_firingmap(i,:)=gauss_norm_pc_firingmap(i,:)./max(gauss_norm_pc_firingmap(i,:));
end

for j=1:pc_num
sorted_norm_pc_firingmap(j,:)=norm_pc_firingmap(sorted_peakheight_bin_ind(j),:);
tri_sorted_norm_pc_firingmap=[sorted_norm_pc_firingmap(j,:) sorted_norm_pc_firingmap(j,:) sorted_norm_pc_firingmap(j,:)];
tri_sorted_norm_pc_firingmap=filtfilt(gauss_filter,smfactor,tri_sorted_norm_pc_firingmap);
gauss_sorted_norm_pc_firingmap(j,:)=tri_sorted_norm_pc_firingmap(81:160);
gauss_sorted_norm_pc_firingmap(j,:)=gauss_sorted_norm_pc_firingmap(j,:)./max(gauss_sorted_norm_pc_firingmap(j,:));
end

[peakheight_gauss peakheight_gauss_bin]=max(gauss_norm_pc_firingmap,[],2);

pfdist_xvec=0.5:1:nbin+0.5;
pfdist_elem=histc(peakheight_bin,pfdist_xvec);
pfdist_elem_gauss=histc(peakheight_gauss_bin,pfdist_xvec);

dirname=savedata(datafolder, ['gausspc' num2str(smfactor)]);
        
dlmwrite([datafolder '_gausspc' num2str(smfactor) '_rep1000_0.1_peak.txt'], peakheight_gauss_bin, '\t');
dlmwrite([datafolder '_gausspc' num2str(smfactor) '_rep1000_0.1_firingmap.txt'], gauss_norm_pc_firingmap, '\t');



cd ..


end

