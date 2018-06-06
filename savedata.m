function  [dirname] = savedata(datafolder,ext)
% this program saves figures in .fig and .pngÅ@format 
% 
%INPUT ARGUMENT
% datafoler: original data folder name, used to create new folder name
% ext: extention to attach when create new folder name 
%
%OUTPUT ARGUMENT 
% dirname: name of new folder
%
% reviewed 11-15-2014

%% create new folder name
savefile=[datafolder '_' ext];               
dirname=[savefile,'_figs'];

if exist(dirname, 'dir') ==0
    mkdir (dirname)
end

cd (dirname)

%% save .fig and .png figures
disp('Saving figures...')

h = findobj('Type','figure');

for  i=1:length(h)
    savename=[savefile,'_Fig',num2str(length(h)-i+1),'.fig'];    
    saveas(h(i),savename,'fig');
end
     
for j=1:length(h)
  savename1=[savefile,'_Fig',num2str(length(h)-j+1),'.png'];
  set(h(j),'PaperPositionMode','auto')
  print(h(j),'-r0','-dpng',savename1);
end         

cd ..

end

