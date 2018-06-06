function [longwalk longrest firstlongwalkind firstlongrestind longwalkind...
    longrestind diffwalkind diffrestind fraclongwalk fraclongrest] = longrunrest_identifier(runrest,position,nframe)
%
% This program identifies long walk and rest periods from run/rest vector.
%
% INPUT ARGUMENT
% runrest: vector of run/rest state resampled at 15Hz
% position: vector of animal position resampled at 15Hz
% nframe:total number of images, usually 8989
%
% OUTPUT ARGUMENT
% longwalk:
% longrest:
% firstlongwalkind: frame indices of the beginning of long walk episodes
% firstlonngrestind: frame indices of the beginning of long rest episodes
% longwalkind: frame indices for long walk
% longrestind: frame indices for long rest
% diffwalkind:
% diffrestind: 
% fraclongwalk: fraction of total time calculated as long walk
% fraclongrest: fraction of total time calculated as long rest
%
% reviewed 11-15-2014

minleng = 30;
longwalkind2 = find(bwareaopen(runrest,minleng));

scrsz = get(0,'ScreenSize');

%% defining walk periods

walkind=find(runrest==1);% find walk data points
walkind_ori=walkind;

if runrest(1)==0
    walkind=[0 walkind];
end

if runrest(nframe)==0
    walkind=[walkind nframe+1];
end
diffwalkind=diff(walkind); 
ncontwalk=find(diffwalkind >2); %find gaps between walk periods

%% define rest periods
restind=find(runrest==0);% find rest data points
restind_ori=restind;

if runrest(1)==1
    restind=[0 restind];
end

if runrest(nframe)==1
    restind=[restind nframe+1];
end

diffrestind=diff(restind);
ncontrest=find(diffrestind >2);%find gaps between rest periods

%% extracting long walk and rest periods

minleng=30;% minimum continuity (frames) 30 = 2 sec 

%% long walk
longwalk=find(diffrestind >= minleng);% find gaps between rests (i.e. walk) >= minleng
firstlongwalkind=restind(longwalk)+1;% one point after 
%% reconstracting long walk indices
longwalkind=[];
for LP1=1:length(firstlongwalkind)
    longwalkind=[longwalkind firstlongwalkind(LP1):firstlongwalkind(LP1)+diffrestind(longwalk(LP1))-2];
end

%% long rest
longrest=find(diffwalkind >= minleng);
firstlongrestind=walkind(longrest)+1;
%% reconstracting long rest indices
longrestind=[];
for LP2=1:length(firstlongrestind)
    longrestind=[longrestind firstlongrestind(LP2):firstlongrestind(LP2)+diffwalkind(longrest(LP2))-2];
end


% longwalkind=longwalkind;
% longrestind=longrestind;

% longwalkind=walkind;
% longrestind=restind;

fraclongwalk=length(longwalkind)/nframe*100;
fraclongrest=length(longrestind)/nframe*100;

disp('---')
disp(sprintf('fraction of long walk is %4.1f', fraclongwalk))
disp(sprintf('fraction of long rest is %4.1f', fraclongrest))
disp('---')

longrunperiod=zeros(1,nframe);
longrunperiod(longwalkind)=1;

longrestperiod=zeros(1,nframe);
longrestperiod(longrestind)=1;

shortrunperiod=zeros(1,nframe);
shortrunperiod(setdiff(walkind_ori,longwalkind))=1;

shortrestperiod=zeros(1,nframe);
shortrestperiod(setdiff(restind_ori,longrestind))=1;

figure('Position',[scrsz(3)*0.1 scrsz(4)*0.1 scrsz(3)*0.8 scrsz(4)*0.8])

subplot(6,1,1);plot(position);title('position-rs');xlim([1 nframe]);ylim ([-1700 100])
subplot(6,1,2);bar(runrest);title('runrest-rs');xlim([1 nframe])
subplot(6,1,3);bar(longrunperiod,'r');title('longrun');xlim([1 nframe])
subplot(6,1,4);bar(shortrunperiod,'r');title('shortrun');xlim([1 nframe])
subplot(6,1,5);bar(longrestperiod,'g');title('longrest');xlim([1 nframe])
subplot(6,1,6);bar(shortrestperiod,'g');title('shortrest');xlim([1 nframe])
xlabel('frame number')

% figure
% all_rest_length=diffwalkind(ncontwalk)
% hist(all_rest_length,100)
% title('distribution of rest duration');xlabel('rest duration in # of frames');ylabel('# of events')

end

