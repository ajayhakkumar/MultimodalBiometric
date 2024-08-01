function fingerseg


% the function read the fingerprint images from dataset/dataset1
% ,enhance the image using ridge orientation
% segment the region of interest (fingerprint region) 
% and store it into mdbase/fingerprint directory


 clc;clear all;
 disp('Finger print Segmentation and Enhancement')
 cpcs = cell(60,1);
 
 for d=1 : 60  
    dname= num2str(d);  
    cpcs{d}=fingerenc(dname);
 end

save('corepoints.mat','cpcs');
  
  
function cpc=fingerenc(dname)

source1 = 'Dataset';
source2 ='fingerprint';

fpath = sprintf('%s\\%s\\%s\\*.tif',source1,source2,dname)
flst=dir(fpath);
mkdir(sprintf('mdbase\\%s',source2),dname);

cpc = cell(length(flst),2);
for n=1 : length(flst)
 fpath2 = sprintf('%s\\%s',flst(n).folder,flst(n).name)
 [sim,cp]=fingerenhance2(fpath2);    
 cpc{n,1} =flst(n).name;  
 cpc{n,2} =cp; 
 dstname = sprintf('mdbase\\%s\\%s\\%s',source2,dname,flst(n).name);
 imwrite(sim,dstname);   
end
