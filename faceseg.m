function faceseg

% the function read the face images from dataset/dataset1
% ,localizae face on the image using face segmentation
% segment the region of interest (face region) 
% and store it into mdbase/face directory


clc; clear all;

detector1 = vision.CascadeObjectDetector('FrontalFaceCART');
detector2 = vision.CascadeObjectDetector('ProfileFace');
detector1.MinSize = [160 160];
detector1.MaxSize = [260 300];
detector2.MinSize = [160 160];
detector2.MaxSize = [260 260];

for d = 1: 60
   dname = num2str(d);   
   fsegment(detector1,detector2,dname);
end



function fsegment(detector1,detector2,dname)


source1 = 'Dataset\database1';
source2 = 'face';

fpath = sprintf('%s\\%s\\%s\\*.jpg',source1,dname,source2);
flst=dir(fpath);


mkdir(sprintf('mdbase\\%s',source2),dname);


for n=1 : length(flst)
    
 fpath = sprintf('%s\\%s',flst(n).folder,flst(n).name);
 dstname = sprintf('mdbase\\%s\\%s\\%s',source2,dname,flst(n).name);
 
 img = imread(fpath);
 bbox = step(detector1,img);
if ( ~isempty(bbox)  )
         
     if (size(bbox,1)==1)                   
       facerg = imcrop(img,bbox);
       sim=rgb2gray(facerg);
       imwrite(sim,dstname);                    
     end

else
     bbox = step(detector2,img);      
      if (size(bbox,1)==1)          
        facerg = imcrop(img,bbox);
        sim=rgb2gray(facerg);
        imwrite(sim,dstname); 
     end
             
end


end

