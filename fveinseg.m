function fveinseg


% the function read the face images from dataset/dataset1
% ,localizae fingervein on the image using thresholding segmentation
% segment the region of interest (fingervein region) 
% and store it into mdbase/fingervein directory


clc;clear all;
disp('Fingervein Segmentation')
     
 for d=1 : 60  
    dname= num2str(d);  
    fveinsegment(dname);
 end
 




function fveinsegment(dname)

source1 = 'Dataset\database1';
source2 = 'fingervein';

fpath = sprintf('%s\\%s\\%s\\*.bmp',source1,dname,source2);
flst=dir(fpath);
mkdir('mdbase\\fingervein\\',dname);

for n=1 : length(flst)
 fpath = sprintf('%s\\%s',flst(n).folder,flst(n).name);
 
im = imread(fpath);
im = im2double(im);
[M,N,C]=size(im);
if( C>1)
 im = rgb2gray(im);    
end

sim=roiseg(im);
dstname = sprintf('mdbase\\fingervein\\%s\\%s',dname,flst(n).name)
imwrite(sim,dstname);


end


function sim=roiseg(im)
    
im = imgaussfilt(im,0.8);
level = graythresh(im);
BW = imbinarize(im,level);

[M,N]=size(BW);
K=true(4,N);
K(1:2,:)=false;

K2=false(4,N);
K2(1:2,:)=true;

S1 = zeros(M,1);
S2 = zeros(M,1);
for i=3 : M -2    
S1(i)=sum(sum( BW(i-2:i+1,:)== K ));         
S2(i)=sum(sum( BW(i-2:i+1,:)== K2 ));         
end

%figure; subplot(2,2,1);imshow(BW);
%subplot(2,2,3);plot(1:M,S1,'*-b');
%subplot(2,2,4);plot(1:M,S2,'*-r');

[v1,m1]=max(S1);

while( m1>100 )
S1(m1)=0;    
[v2,m1]=max(S1);    
end

[v2,m2]=max(S2);
while( (m2-m1)<50 )
 S2(m2)=0;    
[v2,m2]=max(S2);     
end
sim = im(m1:m2,:);
