function mmfeatureExtr


%fingerprint- extract local-lbp features from enchanced fingerprint 

nuser = 60;

if ( exist('fingerprintfeatures.mat','file')~=2 )
   
disp('Starts Features Extraction from finger print image')    

fingerX =cell(nuser,1);fingerY=cell(nuser,1);

dbase='mdbase\fingerprint\';
for d = 1 : nuser
  dname = d;     
  fpath=sprintf('%s%s',dbase,num2str(dname))
  fX=fingerfeature22(dbase,dname);    
  fingerX{d} =fX;    
  ns=size(fX,1);    
  fingerY{d}= ones(ns,1).*d;       
  
end
 
 
save('fingerprintfeatures.mat','fingerX','fingerY');
disp('Features Extraction ends, and saved on fingerprintfeatures.mat')    

end

%Face recognition- extract LBP features from segmented face

if ( exist('facefeatures.mat','file')~=2 )
   
disp('Starts Features Extraction from face image')        
faceX =cell(nuser,1);faceY=cell(nuser,1);

source='mdbase\face';
for d = 1 : nuser
  
  dname = num2str(d);     
  fpath=sprintf('%s\\%s\\*.jpg',source,dname)    
  flst=dir(fpath);  
  fX=[];fY=[];
  for n=1 : length(flst)
    fpath = sprintf('%s\\%s',flst(n).folder,flst(n).name); 
    img = imread(fpath);
    img = imresize(img,[240,240]);
    lbpFeatures = extractLBPFeatures(img,'CellSize',[32 32],'Normalization','None');        
    fX(n,:) = lbpFeatures;
    fY=[fY; d];
  end  
  faceX{d} = fX; 
  faceY{d} = fY;
  
end
save('facefeatures.mat','faceX','faceY');
disp('Face Features Extraction end and  saved on facefeatures.mat')        

end


%fingervein - extract LBP features from  segmented fingervein

if ( exist('fveinfeatures.mat','file')~=2 )
   
fveinX =cell(nuser,1);fveinY =cell(nuser,1);
disp('Fingervein Features Extraction starts');
source='mdbase\fingervein';
for d = 1 : nuser  
  dname = num2str(d);     
  fpath=sprintf('%s\\%s\\*.bmp',source,dname)     
  flst=dir(fpath);  
  fX=[];fY=[];
  for n=1 : length(flst)
    fpath = sprintf('%s\\%s',flst(n).folder,flst(n).name) 
    img = imread(fpath);
    img=imcrop(img,[1 1 320 120]);
    lbpFeatures = extractLBPFeatures(img,'CellSize',[32 32],'Normalization','None');            
    fX(n,:) = lbpFeatures;
    fY(n)= d;
  end  
  fveinX{d} = fX;
  fveinY{d} = fY;  
end
save('fveinfeatures.mat','fveinX','fveinY');
disp('Fingervein Features Extraction end and  saved on fveinfeatures.mat')

end


