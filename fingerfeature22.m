function lbpc=fingerfeature22(dbase,dname)


%dname =1;
%dbase='mdbase\fingerprint\';
%dbase='encfinger\';
dpath= sprintf('%s%s\\*.tif',dbase,num2str(dname));
flist=dir(dpath);

fname ={};
for i=1 : length(flist)
fname{i} = flist(i).name;
end


[ptsc,sptc,Tim]=geopoints(dbase,dname,fname);
lbpc=lbpfeatures2(ptsc,sptc,Tim);   




function hbc=lbpfeatures2(ptsc,sptc,Timc)

nb=360;nfs = 10;
ns = length(ptsc);
hbc = zeros(ns,nb*nfs);
for i=1 : ns    
hbc(i,:)=lbpfeatures(ptsc{i},sptc{i},Timc{i});
end


function lbpc=lbpfeatures(Pts,spt,Tim)

nb=360;nfs = 10;
np=size(Pts,1);
lbpc=zeros(nb,nfs);

for j=1 : np     
  m=Pts(j,1);n=Pts(j,2);    
  img=Tim(m-15:m+16,n-15:n+16);    
  l1=extractLBPFeatures(img,'Upright',false,'Interpolation','Nearest','Normalization','L2');                               
  lbpc(spt(j),:) =l1;  
end

lbpc=reshape(lbpc,1,nb*nfs);


function [ptsc,sptc,eim]=geopoints(dbase,dname,fname)


Tim={};ptsc={};eim={};sptc={};
load('corepoints','cpcs');

L=length(fname);
fs=cell(1,L);
sflag =false(1,L);

for i=1: L

fpath=sprintf('%s%s\\%s',dbase,num2str(dname),fname{i});    

[pts,Tim{i},eim{i}]=morphprocess(fpath);
pts=remove_borderpoints(pts,size(eim{i}));

I=find( strcmp(cpcs{dname}(:,1), fname{i}));
crp=cpcs{dname}{I,2};

%dis = sqrt ( (pts(:,1)-crp(1)).^2 + ( pts(:,2)-crp(2)).^2 );
%if( size(pts,1) > 30 )
 %[~,mi]=sort(dis);
 %pts=pts(mi(1:30),:);
%end

if( ~isempty(crp) )    

distance = sqrt ( (pts(:,1)-crp(1)).^2 + ( pts(:,2)-crp(2)).^2 );
y = ( pts(:,1) - crp(1) );  x=( pts(:,2) - crp(2) );

angle= atan2(x,y);
angle =(angle + pi ).*(360/(2*pi));
fs{i} = [distance angle]; 
ptsc{i} =pts;
sflag(i) = true;    
    
end


end


fs = fs(sflag);
ptsc = ptsc(sflag);
Tim = Tim(sflag);
eim = eim(sflag);



nb=1:size(fs{1},1);  
sptc{1} =ceil(fs{1}(:,2)');
%[nb; ceil(fs{1}(:,2)')]
L = length(fs);

for i=2 : L

 nb=1:size(fs{i},1);  
 ds=pdist2(fs{1},fs{i} );
 [mv,mi]=min(ds,[],1);
 
 I=(mv <= 20);
 nb=nb(I);
 mi = mi(I);
 
 %[nb; mi; ceil(fs{1}(mi,2)')] 
 %sptc{i} =[nb; mi; ceil(fs{1}(mi,2)')];
 
 ptsc{i}=ptsc{i}(nb,:); 
 sptc{i} =ceil(fs{1}(mi,2)');
  
end



function  Pts=remove_borderpoints(Pts,sz)

np=size(Pts,1);
S=[];M=sz(1);N=sz(2);
for j=1 : np     
    m=Pts(j,1);n=Pts(j,2);  
    m1 = m-15; m2 = m+16;
    n1 = n-15; n2 = n+16;
    if ( m1>0 && m2<=M && n1>0 && n2<=N) 
        S=[S;j];
    end
    
end
Pts=Pts(S,:);
Pts=remove_closestpoints(Pts,32);


function  Pts=remove_closestpoints(Pts,th)
np=size(Pts,1);
cpf=false(np,1);
for j=1 : np         
    if ( cpf(j)==false )  
        ds=pdist2(Pts(j,:),Pts);      
        [mv,mi]=sort(ds);    
        mv=mv<th;mv(1)=false;
        cpf(mi(mv))=true;
    end    
end
Pts=Pts(~cpf,:);






function [Pts,Tim,eim]=morphprocess(fpath)

eim=imread(fpath);
se = strel('disk',1);
eim=imclose(eim,se);

% minitia points
Tim = bwmorph(eim,'thin',Inf);

% ridge ends points
Ep = bwmorph(Tim,'endpoints');
[Ie,Je]=find(Ep);
[Ie2,Je2]=edgeindex(Tim,Ie,Je);

% ridge branch points
Bp = bwmorph(Tim,'branchpoints');
[Ib,Jb]=find(Bp);
[Ib2,Jb2]=edgeindex(Tim,Ib,Jb);

% combine ridge end and branch points 
mI = [Ie2; Ib2];
mJ = [Je2; Jb2];

Pts = [mI mJ];


function [Ie,Je]=edgeindex(Tim,Ie,Je)

[M,N]=size(Tim);
Fm =[];m=1;

for k=1 : M
   
 kj=find(Tim(k,:),1,'first');   
 if ( ~isempty(kj) )    
   Fm(m,:) = [k kj];   
   m=m+1;
 end  
 kj=find(Tim(k,:),1,'last');  
 if ( ~isempty(kj) )       
     Fm(m,:) = [k kj];
     m= m+1;
 end    
end


for k=1 : N
   
 kj=find(Tim(:,k),1,'first');   
 if ( ~isempty(kj) )    
   Fm(m,:) = [kj k];   
   m=m+1;
 end  
 kj=find(Tim(:,k),1,'last');  
 if ( ~isempty(kj) )       
     Fm(m,:) = [kj k];
     m= m+1;
 end    
end

S = false(length(Ie),1);
for n=1 : length(Ie)
 S(n)= any(Fm(:,1)==Ie(n)  & Fm(:,2)==Je(n) );    
end

%[length(Ie) length(Je) length(S) sum(S)]
Ie = Ie(~S);Je = Je(~S);
