function mmbioauthclassifier


clc; clear all;

disp('MultiModal Biometric  Recognition');

if ( exist('fingerprintfeatures.mat','file')~=2 )
     disp('fingerprintfeatures.mat features file not found');
     return;
end
load('fingerprintfeatures.mat','fingerX','fingerY');



if ( exist('fveinfeatures.mat','file')~=2 )
     disp('fveinfeatures.mat features file not found');
     return;
end
load('fveinfeatures.mat','fveinX','fveinY');


if ( exist('facefeatures.mat','file')~=2 )
     disp('facefeature.mat file not found');
     return;
end

load('facefeatures.mat','faceX','faceY');



fingerX=removeZeroColumns(fingerX);
fveinX=removeZeroColumns(fveinX);
faceX=removeZeroColumns(faceX);

[X,Y]=combinefeatures(fingerX,fveinX,faceX);

X= X';

disp( 'Features Set Size' );
disp('row - number of samples');
disp('col - number of features');

[ns,nf] = size(X);
disp('Features size');
disp([ns,nf])


[nv,ns]=size(X);
uc = unique(Y);
nc = length(uc);

Y2=zeros(nc,ns);
for i=1 : ns    
 Y2(Y(i),i)=1;     
end

% normalization -   min-max
for n=1 : nf
 mx=max(X(:,n));  
 mn=min(X(:,n));    
 X(:,n)= ( (X(:,n) - mn)./ (mx-mn)  )  ;
end


% normalization -zsocre
%me= mean(X);st = std(X);
%for i=1 : ns
%X(i,:) = (X(i,:)-me)./st  ; 
%end


trainRatio=0.7;valRatio =0.15;
testRatio=0.15;
[trainInd,valInd,testInd] = dividerand(ns,trainRatio,valRatio,testRatio);


disp('Spliting Dataset for ')
disp(['Training :', num2str(trainRatio)] )
disp(['validation :',num2str(valRatio)] )
disp(['Testing :',num2str(testRatio)] )

% training  set 
train_x = X(:,trainInd);  train_y=Y2(:,trainInd);
%validation  set
valid_x = X(:,valInd); valid_y = Y2(:,valInd);
% testing  set 
test_x = X(:,testInd); test_y = Y2(:,testInd);

disp( 'Training Set Size' );
disp('number of variables  -  number of samples');
disp( size(train_x) )

disp( 'Validation Set Size' );
disp('number of variables  -  number of samples');
disp( size(valid_x) )

disp( 'Test Set Size' );
disp('number of variables  -  number of samples');
disp(size(test_x))


if ( exist('psofeatures.mat','file')~=2 )
        
% optimization modified binary pso
np =20; % number of swarm  increase it up to 80
T =5; % number of iteration increase it up to 100

[tbfs,gbest,gfs,gklm] = bpsokelm(train_x,train_y,valid_x,valid_y,np,T);
[alpha,C,sf]=splitsequence(gbest);

figure;plot(1:T,tbfs,'-*');
xlabel('Iteration') ; ylabel('fitness value');
title('Modified Binary PSO Optimization Performance');

save('psofeatures.mat', 'sf','gklm');

end

load('psofeatures.mat', 'sf','gklm');
cp1=performance(sf,train_x,train_y,gklm,'Training');
cp2=performance(sf,valid_x,valid_y,gklm,'Validation');
cp3=performance(sf,test_x,test_y,gklm,'Testing');
save('mbaperformance.mat', 'cp1','cp2','cp3');


function cp=performance(fs,X,Y,gklm,stitle)

X=X(fs,:);
[~, Yc] = max(Y);
Yp=mbaclassify(X,gklm);
Yp=Yp';

fprintf('\n %s -: performance \n',stitle);
cp = classperf(Yc,Yp)


function Yp=mbaclassify(X,gklm)

ns = size(X,2);
T=gklm.W*X + repmat(gklm.bias,1,ns);
H=radbas(T);
p=H'*gklm.B;
Yp =(p>=0.5);
[~, Yp] = max(Yp,[],2);



function [X,Y]=combinefeatures(fX1,fX2,fX3)

nuser = length(fX1);
X = [];Y = [];

for n= 1: nuser
    
 X1 = fX1{n};X2 = fX2{n};
 X3 = fX3{n};
 
 nx1=size(X1,1);nx2=size(X2,1);
 nx3=size(X3,1);
 
 mx = max([nx1,nx2,nx3]);
 rx1=randperm(nx1,nx1);
 if ( nx1< mx )
  rx=ceil( nx1*rand(1,mx-nx1));
  rx1=[rx1 rx];
 end
 rx2=randperm(nx2,nx2);
 if ( nx2< mx )
  rx= ceil(nx2*rand(1,mx-nx2));
  rx2=[rx2 rx];   
 end
 
 rx3=randperm(nx3,nx3);
 if ( nx3< mx )     
  rx=ceil(nx3*rand(1,mx-nx3));
  rx3=[rx3 rx];
 end   
 
 
 cX1= horzcat(X1(rx1,:),X2(rx2,:),X3(rx3,:));  
 X = [X; cX1];  
 Y = [Y; ones(mx,1)*n];
 
end



function fX=removeZeroColumns(fX)

nu = length(fX);
X = [];

for n= 1: nu
 X = [X; fX{n}];   
end

%check zero and nan colunms available
nS = sum(isnan(X));
S=sum(X);
I=(S==0) |  nS>0;

for n= 1: nu
 fX{n} = fX{n}(:,~I);    
end




 
