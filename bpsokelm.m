
function [tbfs,gbest,gfs,gklm]=bpsokelm(tx,ty,vx,vy,np,T)

%T  number of iteration
%np  number of swarm
%n1 number of bits for alpha
%n2 number of bits for C

global n1 n2;
global train_x train_y;
global valid_x valid_y;

train_x = tx;train_y = ty;
valid_x = vx;valid_y = vy;

n1=8;
n2=8;
n3= size(train_x,1);

nb = n1+n2+n3;

pop=init(np,nb);
vel = 2*rand(np,nb) -1;
pbest = pop;
[fs,klms]=popfitness(pop);

lbfs = fs;
[gfs,mi]=max(fs);
gbest = pop(mi,:);
gklm = klms(mi);

wmax=0.9; wmin=0.1;
tbfs = zeros(T,1);

for t=1  : T
    
w = wmax- (t/T)*(wmax - wmin);
[pop,vel]=update(vel,pop,pbest,gbest,w); 
[fs,iklm]=popfitness(pop);   
    
I= (fs > lbfs);
lbfs(I) = fs(I); 
pbest(I) = pop(I);
klms(I) = iklm(I);

[mx,mi]=max(fs);
if(mx > gfs)
 gbest = pop(mi,:);
 gfs  = mx;
 gklm=klms(mi);
end

[gbest2,gfs2,gklm2]=localsearch(gbest,gfs);
if(gfs2 > gfs)
 pop(mi,:) = gbest2;   
 gbest = gbest2;
 gfs  = gfs2;
 gklm=gklm2;
 klms(mi) =gklm;
end

tbfs(t) = gfs;


end



function [gbest,gfs,gklm]=localsearch(gbest,gfs)

K=10; gklm =0;
nb = length(gbest); 
     
for n=1 : K
   
 bse=bsequence(nb); 
 Z=(bse & gbest); 
 [fs,klm] = fitness(Z); 
 if ( fs > gfs)
    gbest = Z;
    gfs = fs;    
    gklm=klm;
 end
    
end


function [fs,klms]=popfitness(x)

np=size(x,1); 
fs = zeros(np,1);
klm = struct('W',[],'bias',[],'B',[],'C',0.0); 
klms =repmat(klm,np,1);
for n=1 : np           
 [fs(n),klms(n)]=fitness(x(n,:));  
end



function [acc,klm]=fitness(xseq)

global train_x;global train_y;
global valid_x;global valid_y;

klm = struct('W',[],'bias',[],'B',[],'C',0.0); 
[alpha,C,fs]=splitsequence(xseq);

X = train_x(fs,:);
vX = valid_x(fs,:);

[nv,ns]=size(X); 
nh=floor(1.0*nv); % hidden neuron

klm.C= C;
klm.W = 2*rand(nh,nv)-1;
klm.bias = rand(nh,1);
T= klm.W*X + repmat(klm.bias,1,ns);
H=radbas(T);
nr=size(H,1);

ID = eye(nr)/C;
KE = ID+H*H';
klm.B= KE\H*train_y';


ns = size(vX,2);
T=klm.W*vX + repmat(klm.bias,1,ns);
H=radbas(T);
p=H'*klm.B;

Yp =(p>=0.5);
[~, Y] = max(valid_y);
[~, Yp] = max(Yp,[],2);
acc = sum(Y==Yp')/ numel(Y);

   

function [xn,vel]=update(vel,x,pbest,gbest,w) 

[np,nb]=size(x);
c1=0.2;c2=0.2;
vmax=5; vmin=-5;

for n=1 : np 
r1 = rand(1,nb);    
r2 = rand(1,nb);
vel(n,:) = w*vel(n,:)  +  c1*r1.*(pbest(n,:) - x(n,:)) + c2*r2.*(gbest - x(n,:)); 
end

I1=( vel(n,:)<vmin ); vel(n,I1)=vmin;
I2=( vel(n,:)>vmax ); vel(n,I2)=vmax;
S =  1./(1 + exp(-vel) );
xn =S>rand(np,nb) ;




function  pop=init(np,nb)

 pop = zeros(np,nb);
 for n=1 :np 
  pop(n,:)=bsequence(nb);       
 end
 

function se=bsequence(n)

se = zeros(1,n);
se(1) = rand;

u =4;
for n=2 : n
  se(n)  =  u * se(n-1)*( 1- se(n-1) ); 
end
se =se>mean(se);


