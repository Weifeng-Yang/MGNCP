%%  All parameters of this function are explained the same as 'main_Run_me' and 'ALGOchoose' functions
function [var,loss,timerun]=MGNMF(var,ngmar,maxiteropt,stopindex,r,btmax)
%% initialization algorithm
loss=[];
rate=0;
timerun=[0];
ngmar=double(tenmat(ngmar,length(size(ngmar))))';
num=length(size(ngmar));
R=size(var{num},2);
rho=1e-8;
var=[];

for i=1:num
    var{i}=rand(size(ngmar,i),R);
end
var{num}=var{num}';

Lapk=MulLLaplace(ngmar);

LK=zeros(1,num);
L=ones(1,num);
tk=1;
bts=[];
wk=zeros(1,num);
varK=var;


returnloss=norm(ngmar,"fro")^2;

for j=1:length(Lapk)
    alpha(j)=1/length(Lapk);
end


loss(1)=computeloss(ngmar,var,Lapk,alpha);






t1=clock;


for i=1:maxiteropt
%% update parameters
fprintf("%d\n",i);
vv=var;

%% Update parameters
for j=1:num
    wk(j)=min(wk(j),btmax);
    vv{j}=var{j}+wk(j)*(var{j}-varK{j});
    varK{j}=var{j};    
    LK(j)=L(j);
    [V,L(j)]=gradMGNMF(var,ngmar,j,r,Lapk,alpha);
    var{j}=PROXn1(V);
    vv{j}=var{j};
end
loss(i+1)=computeloss(ngmar,var,Lapk,alpha);

sumnorm=0;
for j=1:length(var)
    sumnorm=sumnorm+norm(var{j}-varK{j},'fro')^2;
end

%% Judging whether to extrapolate
if(loss(i+1)>loss(i)-rho/2*sumnorm)
    var=varK;
    for j=1:num
    [V,L(j)]=gradMGNMF(var,ngmar,j,r,Lapk,alpha);
    var{j}=PROXn1(V);
    end
    loss(i+1)=computeloss(ngmar,var,Lapk,alpha);
end


bts{i}=wk;
t2=clock;
timerun(i+1)=etime(t2,t1);
fprintf("MGNMF\n");
check1=0;
check2=0;
%% Check if termination condition is met
for j=1:num
%     fprintf("nonzero:%d\n",nnz(var{j}~=0));
    check1=check1+norm(var{j}-varK{j},'fro');
    check2=check2+norm(varK{j},'fro');
end
Res=check1/check2;
% Res=abs(loss(i+1)-loss(i))/returnloss;
fprintf("cri：%d\n",Res);
stop=stopcheck(Res,timerun,stopindex);
if(stop==1)
    fprintf("Number of terminations：%d\n",i);
    pause(4);
    break;
end


tk=(1+sqrt(1+4*tk^2))/2;
for j=1:num
wk(j)=(tk-1)/(tk);
end


end

end


function x=PROXn1(x)
x=double(x);
x(x<0)=0;
end

function ck=checkck(ck)
    if(ck==0) 
        a=0;
       while(a==0)
           a=rand(1);
       end
       ck=a;
     end
end

function [U,L]=gradMGNMF(var,ngmar,n,r,Lapk,alpha)
    if(n==1)
    ZT=var{n+1}; 
    ck=norm(ZT*ZT','fro');
    ck=checkck(ck);
    mar=var{n}*ZT;
    U=var{n}-1/(r*ck)*(mar-ngmar)*ZT'; 
    L=ck;
    else
    Z=var{1};
    ck=norm(Z'*Z,'fro');
    Lapktemp=0;
        for j=1:length(Lapk)
        Lapktemp=alpha(j)*var{end}*Lapk{j}+Lapktemp;
        ck=ck+alpha(j)*norm(Lapk{j},'fro');
        end
    mar=Z*var{n};
    L=ck;
    U=var{n}-1/(r*L)*(Z'*(mar-ngmar)+Lapktemp); 
    end
    
end



function loss=computeloss(ngmar,var,Lapk,alpha)
    loss=compute(var,length(var),ngmar);
    for i=1:length(Lapk)
    loss=loss+alpha(i)*trace(var{end}*Lapk{i}*var{end}');
    end
end

function loss=compute(var,num,ngmar)
    nga=var{1};
    for i=2:num
        nga=nga*var{i};
    end
    loss=norm(ngmar-nga,'fro');
end
