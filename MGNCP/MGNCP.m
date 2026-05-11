%%  All parameters of this function are explained the same as 'main_Run_me' and 'ALGOchoose' functions
function [var,loss,timerun,bts]=MGNCP(var,ngmar,maxiteropt,stopindex,r,alphat)
%% initialization algorithm
loss=[];
timerun=[0];
num=length(size(ngmar));

LK=zeros(1,num);
L=ones(1,num);
tk=1;
Lapk=MulLLaplace(ngmar);
bts=[];
wk=zeros(1,num);
varK=var;
for j=1:num
wk(j)=(tk-1)/(tk);
end

% for i=1:length(Lapk)
%     alpha(i)=1/length(Lapk);
% end
alpha=alphat;

for j=1:num
    fprintf("nonzero:%d\n",nnz(var{j}~=0));
end


loss(1)=computeloss(ngmar,var,Lapk,alpha);


t1=clock;


for i=1:maxiteropt
%% update parameters
fprintf("%d\n",i);
vv=var;

for j=1:num
    wk(j)=min(wk(j),0.9999*sqrt(LK(j)/L(j)));
    %% Update parameters
    vv{j}=var{j}+wk(j)*(var{j}-varK{j});
    varK{j}=var{j};
    LK(j)=L(j);
    [V,L(j)]=gradMGNCP(vv,ngmar,j,num,r,Lapk,alpha);
    var{j}=PROXn1(V);
end
loss(i+1)=computeloss(ngmar,var,Lapk,alpha);


%% Judging whether to extrapolate
if(loss(i+1)>loss(i))
    var=varK;
    for j=1:num
    [V,L(j)]=gradMGNCP(var,ngmar,j,num,r,Lapk,alpha);
    var{j}=PROXn1(V);
    end
    loss(i+1)=computeloss(ngmar,var,Lapk,alpha);
end



%% Check if termination condition is met
fprintf("MGNCP\n");
check1=0;
check2=0;
for j=1:num
%     fprintf("nonzero:%d\n",nnz(var{j}));
    check1=check1+norm(var{j}-varK{j},'fro');
    check2=check2+norm(varK{j},'fro');
end



bts{i}=wk;
t2=clock;
timerun(i+1)=etime(t2,t1);
Res=check1/check2;
fprintf("Rel：%d\n",Res);
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


function loss=computeloss(ngmar,var,Lapk,alpha)
    loss=computeCP(var,ngmar);
    num=length(var);
    for i=1:length(Lapk)
    loss=loss+alpha(i)/2*trace(var{num}'*Lapk{i}*var{num});
    end
end

function [U,L]=gradMGNCP(var,ngmar,n,num,r,Lapk,alpha)
  [Xtemp,temp]=krob2(var,n,num,ngmar);
   if(n==num)
        ck=norm(temp,'fro');
        Lapktemp=0;
        for j=1:length(Lapk)
        Lapktemp=alpha(j)*Lapk{j}*var{num}+Lapktemp;
        ck=ck+alpha(j)*norm(Lapk{j},'fro');
        end
        mar=var{n}*temp-Xtemp+Lapktemp;
        L=ck;
   else
        ck=norm(temp,'fro');
        L=ck;
        mar=var{n}*temp-Xtemp;
   end

   tao=1/(r*ck);
   U=var{n}-tao*mar;
end

function x=PROXn1(x)
x=double(x);
x(x<0)=0;
end


