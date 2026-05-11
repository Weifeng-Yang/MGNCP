clearvars -except 
clc
warning('off');

%% Parameter.
%   index     : The dataset to be used, when index=1, use ConcreteCrackImages dataset
%               when index=2, use pixraw10P dataset
%   r         : Step factor
%   maxiteropt: Maximum iteration alloted to the method
%   trigger   : Whether to enable the indicator array of each method, where
%               when 1∈trigger, enable the MGNMF method
%               when 2∈trigger, enable the GSNMF method
%               when 3∈trigger, enable the 𝓁1-SNCP method
%               when 4∈trigger, enable the SGCP method
%               when 5∈trigger, enable the GSNTD method
%               when 6∈trigger, enable the 𝓁0-SGNCP method
%               when 7∈trigger, enable the MGNCP method
%   alphat    : The graph regularization parameters
%   stopindex : The indicator of the stop condition.  
%               To set the specific termination condition, see the 'stopcheck' function for details.  
%               The default termination condition is: ϵ<1e-5 or rum time>400 seconds 
%% Display
%   Rel       : The difference in the variable value between two iterations.


%% Parameter settings
rng('shuffle')
index=2;
r=1.01;
maxiteropt=8000;
trigger=[1,2,3,4,5,6,7];
% R=10;
alphat=[1,1,1];
stopindex=4;





%% Select dataset
[ngmar,R,Rdims,y]=readfile(index);
num=length(size(ngmar));
N=R;
X1=double(tenmat(ngmar,length(size(ngmar))))';







for j=1:10
%% Init
var=[];
for i=1:num
    var{i}=rand(size(ngmar,i),Rdims(num));
end


%% Solving
for i=1:length(trigger)
[datas{i},varss{i}]=ALGOchoose(var,ngmar,maxiteropt,Rdims,trigger(i),stopindex,r,alphat);
end
vart=var;
datas{length(trigger)+1}=varss;
datas{length(trigger)+2}=vart;
for i=1:length(trigger)
    vars1=datas{length(trigger)+1};
    vars11=vars1{i};
    vartemp=vars11{end};

    if(trigger(i)==1 || trigger(i)==2)
        [acc(j,i),rdx(j,i),NMIs(j,i)]=clustermeans(vartemp{end}',N,y);

    elseif(trigger(i)==5)

        cluster=ttm(tensor(vars11{1}),vartemp{end},num);
        temp=double(tenmat(cluster,num));
        [acc(j,i),rdx(j,i),NMIs(j,i)]=clustermeans(temp,N,y);

    else
         [acc(j,i),rdx(j,i),NMIs(j,i)]=clustermeans(vartemp{end},N,y);
    end

end
   datass{j}=datas;
end






%% Display results
accmean=mean(acc)
rdxmean=mean(rdx)
nmimean=mean(NMIs)
plt0=barplot(trigger,acc,rdx,NMIs); plt0=barplot(trigger,acc,rdx,NMIs);
