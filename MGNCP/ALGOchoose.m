
% Input.
% var,core    : initial matrix and core tensor
% ngmar       : decomposed tensor
% The remaining parameters are explained the same as the 'main_Run_me' function

% Output.
% cores, vars : Decomposition matrix and core tensor resulting from the final iterative result
% loss:       : Array of loss functions generated during iteration
% tr:         : Runtime array during iteration
%             : where 'atz' represents the array of additional extrapolated parameters for IBPG and iPALM


function [data,varss]=ALGOchoose(var,ngmar,maxiteropt,Rdims,flag,stopindex,r,alphat)

if(flag==1)
[vars,loss,tr]=MGNMF(var,ngmar,maxiteropt,stopindex,r,1);
varss{1}=vars;
lossdata=loss;
trdata=tr;


elseif(flag==2)
[vars,loss,tr]=GSNMF(var,ngmar,maxiteropt,stopindex,r,1);
varss{1}=vars;
lossdata=loss;
trdata=tr;



elseif(flag==3)
[vars,loss,tr]=L1SNCP(var,ngmar,maxiteropt,stopindex,r,1e-4,2);
varss{1}=vars;
lossdata=loss;
trdata=tr;



elseif(flag==4)    
[vars,loss,tr]=SGCP(var,ngmar,maxiteropt,stopindex,r,2);
varss{1}=vars;
lossdata=loss;
trdata=tr;



elseif(flag==5) 
[cores,vars,loss,tr]=GSNTD(var,ngmar,maxiteropt,Rdims,stopindex,r,1.7,10);
varss{1}=cores;
varss{2}=vars;
lossdata=loss;
trdata=tr;


elseif(flag==6) 
[vars,loss,tr]=L0SGNCP(var,ngmar,maxiteropt,stopindex,r,1,0.99);  
varss{1}=vars;
lossdata=loss;
trdata=tr;


elseif(flag==7) 
[vars,loss,tr]=MGNCP(var,ngmar,maxiteropt,stopindex,r,alphat); 
varss{1}=vars;
lossdata=loss;
trdata=tr;
   

end


data{1}=lossdata;
data{2}=trdata;



end