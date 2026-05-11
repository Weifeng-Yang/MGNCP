%% Calculate the value of the objective function
function [loss,S]=computeCP(var,ngmar)
X=ktensor(var);
S=tensor(X);
loss=0.5*norm(ngmar-S)^2;
end

