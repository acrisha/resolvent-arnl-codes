%% Function to impose boundary conditions and compute the resolvent for half channel
function H0 = fullChannelBC(LHS,RHS,N)
% The NS equations have been expressed as:
% (-1i*om*M-L)[u;v;w;p] = M[fx;fy;fz;0] -> LHS [u;v;w;p] = RHS [fx;fy;fz;0]

% This function imposes the BCs corresponding to no slip and 
% computes the resolvent

%Inputs:
% LHS: block matrix (-1i*omega*M-L)
% RHS: block matrix M
% N: grid resolution

%Outputs:
% H0: resolvent for no-slip

% Impose no slip
LHS0 = LHS;
RHS0 = RHS;
for ni = [1,N,N+1,2*N,2*N+1,3*N]
    LHS0(ni,:)=0;
    RHS0(ni,:)=0;
    LHS0(ni,ni)=1;
end
H0 = LHS0\RHS0;

end

