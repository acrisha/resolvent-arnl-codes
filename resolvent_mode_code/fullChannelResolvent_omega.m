% last updated Dec 16 2025

%% Primitive-variable resolvent code for channel flow
% Written by Mitul Luhar 01/06/2015
% For details on the primitive-variable resolvent analysis see:
% Luhar, Sharma, McKeon (2014), Opposition Control within the Resolvent 
% Analysis Framework, J. Fluid Mech., 749:597-626
% modified for full channel and finite difference by E. Lenz on 5/9/2025

function [RIG,y,D1,D2,dy,U0,sqW,isqW] = fullChannelResolvent_omega(Re,kx,kz,omega,N,nsvd)
% Inputs

% N:    Wall-normal grid resolution
% nsvd: Number of singular modes to compute
% 
% % Flow Parameters
% Re:   Reynolds number (Re = u_tau h / nu)
% kx:   streamwise wavenumber (normalized by half-height,h)
% kz:   spanwise wavenumber (normalized by half-height,h)
% cP:   Wave speed in plus units (normalized by u_tau)
%omega = cP*kx;  % Radian frequency
cP=omega/kx;

% Set up coordinate system and estimate mean velocity profile
[y,D1,D2,dy,U0,dU0] = channelMeanVel_FD_180mfu();
% y:[0 2]
% D1: First differential
% D2: Second differential
% dy: Integration weights for scaling the resolvent
% U0: Mean velocity profile

% set up the linear NS operator and eventually the resolvent
% A few basic matrices
I   = eye(N);
Z   = zeros(N);

% Calculate mean shear
U0 = diag(U0);
dU0 = diag(dU0);

% Create important block components
ikU0 = 1i*kx*U0;
LAP  = -(kx^2)*I - (kz^2)*I + D2;
block1=.5*(-ikU0+LAP/Re+flip(flip(-ikU0+LAP/Re,1),2));

% Block matrix L representing linearized NS equations
% The last column is the pressure gradient, last row is continuity
L1 = [block1      , -dU0        , Z          , -1i*kx*I]; %u (ax.)
L2 = [Z           , block1      , Z          ,      -D1]; %v (rad.)
L3 = [Z           , Z           , block1     , -1i*kz*I]; %w (az.)
L4 = [-1i*kx*I    ,-D1          ,-1i*kz*I    ,        Z]; %contin.
L  = [L1; L2; L3; L4];

% Mass Matrix
M = [I Z Z Z; Z I Z Z; Z Z I Z; Z Z Z Z];

% The governing equation reads: (-i*om*M - L) [u;v;w;p] = M [fx;fz;fw;0]
LHS = -1i*omega*M-L;
RHS = M;

% Apply boundary conditions and compute resolvent
H0 = fullChannelBC(LHS,RHS,N);

% Scale resolvent and perform SVD
IW   = sqrtm(diag(dy));
iIW  = I/IW;
sqW  = [ IW Z Z Z; Z  IW Z Z; Z Z  IW Z; Z Z Z Z];
isqW = [iIW Z Z Z; Z iIW Z Z; Z Z iIW Z; Z Z Z Z];

% Weighted resolvents
H0W = sqW*H0*isqW;

% Singular value decomposition
% [u0W,s0W,v0W] = svds(H0W,nsvd,'largest','MaxIterations',1000);  % from emmas
[u0W,s0W,v0W] = svds(H0W,nsvd,'largest','MaxIterations',7000);  % lexi's because warning said maxiterations was not enought.

u0 = isqW*u0W; 
s0 = diag(s0W);
v0 = isqW*v0W;

% set phase of first non-zero point based on critical layer
%THIS IS NOW HAPPENING IN THE CHANNELSYM_VALIDATION SCRIPT
% if(cP < max(diag(U0)))
%     ind = find(diag(U0)>cP,1,'first');
% else
%     ind = N;
% end
% phase_shift = -1i*angle(u0(ind,:));
% v0 = v0*diag(exp(phase_shift));

% Because of the l2 norm used to scale the resolvent, we do not have any
% pressure data.  Calculate pressure modes using the un-scaled resolvent.

%u0_check = H0*v0*diag(1./s0); %matches u0
RIG.u = u0;
RIG.s = s0;
RIG.v = v0;

end
