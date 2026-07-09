%%  Code to get the turbulent mean velocity profile, modified for finite difference in wall-normal direction
% last updated May 29



% Outputs
% U0: mean velocity profile in plus-units (i.e. scaled by u_tau) at cell edges
% y: wall-normal edge points
% dy: integration weights for wall normal coordinate
% D1: First differential matrix
% D2: Second differential matrix
% dU0: first derivative of mean velocity profile

%%
function [y,D1,D2,dy,U0, dU0] = channelMeanVel_FD_180mfu()



% yg = p.z_uv + an extra cell -4.e-6 at cell 1. (should be 1x130)
% ym = p.z_uv(1:end-1) (1x128) remove the 2.0004 because out of domain. 
% y = p.z_w (1 x129)



%{

load('channel180_minchan_new.mat','Retau','Umean','utau','yg','y');
Umean = 0.5*(Umean + Umean(end:-1:1))/utau;  %symmetrize about centerline
%I don't know how or when these got rounded off, but make sure they're
%perfectly symmetrized
y=.5*(y+(2-y(end:-1:1))); yg=.5*(y(2:end)+y(1:end-1)); yg=[-yg(1); yg; 2+yg(1)];

% y=y'; yg=yg'; %make sure everything is column vectors % emma version
y=y(:); yg=yg(:); %make sure everything is column vectors  % lexi version


dy=yg(2:end)-yg(1:end-1); %weight of each edge point for SVD

U0= interp1(yg,Umean,y,'spline'); %mean velocity at cell edges (easier for BCs this way)
dU0= (Umean(2:end)-Umean(1:end-1))./(yg(2:end)-yg(1:end-1));  %value at cell edges, slightly more accurate than D1*U0


%}

load('channel180_moser_loaddata.mat','Umean','utau','yg','y');
Umean = 0.5*(Umean + Umean(end:-1:1))/utau;  %symmetrize about centerline
%I don't know how or when these got rounded off, but make sure they're
%perfectly symmetrized
y=.5*(y+(2-y(end:-1:1))); yg=.5*(y(2:end)+y(1:end-1)); yg=[-yg(1); yg; 2+yg(1)];

Umean = Umean';
y=y'; yg=yg'; 


dy=yg(2:end)-yg(1:end-1); %weight of each edge point for SVD

%U0= interp1(yg,Umean,y,'spline'); %mean velocity at cell edges (easier for BCs this way)
U0=Umean; %already at cell edges


%first derivative - central difference, except at wall and centerline
D1vec=1./(y(3:end)-y(1:end-2));
D1vec=[1/(y(2)-y(1)) D1vec 1/(y(end)-y(end-1))];
D1=-diag(D1vec(2:end),-1)+diag(D1vec(1:end-1),1);
D1(1,1)=-D1vec(1); D1(end,end)=D1vec(end);

dU0=(D1*(U0'))';
%dU0= (Umean(2:end)-Umean(1:end-1))./(yg(2:end)-yg(1:end-1));  %value at cell edges, slightly more accurate than D1*U0


%second derivative - central difference, except at walls
D2diagvec=-.5./(y(3:end)-y(1:end-2)).*( 1./(y(3:end)-y(2:end-1)) + 1./(y(2:end-1)-y(1:end-2)) ); 
D2diagvec=[0 D2diagvec 0]; %will add boundaries later, just getting the right dimensions for now
D2overvec=.5./((y(3:end)-y(1:end-2)).*(y(3:end)-y(2:end-1)));
D2overvec=[0 D2overvec];
D2undervec=.5./((y(3:end)-y(1:end-2)).*(y(2:end-1)-y(1:end-2)));
D2undervec=[D2undervec 0];
%fix roundoff issues again
D2overvec=.5*(D2overvec+D2undervec(end:-1:1)); D2undervec=D2overvec(end:-1:1);

D2=diag(D2diagvec)+diag(D2overvec,1)+diag(D2undervec,-1);

%forward/backward difference at walls (this will get overwritten by BC)
D2(1,:)=D2(2,:); D2(end,:)=D2(end-1,:);

end
