%% Interpolate chebvychev to hyperbolic tangent (full channel data):
%% Saves input .mat files for resolvent codes
clc, clear, close all;
% last modified EL, June 2026

% READ ME:
% This file:
% 1. Reads in full channel data (from ConvertPlotMoserDNS_june.m)
% 2. Converts the Moser full channel data (chebychev WN grid) to a
%    hyperbolic tangent stretched WN grid used in LESGO simulations and saves
%    data for running resolvent codes. 

% Inputs needed in main directory to run this script:
%    fullchannelXXmoser_june2026.mat

% Outputs from running this script (if savedata = true):
%    channelXXX_moser_loaddata.mat
%    ^^^ this is what we will load in the resolvent codes 
 
clc, clear, close all;

%% user inputs
% Do you want to save output data file?
savedatafile = true;


% Colors for plotting data (180, 395, 590 cases)
Colorvals = {'red','cyan','blue'};



%%
targetRetau_vec = [180, 395, 590]; 
for i = 1:length(targetRetau_vec)
    retau = targetRetau_vec(i);

    filename = ['fullchannel',num2str(retau),'moser_june2026.mat'];
    moser = load(filename);
    ymoser = moser.full_ydns; %edges
    umoser = moser.full_Udns; %originally at edges
    actualretau=moser.actualretau;


    % hyperbolic tangent interpolation
    WNpts = moser.Ny; % number of wall normal points.
    h = moser.Ly; % channel height
    a = 5.2; % stretch factor (used in ARNL simulations)
    dy0 = h/WNpts; %original length of each segment

    y1 = linspace(0,h,WNpts)'; % w-values of original grid (cell edges)
    y1m = y1+0.5*dy0;          % uv-values of original grid (cell centers)
    y_rnl = h*0.5*(1+tanh(a*(y1/h-0.5))/tanh(0.5*a)); % y5 is the stretched grid
    Umean_rnl = spline(ymoser, umoser, y_rnl);

    figure(1)
    plot(ymoser, umoser,'color',Colorvals{i},'linestyle','none','linewidth',1.5,'marker','o')
    hold on
    plot(y_rnl, Umean_rnl,'color','black','linestyle','-','linewidth',1.5,'marker','.','markersize',8)
    legend('Original', 'Interpolated');
    xlabel('y'); ylabel('U');



    % saving for .mat file
    %Umean = zeros(length(Umean_rnl)+1,1);
    %Umean(1) = 0;
    %Umean(2:end) = Umean_rnl;
    %Umean(end) = 0;

    Umean=Umean_rnl; y=y_rnl; %edges
    ym=.5*(y(2:end)+y(1:end-1)); %centers
    yg=[-ym(1); ym; 2+ym(1)]; %centers, including through-wall ghosts

    %yg = zeros(length(y_rnl)+1,1);
    %ghost = y_rnl(1);
    %yg(1) = ghost;
    %yg(2:end) = y_rnl;
    %ym = y_rnl(1:end-1);
    %y = y_rnl;    % this is not entirely correct. this is stored on faces (which we do not have that information)


    utau = moser.utau;
    Lx = moser.Lx;
    Ly = moser.Ly;
    Lz = moser.Lz;
    Nx = moser.Nx;
    Ny = moser.Ny;
    Nz = moser.Nz;


if savedatafile
    filesavename = ['channel',num2str(retau),'_moser_loaddata.mat'];
    save(filesavename, 'retau','utau','Umean','yg','ym','y','Lx','Ly','Lz','Nx','Ny','Nz','actualretau');
end

end












