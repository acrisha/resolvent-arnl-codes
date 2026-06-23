%% Plots DNS from Moser database text files
%% Converts half channel data to full channel dataset (with mirroring)
%% Saves input .mat files for resolvent codes
clc, clear, close all;
% last modified A.C. Risha, JHU, June 2026

% READ ME:
% This file:
% 1. Outlines the parameters for the Moser DNS channel data sets for Retau
%    180, 395, and 590
% 2. Reads in chanXX_means.txt files which are half channel data sets from https://turbulence.oden.utexas.edu/data/MKM/ text files
% 3. Mirrors the half channel data sets into full data sets and saves full
%    channel data as fullchannellXXmoser_june2026.mat files (if flag is turned
%    on)
% 4. Converts the Moser full channel data (chebychev WN grid) to a
%    hyperbolic tangent stretched WN grid used in LESGO simulations and saves
%    data for running resolvent codes. 

% Inputs needed in main directory to run this script:
%    chanXXX_means.txt

% Outputs from running this script (if savedata = true):
%    fullchannelXXmoser_june2026.mat
 

%% user inputs
% Do you want to save output data file?
savedatafile = true;

% Colors for plotting data (180, 395, 590 cases)
Colorvals = {'red','cyan','blue'};

%% full channel parameters from Moser, Kim, & Mansour 1999:
targetRetau_vec = [180, 395, 590]; %targetRetau_vec=targetRetau_vec(2);
actualRetau_vec = [178.13, 392.24, 587.19]; %actualRetau_vec=actualRetau_vec(2);
Nx_vec = [128, 256, 384]; % streamwise grid points 
Ny_vec = [129, 257, 257]; % wall-normal grid points
Nz_vec = [128, 192, 384]; % spanwise grid points
%Nx_vec=Nx_vec(2); Ny_vec=Ny_vec(2); Nz_vec=Nz_vec(2);
% Note: Ny grid is based on https://turbulence.oden.utexas.edu/data/MKM/ text files (not table in paper)
delta = 1;   % half channel height
utau = 1;    % friction velocity
Lx_vec = [4*pi*delta, 2*pi*delta, 2*pi*delta]; % streamwise domain
Ly_vec = [2, 2, 2];                            % wall-normal domain
Lz_vec = [4/3*pi*delta, pi*delta, pi*delta];   % spanwise domain
%Lx_vec=Lx_vec(2); Ly_vec=Ly_vec(2); Lz_vec=Lz_vec(2);


%% 
for i = 1:length(targetRetau_vec)

    retau = targetRetau_vec(i);
    actualretau = actualRetau_vec(i);
    filename = ['chan',num2str(retau),'_means.txt'];
    % y/delta, y^+,  U, dU/dy,  W, P
    dns = importdata(filename);

    ydns = dns.data(:,1);
    yplusdns = dns.data(:,2);
    Udns = dns.data(:,3);

    % plot half channel data set
    figure(1)
    title('plotting raw data from Moser database')
    semilogx(yplusdns,Udns,'color',Colorvals{i},'linestyle','none','linewidth',1.5,'marker','o')
    hold on
    xlabel('y^+')
    ylabel('U^+')
    axis square



    % calculate the full channel dataset
    full_Udns = [Udns; flipud(Udns(1:end-1))];
    full_ydns = [ydns; 2*ydns(end) - flipud(ydns(1:end-1))];
    full_yplusdns = [yplusdns; 2*yplusdns(end) - flipud(yplusdns(1:end-1))];

    filesavename = ['fullchannel',num2str(retau),'moser_june2026.mat'];
    Nx = Nx_vec(i);
    Ny = Ny_vec(i);
    Nz = Nz_vec(i);

    Lx = Lx_vec(i);
    Ly = Ly_vec(i);
    Lz = Lz_vec(i);

    if savedatafile
        save(filesavename,'full_Udns','full_ydns','full_yplusdns','actualretau','retau','utau','delta','Lx','Ly','Lz','Nx','Ny','Nz');
    end


    % plot full channel dataset (raw data from Moser)
    figure(2)
    subplot(1,2,1)
    title('plotting full channel data (raw data mirrored) ')
    plot(full_ydns,full_Udns,'color',Colorvals{i},'linestyle','none','linewidth',1.5,'marker','o')
    hold on
    xlabel('y')
    ylabel('U^+')
    axis square
    
    subplot(1,2,2)
    title('plotting full channel data (raw data mirrored) ')
    semilogx(full_yplusdns,full_Udns,'color',Colorvals{i},'linestyle','none','linewidth',1.5,'marker','o')
    hold on
    xlabel('y^+')
    ylabel('U^+')
    axis square



end
figure(1); legend('Retau180','Retau395','Retau590','location','northwest','fontsize',14)
figure(2);subplot(1,2,1);legend('Retau180','Retau395','Retau590','location','northwest','fontsize',14)



