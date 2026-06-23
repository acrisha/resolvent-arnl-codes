%% plotting of scatter plots
clc, clear, close all;

% This scripts reads in savekeepvalsXXX.mat
% keeper_vals = size(XXXXXXXXXX x 10)
% kx1, kz1, kx2, kz2, isval1, isval2, isval3, alpha, kx, kz

% savekxplotXXXX.mat
% kx_plot = size(XX x XX)

% 



%% user inputs
% which files / data do you want to load in?
filename1 = "savekxplotRe180_Nx6_Nz4_svd1_omegacase0_filteredmean.mat";
filename2 = "savekeepvalsRe180_Nx6_Nz4_svd1_omegacase0_filteredmean.mat" ;
filename3 = "channel180_moser_loaddata.mat";

% used for wave number calculations
gpx = 6;  % grid points in streamwise (based on filenameX input values)
gpz = 4;   %                spanwise 
delta = 1; % half channel height


% what do we want to plot?
% kx' vx kx'' (case 1)
% kx vx kz    (case 2)
plotchoice = 2;

% figure settings
manual_bounds = true;   % if we want to make the plot look pretty and hard code the bounds
savefigure = false;      % if we want to save figure
FS = 20;  % fontsize
LW = 1.5; % linewidth



%% load and organize data
load(filename1)
load(filename2)
params = load(filename3);

LxDNS = params.Lx;
LzDNS = params.Lz;
val = 9;
valskip = 2;

kx1_ind = keeper_vals(:,1);
kz1_ind = keeper_vals(:,2);
kx2_ind = keeper_vals(:,3);
kz2_ind = keeper_vals(:,4);
eigval1 = keeper_vals(:,5);
eigval2 = keeper_vals(:,6);
eigval3 = keeper_vals(:,7);
alpha = keeper_vals(:,8);
kx_ind = keeper_vals(:,9);
kz_ind = keeper_vals(:,10);

%% discrete selection of rows
% this is used for selecting chosing what to plot (i.e., all svd = 1 etc)
% rows2keep = (eigval1 == 1) & (eigval2 == 1) & (eigval3 == 1);
% keeper_rows = [kx1_ind(rows2keep,:) kz1_ind(rows2keep,:) kx2_ind(rows2keep,:) kz2_ind(rows2keep,:) kx_ind(rows2keep,:) kz_ind(rows2keep,:) alpha(rows2keep,:)];


%% prepare plotting variables
determine keeper rows for plotting based on what we want to plot
switch plotchoice
    case 1 % kx' vx kx'' (case 1)
        disp('Case 1: Plotting kx'' vx kx''''')
        keeper_rows = [kx1_ind kx2_ind alpha];
    case 2 % kx vx kz    (case 2)
        disp('Case 2: Plotting kx vx kz')
        keeper_rows = [kx_ind kz_ind alpha];
    otherwise
        disp('Error: not plotting anything ...')
        keyboard
end

disp('sorting pairs, determining max alpha, making result matrix ....')

% find the unique pairs of wave number combinations
pairs = [keeper_rows(:,1), keeper_rows(:,2)];  %kz1 kz2
[unique_pairs, ~, idx] = unique(pairs, 'rows');

% determine the maximum alpha of the repeated pairs
alpha_keeper = keeper_rows(:,3);
max_alpha_of_pairs = accumarray(idx, alpha_keeper, [], @max);

% create an array with the unique pair and alpha combinations
result = [unique_pairs, max_alpha_of_pairs];
alpha = result(:,3)/max(result(:,3),[],'all'); % normalized alpha we plot


% determine plot parameters for each case
kprime_index = result(:,1);
kdoubprime_index = result(:,2);


switch plotchoice
    case 1 % kx' vx kx'' (case 1)
        xlabelstring = {'$K_1''$'};
        ylabelstring = {'$k_1''''$'};
        cmap = slanCM('OrRd');
        coefficienta = 50;  % for visualization on scatter plot
        coefficientb = 800;
        kprime = kprime_index * 2*pi * delta / LxDNS;
        kdoubprime = kdoubprime_index * 2*pi * delta / LxDNS;


        figure
        pcolor(-gpx:gpx,-gpx:gpx,kx_plot/max(kx_plot,[],'all')); colorbar
        c = colorbar; c.TickLabelInterpreter = 'latex'; clim([0 1]);
        colormap(cmap);
        ylabel(c, '$\frac{\alpha}{\alpha_{max}}$','Interpreter','latex','Rotation',0,'FontSize',24);
        ax = gca; ax.FontSize = FS; ax.TickDir = 'in'; ax.LineWidth = LW; ax.TickLabelInterpreter = 'Latex';
        set(gca,'YDir','normal');
        axis square
        xlabel('$-n_1'':n_1''$','FontSize',FS,'interpreter','latex')
        ylabel('$-n_1'''':n_1''''$','FontSize',FS,'interpreter','latex')
        title('Quick Check: integers: $n_1 = n_1'' + n_1''''$')

    case 2 % kx vx kz    (case 2)
        xlabelstring = {'$k_3$'};
        ylabelstring = {'$k_1$'};
        cmap = slanCM('GnBu');
        coefficienta = 50;
        coefficientb = 1200;
        kprime = kprime_index * 2*pi * delta / LxDNS;
        kdoubprime = kdoubprime_index * 2*pi * delta / LzDNS;

    otherwise
        disp('Error: not plotting anything ...')
        keyboard
end

%% plotting
hfig = figure;
s1fig = subplot(1,2,1);
sizes = coefficienta * (alpha) + 1; % plus one is just because sizes cannot have 0s. Just for visual purposes
scatter(kdoubprime,kprime, sizes, alpha, 's', 'filled'); % Scatter plot with square markers
c = colorbar; c.TickLabelInterpreter = 'latex'; clim([0 1]);
colormap(cmap);
ylabel(c, '$\frac{\alpha}{\alpha_{max}}$','Interpreter','latex','Rotation',0,'FontSize',24);
ax = gca; ax.FontSize = FS; ax.TickDir = 'in'; ax.LineWidth = LW; ax.TickLabelInterpreter = 'Latex';
set(gca,'YDir','normal');
axis square
xlabel(xlabelstring,'FontSize',FS,'interpreter','latex')
ylabel(ylabelstring,'FontSize',FS,'interpreter','latex')
grid off;
xtickangle(0)
ax.XAxis.Label.Rotation = 0;
ax.YAxis.Label.Rotation = 0;

s2fig = subplot(1,2,2);
sizes = coefficientb * (alpha) + 1;
scatter(kdoubprime,kprime, sizes, alpha, 's', 'filled','MarkerEdgeColor','black'); % Scatter plot with square markers
c = colorbar; c.TickLabelInterpreter = 'latex'; clim([0 1]);
colormap(cmap);
ylabel(c, '$\frac{\alpha}{\alpha_{max}}$','Interpreter','latex','Rotation',0,'FontSize',24);
ax = gca; ax.FontSize = FS; ax.TickDir = 'in'; ax.LineWidth = LW; ax.TickLabelInterpreter = 'Latex';
set(gca,'YDir','normal');
axis square
xlabel(xlabelstring,'FontSize',FS,'interpreter','latex')
ylabel(ylabelstring,'FontSize',FS,'interpreter','latex')
grid off;
xtickangle(0)
ax.XAxis.Label.Rotation = 0;
ax.YAxis.Label.Rotation = 0;



%% manually hard code bounds for visualization puropses
if manual_bounds 

% hard coded bounds
switch plotchoice
    case 1 % kx' vx kx'' (case 1)
        hfig;
        s1fig;
        xlim([-64 64]); ylim([-64 64]);
        xticks(-64:16:64); yticks(-64:16:64);
        text(s1fig, -0.23, 0.99, '(a)', 'Units', 'normalized', 'FontSize', FS, 'Interpreter','latex', 'FontWeight','bold');

        s2fig;
        yline(-val:2:val,'LineWidth',0.75)
        xline(-val:2:val,'LineWidth',0.75)
        xlim([-9 9]); ylim([-9 9]);
        xticks(-8:2:8); yticks(-8:2:8);
        text(s2fig, -0.23, 0.99, '(b)', 'Units', 'normalized', 'FontSize', FS, 'Interpreter','latex', 'FontWeight','bold');


    case 2 % kx vx kz    (case 2)
       hfig;
       s1fig;
        xlim([-128 128]);
        ylim([-64 64]);
        xticks(-128:32:128); yticks(-64:16:64);
        text(s1fig, -0.23, 0.99, '(a)', 'Units', 'normalized', 'FontSize', FS, 'Interpreter','latex', 'FontWeight','bold');

        s2fig;
        yline(-9:2:9,'LineWidth',0.75)
        xline(-18:4:18,'LineWidth',0.75)
        xlim([-18 18]);
        ylim([-9 9]);
        xticks(-16:4:16);
        yticks(-8:2:8);
        text(s2fig, -0.23, 0.99, '(b)', 'Units', 'normalized', 'FontSize', FS, 'Interpreter','latex', 'FontWeight','bold');


    otherwise
        disp('Error: not plotting anything ...')
        keyboard
end


end

%% set figure size
set(hfig, 'Units', 'inches', 'Position', [1 1 16 6]);

%% save figure
if savefigure
    % Remove extra white space around subplots
    set(hfig, 'PaperUnits', 'inches');
    set(hfig, 'PaperPositionMode', 'auto');
    set(hfig, 'PaperSize', [16 6]);  % Match the figure size exactly

    % OPTIONAL: Make subplots tighter (if needed)
    % tight_layout = true;
    % if tight_layout
    %     % Works in newer MATLAB versions
    %     tiledlayout(2,2,'Padding','compact','TileSpacing','compact');
    % end

    % Save as PDF with no additional white space
    print(hfig, 'tester.pdf', '-dpdf', '-painters');
end



