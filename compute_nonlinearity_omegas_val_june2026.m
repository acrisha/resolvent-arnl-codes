clear, close all;
% last updated June 22, 2026
% Run with validationflag=1 to check that we are calculating the nonlinear term correctly
% Run with omegacase=0,1,or 2 for different combinations of omega=0 and omega=omega_critical
% note: u, v, and w are all stored at cell edges
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% user inputs
% what Reynolds number channel we want to look at?
Retau = 180;   % 180, 395, or (eventually) 590 (from Moser database)

% do we want to validate symmetry and projection are computed correctly?
validationflag = 0;    %0 to skip validation, 1 to run validation


% what omega case do we want to run: (where omega = omega' + omega'')
% 0: all three omega = 0                 
% 1: for first response mode omega=omega_critical, second response mode omega=0 (forcing mode omega=omega_critical)                    
% 2: for first response mode omega=0, second response mode omega=omega_critical (forcing mode omega=omega_critical)
omegacase = 0; 



%% load in channel parameters
load(['channel',num2str(Retau),'_moser_loaddata.mat'],'Lx','Lz','Nx','Nz','Ny','y','ym');
% full channel number of gridpoints
nx = Nx; nz = Nz; ny = Ny; clearvars Nx Ny Nz; % clear for consitency (avoid overwriting)

Nx = nx/2; % streamwise
Nz = nz/2; % spanwise
Ny = ny;   % wall normal
nsvd = 2;  % max amount of modes we can keep, Nkeep is number of modes to keep, preferred divisibility by 2
Nkeep = nsvd;

% % used for debugging
% Nx = 6; % streamwise
% Nz = 2; % spanwise
% Ny = ny;   % wall normal
% nsvd =  2;  % max amount of modes we can keep, Nkeep is number of modes to keep, preferred divisibility by 2
% Nkeep = nsvd;


%% load in modes from resolvent calculator code
if omegacase == 0
    load(['modes',num2str(Retau),'_moser_june2026_omegazero.mat'],'fvec','uvec','sval','D1','sqW');
    uvec1=uvec; fvec1=fvec; sval1=sval;
    uvec2=uvec; fvec2=fvec; sval2=sval;
    uvec3=uvec; fvec3=fvec; sval3=sval;
    clearvars uvec fvec sval

elseif omegacase == 1
    % critical omega
    load(['modes',num2str(Retau),'_moser_june2026_omegacrit.mat'],'fvec','uvec','sval','D1','sqW');
    uvec1=uvec; fvec1=fvec; sval1=sval; %critical omega
    uvec3=uvec; fvec3=fvec; sval3=sval; %critical omega

    % omega = 0
    load(['modes',num2str(Retau),'_moser_june2026_omegazero.mat'],'fvec','uvec','sval','D1','sqW');
    uvec2=uvec; fvec2=fvec; sval2=sval; %omega=0
    clearvars uvec fvec sval

elseif omegacase == 2
    % critical omega
    load(['modes',num2str(Retau),'_moser_june2026_omegacrit.mat'],'fvec','uvec','sval','D1','sqW');
    uvec2=uvec; fvec2=fvec; sval2=sval; %critical omega
    uvec3=uvec; fvec3=fvec; sval3=sval; %critical omega
    
    % omega = 0
    load(['modes',num2str(Retau),'_moser_june2026_omegazero.mat'],'fvec','uvec','sval','D1','sqW');
    uvec1=uvec; fvec1=fvec; sval1=sval; %omega=0
    clearvars uvec fvec sval
end



%% initialization of variables
filesavename1 = ['saveNLprojRe',num2str(Retau),'_Nx',num2str(Nx),'_Nz',num2str(Nz),'_svd',num2str(Nkeep),'_omegacase',num2str(omegacase),'.mat'];
filesavename2 = ['savekeepvalsRe',num2str(Retau),'_Nx',num2str(Nx),'_Nz',num2str(Nz),'_svd',num2str(Nkeep),'_omegacase',num2str(omegacase),'.mat'];
filesavename3 = ['savekxplotRe',num2str(Retau),'_Nx',num2str(Nx),'_Nz',num2str(Nz),'_svd',num2str(Nkeep),'_omegacase',num2str(omegacase),'.mat'];

sqW=sqW(1:Ny*3,1:Ny*3); %remove the pressure part of the mode
proj = zeros(2*Nx+1,2*Nz+1,2*Nx+1,Nz+1,nsvd,nsvd,nsvd);
Errvec=zeros(2*Nx+1,2*Nz+1,2*Nx+1,Nz+1,Nkeep); %this will be zero if not running validation
counter = 0; % only used if running with validationflag





%% computing nonlinearity with triadic relationships
% loop over all wavenumbers kx1, kz1, kx2, kz2

disp('starting NL calculation')
tic
for kkx1 = -Nx:Nx
    disp(kkx1)

    for kkz1 = -Nz:Nz
        for kkx2 = -Nx:Nx
            for kkz2 = 0:Nz
                % disp(['(kx1,kz1) = (', num2str(kkx1),',',num2str(kkz1),') , (kx2,kz2) = (', num2str(kkx2),',',num2str(kkz2),')'])


                % check for any nonzero wave numbers (exluding kz',kz'' = 0)
                % if kkx1 ~= 0 && kkx2 ~= 0 && kkz1 ~= 0 && kkz2 ~= 0 && (kkx1 + kkx2 ~= 0 || kkz1 + kkz2 ~= 0)
                % if abs(kkx1) ~= 0 && abs(kkx2) ~= 0 && (kkx1 + kkx2 ~= 0 && kkz1 + kkz2 ~= 0)
                if  (kkx1 + kkx2 ~= 0 && kkz1 + kkz2 ~= 0 && kkx1 ~= 0 && kkx2 ~= 0)


                    % check to make sure wave number combinations are within domain bounds.
                    if (abs(kkx1 + kkx2) <= Nx) && (abs(kkz1 + kkz2) <= Nz)

                        % get actual wavenumbers
                        kx1 = kkx1*2*pi/Lx;
                        kz1 = kkz1*2*pi/Lz;

                        kx2 = kkx2*2*pi/Lx;
                        kz2 = kkz2*2*pi/Lz;

                        % load modes
                        for isval1 = 1:nsvd
                            [u1, v1, w1] = nmode(uvec1,kkx1,kkz1,isval1,Ny); %column vectors
                            sigma1=sval1(abs(kkx1)+1,abs(kkz1)+1,isval1);

                            for isval2 = 1:nsvd
                                [u2, v2, w2] = nmode(uvec2,kkx2,kkz2,isval2,Ny);
                                sigma2=sval2(abs(kkx2)+1,abs(kkz2)+1,isval2);


                                fu = (1i*(kx1+kx2)* (u1.*u2)) + (D1 * (u1.*v2)) + (1i*(kz1+kz2) * (u1.*w2));
                                fv = (1i*(kx1+kx2)* (v1.*u2)) + (D1 * (v1.*v2)) + (1i*(kz1+kz2) * (v1.*w2));
                                fw = (1i*(kx1+kx2)* (w1.*u2)) + (D1 * (w1.*v2)) + (1i*(kz1+kz2) * (w1.*w2));



                                innerprod=zeros(nsvd,1);
                                for isval3=1:nsvd
                                    [u3, v3, w3] = nmode(fvec3,kkx1+kkx2,kkz1+kkz2,isval3,Ny);
                                    innerprod(isval3)=[u3; v3; w3]'*(sqW.^2)*[fu; fv; fw];
                                end


                                %VALIDATION
                                if validationflag
                                    kkx1f=-kkx1-kkx2; kkz1f=-kkz1-kkz2;
                                    kkx2f=kkx2; kkz2f=kkz2;
                                    isval3f=isval1;
                                    [u3f, v3f, w3f] = nmode(uvec1,kkx1f+kkx2f,kkz1f+kkz2f,isval3f,Ny);
                                    kx1f = kkx1f*2*pi/Lx;
                                    kz1f = kkz1f*2*pi/Lz;

                                    innerprodf=zeros(size(innerprod));

                                    if sum(abs(u3f-conj(u1)))
                                        keyboard
                                    end


                                    for isval1f=1:nsvd
                                        [u1f, v1f, w1f] = nmode(fvec3,kkx1f,kkz1f,isval1f,Ny);
                                        sigma1f=sval3(abs(kkx1f)+1,abs(kkz1f)+1,isval1f);

                                        %                                     fuf = u1f.* (1i*kx2 * u2 ) + v1f.* (D1 * u2 ) + w1f.* (1i*kz2 * u2)+...
                                        %                                           u2 .* (1i*kx1f* u1f) + v2 .* (D1 * u1f) + w2 .* (1i*kz1f* u1f);
                                        %                                     fvf = u1f.* (1i*kx2 * v2 ) + v1f.* (D1 * v2 ) + w1f.* (1i*kz2 * v2)+...
                                        %                                           u2 .* (1i*kx1f* v1f) + v2 .* (D1 * v1f) + w2 .* (1i*kz1f* v1f);
                                        %                                     fwf = u1f.* (1i*kx2 * w2 ) + v1f.* (D1 * w2 ) + w1f.* (1i*kz2 * w2)+...
                                        %                                           u2 .* (1i*kx1f* w1f) + v2 .* (D1 * w1f) + w2 .* (1i*kz1f* w1f);

                                        %                                     fuf = u2 .* (1i*kx1f* u1f) + v2 .* (D1 * u1f) + w2 .* (1i*kz1f * u1f);
                                        %                                     fvf = u2 .* (1i*kx1f* v1f) + v2 .* (D1 * v1f) + w2 .* (1i*kz1f * v1f);
                                        %                                     fwf = u2 .* (1i*kx1f* w1f) + v2 .* (D1 * w1f) + w2 .* (1i*kz1f * w1f);

                                        fuf = (1i*(kx1f+kx2)* (u1f.*u2)) + (D1 * (u1f.*v2)) + (1i*(kz1f+kz2) * (u1f.*w2));
                                        fvf = (1i*(kx1f+kx2)* (v1f.*u2)) + (D1 * (v1f.*v2)) + (1i*(kz1f+kz2) * (v1f.*w2));
                                        fwf = (1i*(kx1f+kx2)* (w1f.*u2)) + (D1 * (w1f.*v2)) + (1i*(kz1f+kz2) * (w1f.*w2));

                                        innerprodf(isval1f)=[u3f; v3f; w3f]'*(sqW.^2)*[fuf; fvf; fwf];


                                        [u3r,~,~] = nmode(fvec3,kkx1+kkx2,kkz1+kkz2,isval1f,Ny);
                                        if sum(abs(u1f-conj(u3r)))
                                            keyboard
                                        end



                                    end

                                    if ( sum( (abs(innerprod+innerprodf)>2e-6) )~=(Nkeep/2) )&& (kkz2 && kkz1) %makes sure symmetric modes are cancelling out as expected
                                        keyboard
                                    end

                                    %what determines whether they need to be
                                    %added or subtracted?
                                    diffinnerprod=zeros(size(innerprod));
                                    for j=1:length(diffinnerprod)
                                        if abs(innerprod(j)+innerprodf(j)) < abs(innerprod(j)-innerprodf(j))
                                            diffinnerprod(j)=abs(innerprod(j)+innerprodf(j)) * (abs(innerprod(j)+innerprodf(j))>2e-6); %if it's less than 10^-6, set it to machine zero. This gets around large errors at very small magnitudes of innerprod
                                        else
                                            diffinnerprod(j)=abs(innerprod(j)-innerprodf(j)) * (abs(innerprod(j)-innerprodf(j))>2e-6);
                                        end


                                    end
                                    % if norm(abs(innerprod+innerprodf)) < norm(abs(innerprod-innerprodf))
                                    %     diffinnerprod=abs(innerprod+innerprodf) .* ( abs(innerprod+innerprodf)>1e-6 ); %if below 1e-10, set to zero
                                    % else
                                    %     diffinnerprod=abs(innerprod-innerprodf) .* ( abs(innerprod-innerprodf)>1e-6 ); %if below 1e-10, set to zero
                                    % end
                                    Errvec(kkx1+Nx+1,kkz1+Nz+1,kkx2+Nx+1,kkz2+1,:)=diffinnerprod./(.5*(abs(innerprod)+abs(innerprodf)));

                                    if any(Errvec(kkx1+Nx+1,kkz1+Nz+1,kkx2+Nx+1,kkz2+1,:) >0.5)
                                        %keyboard
                                        counter=counter+.25;
                                    end
                                end %end validation if-statement

                                % save projection
                                proj(kkx1+Nx+1,kkz1+Nz+1,kkx2+Nx+1,kkz2+1,isval1,isval2,:) = ...
                                    innerprod*abs(sigma1*sigma2);

                                % stored as: proj(kx',kz',kx'',kz'',my',my'',my)


                            end
                        end
                    end
                end
            end

        end
    end
end
toc

disp('finished  NL calculation')
save('-v7.3',filesavename1,'proj','Errvec','counter','omegacase');
disp('projection saved.')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plotting and post processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('sorting alpha.')
[alpha,sort_index] = sort(abs(proj(:)),1,'descend');

keeper_vals = zeros(length(alpha),8);

for i=1:length(alpha)
    [kx1, kz1, kx2, kz2, isval1, isval2, isval3] = ind2sub(size(proj), sort_index(i));
    keeper_vals(i,:) = [kx1-Nx-1, kz1-Nz-1, kx2-Nx-1, kz2-1, isval1, isval2, isval3, alpha(i)];
end

kx = keeper_vals(:,1) + keeper_vals(:,3);
kz = keeper_vals(:,2) + keeper_vals(:,4);

keeper_vals = [keeper_vals, kx, kz];
save('-v7.3',filesavename2,'keeper_vals')
keeper_valsorig = keeper_vals;



%% pulls top triadic interaction for each wave number combination despite the kz value
disp('starting to bin for plotting')
count1 = 0;

for ii = -Nx:Nx
    count2 = 0;
    count1 = count1 + 1;

    for jj = -Nx:Nx
        count2 = count2 + 1;
        ind = keeper_vals(:,1) == ii;
        keep_temp = keeper_vals(ind,:);
        clear ind
        ind = keep_temp(:,3) == jj;
        final_temp = keep_temp(ind,:);
        kx_plot(count1,count2) = final_temp(1,8);

    end
end
save('-v7.3',filesavename3,'kx_plot')



%%
x_plot = -Nx:Nx;
y_plot = -Nx:Nx;
kxdeg = 2*pi*x_plot / Lx;


figure;
imagesc(x_plot, y_plot, kx_plot/max(kx_plot,[],'all'));
axis xy;  % Flip Y-axis to match matrix layout
ax = gca; ax.FontSize = 18; ax.TickDir = 'in'; ax.LineWidth = 1.5; ax.TickLabelInterpreter = 'Latex';
c = colorbar; 
% clim([0 1.5]); 
c.TickLabelInterpreter = 'Latex'; %c.Ticks = ([ 0 0.25 0.5 0.75 1 1.25 1.5 1.75]);
c.Label.String = '$\alpha$'; c.Label.Interpreter = 'Latex'; c.Label.Rotation = 360;
cmap = cbrewer2('seq', 'Blues', 100); % cmap = flipud(cmap);
colormap(cmap);

xlabel('$k_x''$','fontsize',20,'Interpreter','Latex')
ylabel('$k_x''''$','fontsize',20,'Interpreter','Latex')

% xticks(x_plot);
% yticks(x_plot);
% xticklabels(kxdeg)
% yticklabels(kxdeg)

% title('Nonlinear interaction coefficient where $k_x = k_x'' + k_x''''$','fontsize',16,'Interpreter','Latex')

