clc, clear; close all;
% Compute resolvent modes for a span of kx, kz, for a critical omega
% x = streamwise, y = wall-normal, z = spanwise


%% user inputs
% what Reynolds number channel we want to look at?
Retau_target = 180; % 180, 395

% what omega case do we want to run
% 0 for omega = 0
% 1 for omega = critical value
omegacase = 0; 


% use for calculating critical frequency (if omegacase = 1)
yplus_crit = 15;          % innerscaled wall normal grid point of interest 
lambdaxplus_crit = 1000;  % inner scaled streamwise wavelength  aqssociated with wall normal location 
% ^^ (can use premultiplied energy spectra or literature to determine these)


checkdivergence = false; % do you want to check divergence and output in command window? % normally used for debugging



%% load in full channel parameters 
load(['channel',num2str(Retau_target),'_moser_loaddata.mat']); nx=Nx; nz=Nz; ny = Ny;  

% ensure Umean is saved at cell edges!!
Umean = 0.5*(Umean + Umean(end:-1:1))/utau;  %symmetrize about centerline

% save original profile parameters (these will be reloaded in each time)
yg_original = yg;
y_original = y;
Umean_original = Umean;


% use for actual mode calculations
Nx = nx/2;   % only retain half of streamwise grid points
Nz = nz/2;   % only retain half of spanwise grid points
Ny = ny;     % total number of edge points in wall-normal direction, including walls   % should be same at RNLGO p.z_w (129)
Nkeep = 16;   % max amount of modes we can keep, Nkeep is number of modes to keep, prefered DIVISIBLE BY 4


% % this is used for debugging!!!!! (comment out for running full mode calc)
% Nx = 4;
% Nz = 2;
% Ny = ny;
% Nkeep = 2;


%% print to command window the cases we are running
fprintf('Retau = %i, omegacase = %i \n',Retau_target,omegacase)
if omegacase == 0
    disp('Calculating omega = 0')
else
    disp('Calculating critical omega')
end
fprintf('Nx = %i, Nz = %i, svd = %i \n\n',Nx,Nz,Nkeep)




%% initiialize storge matrices
sval = zeros(Nx+1,Nz+1, Nkeep);     
uvec = zeros(Nx+1,Nz+1,4*Ny, Nkeep);
fvec = zeros(Nx+1,Nz+1,4*Ny, Nkeep);
uvec_sym = zeros(Nx+1,Nz+1,4*Ny, Nkeep);
fvec_sym = zeros(Nx+1,Nz+1,4*Ny, Nkeep);



%% determine omegas
if omegacase == 0
    omega = 0;

elseif omegacase == 1

    % OLD VERSION (hardcoded values)
    % %find the critical omega for y+=15, lambda+=1000
    % ygplus=yg*retau; %nearest grid point is y+=14.5, which is entry 28
    % critwavespeed=Umean(28);
    % Lxplus=Lx*retau; %length of box in plus units
    % kxcrit=(Lxplus/1000)*2*pi/Lx; %should be about 1
    % omegacrit=critwavespeed*kxcrit; %resolvent code is set up so that wave speed is in plus units but kx is scaled by Lx ... need to figure out why


    %find the critical omega (based on parameter input aboove)
    ygplus = yg*retau;                                      % inner scale wall normal grid
    ygindex_crit = find(yplus_crit>=ygplus, 1, 'last' );    % matlab index for target wall normal location
    U_crit = Umean_original(ygindex_crit);               % critical wavespeed

    Lxplus=Lx*retau; %length of box in plus units
    kx_crit=(Lxplus/lambdaxplus_crit)*2*pi/Lx; %should be about 1

    omega_crit=U_crit*kx_crit; %resolvent code is set up so that wave speed is in plus units but kx is scaled by Lx ... need to figure out why

    omega = omega_crit;
    fprintf('critical omega = %f\n',omega)
end




%% compute resolvent modes
for kkx = 0:Nx
      disp(['kx=',num2str(kkx)]);

  for kkz = 0:Nz

    if (kkx + kkz ~= 0) % skip the (0,0) mode
      % disp(['kx=',num2str(kkx),', kz=',num2str(kkz)]);

      kx = kkx*2*pi/Lx;
      kz = kkz*2*pi/Lz;

      [RIG,y,D1,~,~,~,sqW,~] = fullChannelResolvent_omega(retau,kx,kz,omega,Ny,Nkeep,Umean_original,utau,yg_original,y_original);

      sval(kkx+1,kkz+1,:) = RIG.s;       % stored as if it's been through fftshift()
      uvec(kkx+1,kkz+1,:,:) = RIG.u;
      fvec(kkx+1,kkz+1,:,:) = RIG.v;

      % symmetrize
      for my = 1:Nkeep/2
          [uvec_sym_part, fvec_sym_part] = channelSym_validation(squeeze(uvec(kkx+1,kkz+1,:,(2*my-1):(2*my))), squeeze(fvec(kkx+1,kkz+1,:,(2*my-1):(2*my))), Ny,kx,kz,D1,y,sqW);
          uvec_sym(kkx+1,kkz+1,:,(2*my-1):(2*my))=uvec_sym_part;  % psi
          fvec_sym(kkx+1,kkz+1,:,(2*my-1):(2*my))=fvec_sym_part;  % phi
          %compare to unsymmetrized modes
          uvec_part=squeeze(uvec(kkx+1,kkz+1,:,(2*my-1):(2*my)));
          fvec_part=squeeze(fvec(kkx+1,kkz+1,:,(2*my-1):(2*my)));
          if 0
              close all; figure(); subplot(1,4,1) 
              plot(real(uvec_sym_part(1:129,1)),y,'b'); hold on; plot(real(uvec_sym_part(1:129,2)),y,'r');
              plot(real(uvec_part(1:129,1)),y,'b:'); hold on; plot(real(uvec_part(1:129,2)),y,'r:');
              title('Response u');
              subplot(1,4,2);
              plot(real(uvec_sym_part(259:387,1)),y,'b'); hold on; plot(real(uvec_sym_part(259:387,2)),y,'r');
              plot(real(uvec_part(259:387,1)),y,'b:'); plot(real(uvec_part(259:387,2)),y,'r:');
              title('Response w');
              subplot(1,4,3) 
              plot(real(fvec_sym_part(1:129,1)),y,'b'); hold on; plot(real(fvec_sym_part(1:129,2)),y,'r');
              plot(real(fvec_part(1:129,1)),y,'b:'); hold on; plot(real(fvec_part(1:129,2)),y,'r:');
              title('Forcing u');
              subplot(1,4,4);
              plot(real(fvec_sym_part(259:387,1)),y,'b'); hold on; plot(real(fvec_sym_part(259:387,2)),y,'r');
              plot(real(fvec_part(259:387,1)),y,'b:'); hold on; plot(real(fvec_part(259:387,2)),y,'r:');
              title('Forcing w');
              keyboard;
          end

          %uvec_part_prev=uvec_part; fvec_part_prev=fvec_part; uvec_sym_part_prev=uvec_sym_part; fvec_sym_part_prev=fvec_sym_part;
      end %symmetrization

      if checkdivergence
          for nind=1:Nkeep
              %check divergence
              Nycenter=ceil(Ny/2); Nycentervec=[Nycenter Nycenter+1];
              divorig_uvec = 1i*kx*squeeze(uvec(kkx+1,kkz+1,1:Ny,nind)) + D1*squeeze(uvec(kkx+1,kkz+1,(Ny+1):(2*Ny),nind)) + 1i*kz*squeeze(uvec(kkx+1,kkz+1,(2*Ny+1):(3*Ny),nind));
              divorig_fvec = 1i*kx*squeeze(fvec(kkx+1,kkz+1,1:Ny,nind)) + D1*squeeze(fvec(kkx+1,kkz+1,(Ny+1):(2*Ny),nind)) + 1i*kz*squeeze(fvec(kkx+1,kkz+1,(2*Ny+1):(3*Ny),nind));
              divsym_uvec = 1i*kx*squeeze(uvec_sym(kkx+1,kkz+1,1:Ny,nind)) + D1*squeeze(uvec_sym(kkx+1,kkz+1,(Ny+1):(2*Ny),nind)) + 1i*kz*squeeze(uvec_sym(kkx+1,kkz+1,(2*Ny+1):(3*Ny),nind));
              divsym_fvec = 1i*kx*squeeze(fvec_sym(kkx+1,kkz+1,1:Ny,nind)) + D1*squeeze(fvec_sym(kkx+1,kkz+1,(Ny+1):(2*Ny),nind)) + 1i*kz*squeeze(fvec_sym(kkx+1,kkz+1,(2*Ny+1):(3*Ny),nind));
              %skip centerpoints of symmetric modes
              divorig_uvec(Nycentervec)=[]; divorig_fvec(Nycentervec)=[]; divsym_uvec(Nycentervec)=[]; divsym_fvec(Nycentervec)=[];
              disp(['kx=',num2str(kkx),', kz=',num2str(kkz),', nind=',num2str(nind)]);
              disp(['orig uvec: ',num2str(norm(divorig_uvec)),', orig fvec: ',num2str(norm(divorig_fvec)),', uvec: ',num2str(norm(divsym_uvec)),', fvec: ',num2str(norm(divsym_fvec))])
              if norm(divorig_uvec)>1e-11 || norm(divorig_fvec)>2e-9 || norm(divsym_uvec)>1e-11 || norm(divsym_fvec)>2e-9
                  keyboard;
              end
          end
      end
    end %if statement for skipping (0,0) mode

  end %kkz
  % disp(['finishing kkx=',num2str(kkx)]);  
end %kkx


uvec=uvec_sym; fvec=fvec_sym;

%% save resolvent modes
if omegacase == 0
    filename = ['modes',num2str(Retau_target),'_moser_june2026_omegazero.mat'];
    save('-v7',filename,'sval','uvec','fvec','omega','y','Nx','Ny','Nz','D1','sqW','Nkeep');


elseif omegacase == 1
    filename = ['modes',num2str(Retau_target),'_moser_june2026_omegacrit',num2str(omega),'.mat'];
    save('-v7',filename,'sval','uvec','fvec','omega','y','Nx','Ny','Nz','D1','sqW','Nkeep','yplus_crit','lambdaxplus_crit');
end


