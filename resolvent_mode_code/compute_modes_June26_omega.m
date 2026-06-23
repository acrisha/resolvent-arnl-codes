clear; close all;
% Compute resolvent modes for a span of kx, kz, for a critical omega
% x = streamwise, y = wall-normal, z = spanwise

% load in minimal full channel parameters 
load('channel180_minchan_new.mat'); 
Umean = 0.5*(Umean + Umean(end:-1:1))/utau;  %symmetrize about centerline


Nx = nx/2;   % only retain half of streamwise grid points
Nz = nz/2;   % only retain half of spanwise grid points
Ny = ny;     % total number of edge points in wall-normal direction, including walls   % should be same at RNLGO p.z_w (129)




Nkeep = 16;   % max amount of modes we can keep, Nkeep is number of modes to keep, prefered DIVISIBLE BY 4
% this is used for debugging:
% Nx = 4;
% Nz = 4;
% Ny = ny;

% initiialize storge matrices
sval = zeros(Nx+1,Nz+1, Nkeep);     
uvec = zeros(Nx+1,Nz+1,4*Ny, Nkeep);
fvec = zeros(Nx+1,Nz+1,4*Ny, Nkeep);
uvec_sym = zeros(Nx+1,Nz+1,4*Ny, Nkeep);
fvec_sym = zeros(Nx+1,Nz+1,4*Ny, Nkeep);
 

%find the critical omega for y+=15, lambda+=1000
ygplus=yg*Retau; %nearest grid point is y+=14.5, which is entry 28
critwavespeed=Umean(28);
Lxplus=Lx*Retau; %length of box in plus units
kxcrit=(Lxplus/1000)*2*pi/Lx; %should be about 1

omegacrit=critwavespeed*kxcrit; %resolvent code is set up so that wave speed is in plus units but kx is scaled by Lx ... need to figure out why

% compute resolvent modes
for kkx = 0:Nx
  for kkz = 0:Nz

    if (kkx + kkz ~= 0) % skip the (0,0) mode

        
      kx = kkx*2*pi/Lx;
      kz = kkz*2*pi/Lz;

      %cP = 0;           % wave speed in plus units - this was needed in the original function, but we modified it to take omega as an input directly
      N = Ny;           % number of wall normal discretization points
      
      [RIG,y,D1,~,~,~,sqW,~] = fullChannelResolvent_omega(Retau,kx,kz,omegacrit,N,Nkeep);


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


        % for nind=1:Nkeep
        %     %check divergence
        %     Nycenter=ceil(Ny/2); Nycentervec=[Nycenter Nycenter+1];
        %     divorig_uvec = 1i*kx*squeeze(uvec(kkx+1,kkz+1,1:Ny,nind)) + D1*squeeze(uvec(kkx+1,kkz+1,(Ny+1):(2*Ny),nind)) + 1i*kz*squeeze(uvec(kkx+1,kkz+1,(2*Ny+1):(3*Ny),nind));
        %     divorig_fvec = 1i*kx*squeeze(fvec(kkx+1,kkz+1,1:Ny,nind)) + D1*squeeze(fvec(kkx+1,kkz+1,(Ny+1):(2*Ny),nind)) + 1i*kz*squeeze(fvec(kkx+1,kkz+1,(2*Ny+1):(3*Ny),nind));
        %     divsym_uvec = 1i*kx*squeeze(uvec_sym(kkx+1,kkz+1,1:Ny,nind)) + D1*squeeze(uvec_sym(kkx+1,kkz+1,(Ny+1):(2*Ny),nind)) + 1i*kz*squeeze(uvec_sym(kkx+1,kkz+1,(2*Ny+1):(3*Ny),nind));
        %     divsym_fvec = 1i*kx*squeeze(fvec_sym(kkx+1,kkz+1,1:Ny,nind)) + D1*squeeze(fvec_sym(kkx+1,kkz+1,(Ny+1):(2*Ny),nind)) + 1i*kz*squeeze(fvec_sym(kkx+1,kkz+1,(2*Ny+1):(3*Ny),nind));
        %     %skip centerpoints of symmetric modes
        %     divorig_uvec(Nycentervec)=[]; divorig_fvec(Nycentervec)=[]; divsym_uvec(Nycentervec)=[]; divsym_fvec(Nycentervec)=[];
        %     disp(['kx=',num2str(kkx),', kz=',num2str(kkz),', nind=',num2str(nind)]);
        %     disp(['orig uvec: ',num2str(norm(divorig_uvec)),', orig fvec: ',num2str(norm(divorig_fvec)),', uvec: ',num2str(norm(divsym_uvec)),', fvec: ',num2str(norm(divsym_fvec))])
        %     if norm(divorig_uvec)>1e-11 || norm(divorig_fvec)>2e-9 || norm(divsym_uvec)>1e-11 || norm(divsym_fvec)>2e-9
        %         keyboard;
        %     end
        % end

    end %if statement for skipping (0,0) mode

  end %kkz
  disp(['finishing kkx=',num2str(kkx)]);  
end %kkx


uvec=uvec_sym; fvec=fvec_sym;
save('-v7',['modes_180_minchan_June26_omegacrit.mat'],'sval','uvec','fvec','omegacrit','y','Nx','Ny','Nz','D1','sqW','Nkeep');

