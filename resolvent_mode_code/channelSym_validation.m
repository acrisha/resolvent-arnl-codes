function [psi_sym, phi_sym] = channelSym_validation(psi, phi,Ny,kx,kz,D1,y,sqW)

% Enforces symmetry of resolvent forcing and response modes about channel
% centerline (y = 0).

cntr = floor( (Ny+1)/2 ); % index closest to channel center

psi_sym = psi;
phi_sym = phi;

[maxuval, ind] = max( abs(psi(1:cntr, 1)) ); % location of max amplitude u

% if maximum is at centerline, then use the previous point
if ind == cntr
    ind = ind-1;
end

psi1 = psi(ind, 1) - psi(Ny+1-ind, 1);
psi2 = psi(ind, 2) - psi(Ny+1-ind, 2);

if maxuval<1e-6
    %try using w instead?
    disp('symmetrizing by w');
    [maxwval, ind] = max( psi((2*Ny+1):(2*Ny+cntr), 1) ); % location of max amplitude w
    if maxwval ~= psi(2*Ny+ind,1)
        keyboard
    end
    psi1 = psi(2*Ny+ind, 1) - psi(3*Ny+1-ind, 1);
    psi2 = psi(2*Ny+ind, 2) - psi(3*Ny+1-ind, 2);
end

r1 = 1/sqrt(1 + abs(psi1/psi2)^2);
r2 = sqrt(1 - r1^2);

phase1 = angle(psi1);
phase2 = angle(psi2);

theta2 = 0;
theta1 = theta2 + phase2 - phase1 + pi;

a1 = r1*exp(1i*theta1);
a2 = r2*exp(1i*theta2);

% make first mode in pair even, second mode odd about center
psi_sym(:, 1) = a1*psi(:,1) + a2*psi(:,2);


psi_sym(     (1:cntr)   , 2) =  psi_sym(     (1:cntr)   , 1);
psi_sym(  Ny+(1:cntr)   , 2) =  psi_sym(  Ny+(1:cntr)   , 1);
psi_sym(2*Ny+(1:cntr)   , 2) =  psi_sym(2*Ny+(1:cntr)   , 1);
psi_sym(3*Ny+(1:cntr)   , 2) =  psi_sym(3*Ny+(1:cntr)   , 1);

psi_sym(     (cntr+1:Ny), 2) = -psi_sym(     (cntr+1:Ny), 1);
psi_sym(  Ny+(cntr+1:Ny), 2) = -psi_sym(  Ny+(cntr+1:Ny), 1);
psi_sym(2*Ny+(cntr+1:Ny), 2) = -psi_sym(2*Ny+(cntr+1:Ny), 1);
psi_sym(3*Ny+(cntr+1:Ny), 2) = -psi_sym(3*Ny+(cntr+1:Ny), 1);

psi_sym([cntr Ny+cntr 2*Ny+cntr 3*Ny+cntr],2)=0;

% do the same for phi
phi_sym(:,1) = a1*phi(:,1) + a2*phi(:,2);


%how did we get away with this before???
phi_sym(     (1:cntr)   , 2) =  phi_sym(     (1:cntr)   , 1);
phi_sym(  Ny+(1:cntr)   , 2) =  phi_sym(  Ny+(1:cntr)   , 1);
phi_sym(2*Ny+(1:cntr)   , 2) =  phi_sym(2*Ny+(1:cntr)   , 1);
phi_sym(3*Ny+(1:cntr)   , 2) =  phi_sym(3*Ny+(1:cntr)   , 1);

phi_sym(     (cntr+1:Ny), 2) = -phi_sym(     (cntr+1:Ny), 1);
phi_sym(  Ny+(cntr+1:Ny), 2) = -phi_sym(  Ny+(cntr+1:Ny), 1);
phi_sym(2*Ny+(cntr+1:Ny), 2) = -phi_sym(2*Ny+(cntr+1:Ny), 1);
phi_sym(3*Ny+(cntr+1:Ny), 2) = -phi_sym(3*Ny+(cntr+1:Ny), 1);

phi_sym([cntr Ny+cntr 2*Ny+cntr 3*Ny+cntr],2)=0;


% set phase to zero at peak of (rescaled) u
if 1
    val1 = max( psi_sym(1:cntr, 1) );
    angle1 = exp(-1i*angle(val1));
    
    psi_sym(:, 1) = psi_sym(:, 1)*angle1;
    phi_sym(:, 1) = phi_sym(:, 1)*angle1;
    
    
    % do the same for the second mode
    val2 = max( psi_sym(1:cntr, 2) );
    angle2 = exp(-1i*angle(val2));
    
    psi_sym(:, 2) = psi_sym(:, 2)*angle2;
    phi_sym(:, 2) = phi_sym(:, 2)*angle2;
end

testorthopsi=psi_sym' * (sqW.^2) * psi_sym;
testorthophi=phi_sym' * (sqW.^2) * phi_sym;
if max(abs(testorthopsi-eye(2)),[],'all')>1e-5 || max(abs(testorthophi-eye(2)),[],'all')>1e-5
    if ~(kx==0 || kz==0)
     % keyboard;
    end
end







%make sure modes conserve mass and are symmetric

% flipvec = @(vec) [vec(Ny:-1:1); -vec((2*Ny):-1:(Ny+1)); vec((3*Ny):-1:(2*Ny+1)); vec((4*Ny):-1:(3*Ny+1))];

% psisymdiff=zeros(size(psi_sym));
% psisymdiff(:,1)=psi_sym(:,1)-flipvec(psi_sym(:,1));
% psisymdiff(:,2)=psi_sym(:,2)+flipvec(psi_sym(:,2));
% if norm(psisymdiff(:,1))>1e-5 || norm(psisymdiff(:,2))>1e-5
%     psisymdiff(:,1)=psi_sym(:,1)+flipvec(psi_sym(:,1));
%     psisymdiff(:,2)=psi_sym(:,2)-flipvec(psi_sym(:,2));
%     if norm(psisymdiff(:,1))>1e-5 || norm(psisymdiff(:,2))>1e-5
%         keyboard
%     end
% end
% 
% phisymdiff=zeros(size(phi_sym));
% phisymdiff(:,1)=phi_sym(:,1)-flipvec(phi_sym(:,1));
% phisymdiff(:,2)=phi_sym(:,2)+flipvec(phi_sym(:,2));
% if norm(phisymdiff(:,1))>1e-5 || norm(phisymdiff(:,2))>1e-5
%     phisymdiff(:,1)=phi_sym(:,1)+flipvec(phi_sym(:,1));
%     phisymdiff(:,2)=phi_sym(:,2)-flipvec(phi_sym(:,2));
%     if norm(phisymdiff(:,1))>1e-5 || norm(phisymdiff(:,2))>1e-5
%         keyboard
%     end
% end


divorig_psi1=1i*kx*psi(1:Ny,1)+D1*psi((Ny+1):2*Ny,1)+1i*kz*psi((2*Ny+1):(3*Ny),1);
divorig_psi2=1i*kx*psi(1:Ny,2)+D1*psi((Ny+1):2*Ny,2)+1i*kz*psi((2*Ny+1):(3*Ny),2);
divorig_phi1=1i*kx*phi(1:Ny,1)+D1*phi((Ny+1):2*Ny,1)+1i*kz*phi((2*Ny+1):(3*Ny),1);
divorig_phi2=1i*kx*phi(1:Ny,2)+D1*phi((Ny+1):2*Ny,2)+1i*kz*phi((2*Ny+1):(3*Ny),2);

divsym_psi1=1i*kx*psi_sym(1:Ny,1)+D1*psi_sym((Ny+1):2*Ny,1)+1i*kz*psi_sym((2*Ny+1):(3*Ny),1);
divsym_psi2=1i*kx*psi_sym(1:Ny,2)+D1*psi_sym((Ny+1):2*Ny,2)+1i*kz*psi_sym((2*Ny+1):(3*Ny),2);
divsym_phi1=1i*kx*phi_sym(1:Ny,1)+D1*phi_sym((Ny+1):2*Ny,1)+1i*kz*phi_sym((2*Ny+1):(3*Ny),1);
divsym_phi2=1i*kx*phi_sym(1:Ny,2)+D1*phi_sym((Ny+1):2*Ny,2)+1i*kz*phi_sym((2*Ny+1):(3*Ny),2);

divvec=[norm(divorig_psi1([1:((Ny-1)/2) ((Ny+3)/2):end])) norm(divorig_psi2([1:((Ny-1)/2) ((Ny+3)/2):end])) ...
        norm(divorig_phi1([1:((Ny-1)/2) ((Ny+3)/2):end])) norm(divorig_phi2([1:((Ny-1)/2) ((Ny+3)/2):end])) ...
         norm(divsym_psi1([1:((Ny-1)/2) ((Ny+3)/2):end]))  norm(divsym_psi2([1:((Ny-1)/2) ((Ny+3)/2):end])) ...
         norm(divsym_phi1([1:((Ny-1)/2) ((Ny+3)/2):end]))  norm(divsym_phi2([1:((Ny-1)/2) ((Ny+3)/2):end]))];

if any(divvec>1e-6)
    divvec
    keyboard
end


if 0
    close all; figure(); subplot(1,4,1) 
    plot(real(psi_sym(1:Ny,1)),y,'b'); hold on; plot(real(psi_sym(1:Ny,2)),y,'r');
    plot(real(psi(1:Ny,1)),y,'b:'); hold on; plot(real(psi(1:Ny,2)),y,'r:');
    title('Response u');
    subplot(1,4,2);
    plot(real(psi_sym((2*Ny+1):(3*Ny),1)),y,'b'); hold on; plot(real(psi_sym((2*Ny+1):(3*Ny),2)),y,'r');
    plot(real(psi((2*Ny+1):(3*Ny),1)),y,'b:'); plot(real(psi((2*Ny+1):(3*Ny),2)),y,'r:');
    title('Response w');
    subplot(1,4,3) 
    plot(real(phi_sym(1:Ny,1)),y,'b'); hold on; plot(real(phi_sym(1:Ny,2)),y,'r');
    plot(real(phi(1:Ny,1)),y,'b:'); hold on; plot(real(phi(1:Ny,2)),y,'r:');
    title('Forcing u');
    subplot(1,4,4);
    plot(real(phi_sym((2*Ny+1):(3*Ny),1)),y,'b'); hold on; plot(real(phi_sym((2*Ny+1):(3*Ny),2)),y,'r');
    plot(real(phi((2*Ny+1):(3*Ny),1)),y,'b:'); hold on; plot(real(phi((2*Ny+1):(3*Ny),2)),y,'r:');
    title('Forcing w');
end








end  %end function
