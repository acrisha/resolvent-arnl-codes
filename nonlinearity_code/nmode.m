function [u, v, w] = nmode(vec,ix,iz,is,Ny)
%retrieve the specified mode from a large matrix containing all modes

aix = abs(ix);
aiz = abs(iz);

if (ix >=0) && (iz >= 0)
    u = squeeze(vec(aix+1,aiz+1,1:Ny,is));
    v = squeeze(vec(aix+1,aiz+1,Ny+(1:Ny),is));
    w = squeeze(vec(aix+1,aiz+1,2*Ny+(1:Ny),is));
elseif (ix >= 0) && (iz < 0)
    u =  squeeze(vec(aix+1,aiz+1,1:Ny,is));
    v =  squeeze(vec(aix+1,aiz+1,Ny+(1:Ny),is));
    w = -squeeze(vec(aix+1,aiz+1,2*Ny+(1:Ny),is));
elseif (ix < 0) && (iz >= 0)
    u = conj( squeeze(vec(aix+1,aiz+1,1:Ny,is)));
    v = conj( squeeze(vec(aix+1,aiz+1,Ny+(1:Ny),is)));
    w = conj(-squeeze(vec(aix+1,aiz+1,2*Ny+(1:Ny),is)));   
elseif (ix < 0) && (iz < 0)
    u = conj(squeeze(vec(aix+1,aiz+1,1:Ny,is)));
    v = conj(squeeze(vec(aix+1,aiz+1,Ny+(1:Ny),is)));
    w = conj(squeeze(vec(aix+1,aiz+1,2*Ny+(1:Ny),is)));   
end

