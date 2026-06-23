% validation script for nonlinearity code.
% script is called from comput_nonlinearity_omeaga_val_opt_june2026.m

%VALIDATION
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
