%% solve the stress balance equations to obtain speed values (U)
function [U,dUdx,vm,T] = U_convergence_varyingE(x,U,U0,dUdx,dhdx,H,A,E,N,W,dx,c,ice_end,n,m,beta,rho_i,rho_sw,g)

b=1;
while b
    
    % Calculate the thickness on the staggered grid for use in stress balance equations
    Hm = (H(1:end-1) + H(2:end))./2; %use forward differences
    Hm = [Hm,0];
    Hm(Hm<0) = 0; %cannot have negative values
    
    %calculate the linearization terms & effective viscosity required for 
    %inversion of the stress coefficient matrix
    if n == 3
        gamma=zeros(1,length(x)); % pre-allocate gamma
        for k=1:length(x)
            gamma(k) = U(k).^((1-n)/n); %linearization term for lateral resistance
        end
        gamma(1) = gamma(2); %set linearization term at the divide (U(1) = 0)
        gamma(gamma>1e+06) = 1e+06; %set the limit so gamma does not approach infinity (minimum U = 1e-09 m s^-1)
        
        %get the effective viscosity on the staggered grid for the longitudinal stress calculation
        Am = (A(1:end-1) + A(2:end))./2; %rate factor on the staggered grid from forward differencing
        Am = [Am,0];
        Um = (U(1:end-1) + U(2:end))./2; %speed on the staggered grid from forward differencing
        Um = [Um,0];
        dUmdx = [(Um(2:end)-Um(1:end-1))./(x(2:end)-x(1:end-1)) 0]; %strain rates on the staggered grid

        vm=zeros(1,length(x)); % pre-allocate vm
        for k=1:length(x)
            vm(k) = ((E(k)*Am(k)).^(-1/n)).*(abs(dUmdx(k))).^((1-n)/n); %effective viscosity (Pa s)
        end
        vm(1) = vm(2); 
        vm(vm>8e+16) = 8e+16; %set a maximum value for very low strain rates
                
        if m > 1
            eta=zeros(1,length(x)); % pre-allocate eta
            for k=1:length(x)
                eta(k) = U(k).^((1-m)/m); %linearization term for basal resistance
            end
            eta(1) = eta(2); %set linearization term at the divide (U(1) = 0)
           
            %set the limit so eta does not approach infinity (minimum U = 1e-09 m s^-1)
            if m == 2
                eta(eta>3.16e+04) = 3.16e+04;
            end
            if m == 3
                eta(eta>1e+06) = 1e+06;
            end
        else
            eta = ones(1,length(x)); %if m=1, the basal resistance term does not need to be linearized
        end
    else
        disp(['Adjust maximum value for the lateral resistance linearization term (gamma)']);
    end
    
    %set-up coefficient vectors for the linearized stress terms over the calving front
    %[C(k)*U(k-1)+E(k)*U(k)+G(k)*U(k+1)=Td]
   for k=1:c-1
        if k == 1
            %set the upper boundary coefficients
            G_minus(1) = 0;
            G(1) = 1;
            G_plus(1) = 0;
        else
            G_minus(k) = (2./(dx.^2)).*Hm(k-1).*vm(k-1); %for U(k-1)
            G(k) = (-2./(dx.^2)).*(Hm(k-1).*vm(k-1)+Hm(k).*vm(k))-...
                (beta(k).*(N(k)/(rho_i*g)).*eta(k))-...
                (((2*gamma(k).*H(k))./W(k)).*((5/(E(k)*A(k).*W(k))).^(1/n))); %for U(k)
            G_plus(k) = (2./(dx.^2)).*Hm(k).*vm(k); %for U(k+1)
            T(k) = (rho_i.*g.*H(k).*dhdx(k)); %gravitational driving stress
        end
   end
    % use driving stress at the upper boundary that results in 90% observed
    % velocity 
    T(1) = U0(1)*0.9.*G(1);
    
    %apply the hydrostatic equilibrium boundary condition from the calving
    %front (c) to the end of the ice-covered domain (terminus)
    for j=1:(ice_end-c)+1
        G_minus(c+(j-1)) = -1;
        G(c+(j-1)) = 1;
        G_plus(c+(j-1)) = 0;
        T(c+(j-1)) = (E(k)*A(c+(j-1)).*(((rho_i.*g./4).*(H(c+(j-1)).*(1-(rho_i./rho_sw)))).^n)).*dx;
    end
    T(ice_end) = 0;
        
    %remove any NaNs from the coefficient vectors
    G_minus(isnan(G_minus)) = 0;
    G(isnan(G)) = 0;
    G_plus(isnan(G_plus)) = 0;
    
    %create a sparse tri-diagonal matrix for the velocity coefficient vectors
    M = diag(G_minus(2:ice_end),-1) + diag(G(1:ice_end)) + diag(G_plus(1:ice_end-1),1);
    M = sparse(M);
    
    %make sure Td is a column vector for the inversion 
    if size(T) == [1,ice_end]
        T=T';
    end
    
    %use the backslash operator to perform the matrix inversion to solve for ice velocities
    Un = M\T(1:ice_end); %velocity (m s^-1)
    
    %remove NaNs and apply the ice divide bounday condition
    Un(isnan(Un)) = 0;
    Un(Un<0) = 0;
    
    %set velocity at the (dummy) terminus
    Un(ice_end+1:length(x)) = Un(c)*ones(1,length(x(ice_end+1:length(x))));
        
    %make sure Un is a row vector so it can be compared with U
    if size(Un) == [length(x),1]
        Un=Un';
    end    
    
    %calculate new strain rates (equal to the velocity gradient)
    dUndx = [(Un(2:end)-Un(1:end-1))./(x(2:end)-x(1:end-1)) 0]; %strain rate (s^-1)
    
    %make sure dUdxn is a row vector
    if size(dUndx) == [length(x),1]
        dUndx=dUndx';
    end
    
    %check if the difference in speed between iteratons (U vs. Un) meets a set tolerance
    if abs(sum(U)-sum(Un))<0.1*abs(sum(U)) %determine if U has converged sufficiently
        %use sufficiently converged values for speeds & strain rates
            U = Un; 
            dUdx = dUndx;
            return %break the U iterations
    end
    
    %if not sufficiently converged, set Un to U and solve the stress balance matrix again
    U = Un;
    dUdx = dUndx;
    
    %terminate the U iteration loop if convergence takes too long
    if str2double(num2str(b)) > str2double(num2str(50))
        b = b;
        return
    end
            
    %loop through
    b = b+1;
    
end

end

