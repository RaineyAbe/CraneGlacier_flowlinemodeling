%% solve the stress balance equations to obtain speed values (U)
function [U,dUdx,vm,Un,dUndx,M,T] = U_convergence(x,U,dUdx,H,h,A,E,N,W,dx,c,n,m,beta,rho_i,rho_sw,g,sigma_b,i)

b=1; 

while b
    
    % Calculate H on a staggered grid
    Hm(1:c-1) = (H(2:c) + H(1:c-1))./2; % forward difference
    Hm(c) = (H(c)+H(c-1))./2; % backward difference at c
    
    %calculate the linearization terms & effective viscosity required for 
    %inversion of the stress coefficient matrix
    if n == 3
        gamma=zeros(1,c); % pre-allocate gamma
        for k=1:c
            gamma(k) = U(k).^((1-n)/n); % linearization term for lateral resistance
        end
        gamma(1) = gamma(2); % set linearization term at the divide (U(1) = 0)
        gamma(gamma>1e+06) = 1e+06; % set the limit so gamma does not approach infinity (minimum U = 1e-09 m s^-1)
        
        % get A, U, & the effective viscosity on the staggered grid for the 
        % longitudinal stress calculation
        Am(1:c-1) = (A(2:c)+A(1:c-1))./2; % forward difference
        Am(c) = (A(c)+A(c-1))./2; % backward difference at c
        
        Um(1:c-1) = (U(2:c)+U(1:c-1))./2; % forward difference
        Um(c) = (U(c)+U(c-1))./2; % backward difference at c
        
        dUmdx(1:c-1) = (U(1:c-1)-U(2:c))/(x(1:c-1)-x(2:c)); % forward difference
        dUmdx(c) = (U(c-1)-U(c))/(x(c-1)-x(c)); % backward difference at c
        
        dUdx(1) = (U(2)-U(1))./(x(2)-x(1)); % forward difference
        dUdx(2:c-1) = (U(3:c)-U(1:c-2))./(x(3:c)-x(1:c-2)); % central difference
        dUdx(c) = (U(c-1)-U(c))/(x(c-1)-x(c)); % backward difference at c
        
        vm = ((E.*Am).^(-1/n)).*(abs(dUmdx)).^((1-n)/n);
        vm(vm>8e+16) = 8e+16; %set a maximum value for very low strain rates
        
        if m > 1
            eta=zeros(1,c); % pre-allocate eta
            for k=1:c
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
            eta = ones(1,c); %if m=1, the basal resistance term does not need to be linearized
        end
    else
        disp(['Adjust maximum value for the lateral resistance linearization term (gamma)']);
    end
    
    %set-up coefficient vectors for the linearized stress terms over the calving front
    %[C(k)*U(k-1)+E(k)*U(k)+G(k)*U(k+1)=Td] 
    % upper boundary condition
    G(1) = -(beta(1).*N(1).*eta(1))-...
            (((2*gamma(1).*H(1))./W(1)).*((5/(E*A(1).*W(1))).^(1/n))); 
    T(1) = U(1).*0.9/G(1); %*(rho_i.*g.*H(1).*(h(1)-h(2))./(x(1)-x(2)));     
    % coefficients up to calving front
    G_minus(2:c-1) = (2./(dx(2:c-1).^2)).*Hm(1:c-2).*vm(1:c-2); %for U(k-1)
    G(2:c-1) = (-2./(dx(2:c-1).^2)).*(Hm(1:c-2).*vm(1:c-2)+Hm(2:c-1).*vm(2:c-1))-...
            (beta(2:c-1).*N(2:c-1).*eta(2:c-1))-...
            (((2.*gamma(2:c-1).*H(2:c-1))./W(2:c-1)).*((5./(E.*A(2:c-1).*W(2:c-1))).^(1/n))); %for U(k)
    G_plus(2:c-1) = (2./(dx(2:c-1).^2)).*Hm(2:c-1).*vm(2:c-1); %for U(k+1)
    T(2:c-1) = (rho_i.*g.*H(2:c-1).*(h(1:c-2)-h(3:c))./(x(1:c-2)-x(3:c))); %gravitational driving stress     
    % calving front condition
    G_minus(c) = -1;
    G(c) = 1;
    G_plus(c) = 0;
    T(c) = (E*A(c).*(((rho_i.*g./4).*((H(c).*(1-(rho_i./rho_sw))-sigma_b./(rho_i.*g)))).^n)).*dx(c); 
    %remove any NaNs from the coefficient vectors
    G_minus(isnan(G_minus)) = 0;
    G_plus(isnan(G_plus)) = 0;
    T(isnan(T)) = 0;
    
    %create a sparse tri-diagonal matrix for the velocity coefficient vectors
    M = diag(G_minus(2:c),-1) + diag(G(1:c)) + diag(G_plus(1:c-1),1);
    M = sparse(M);
    
    %make sure Td is a column vector for the inversion 
    if size(T) == [1,c]
        T=T';
    end
    
    %use the backslash operator to perform the matrix inversion to solve for ice velocities
    Un = M\T; % velocity (m s^-1)
    Un(isnan(Un)) = 0;
    Un(Un<0) = 0;   
    
    % calculate new strain rates (1/s)
    clearvars dUndx 
    dUndx(1) = (U(1)-U(2))/(x(1)-x(2)); % forward difference
    dUndx(2:c-1) = (U(1:c-2)-U(3:c))./(x(1:c-2)-x(3:c)); % central difference
    dUndx(c) = (U(c)-U(c-1))./(x(c)-x(c-1)); % backward difference

    %make sure Un is a row vector so it can be compared with U
    if size(Un) == [c,1]
        Un=Un';
    end 
    
    %make sure dUdxn is a row vector
    if size(dUndx) == [c,1]
        dUndx=dUndx';
    end
    
    %check if the difference in speed between iteratons (U vs. Un) meets a set tolerance
    if i==1
        U = Un;
        dUdx = dUndx;
        return; % break the U iterations
    elseif abs(sum(U)-sum(Un))<0.1*abs(sum(U)) %determine if U has converged sufficiently
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
        return
    end
            
    %loop through
    b = b+1;
    
end

end

