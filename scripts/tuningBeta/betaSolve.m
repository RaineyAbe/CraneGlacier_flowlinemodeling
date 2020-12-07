%% solve the stress balance equations to obtain the basal roughness factor (beta)
    % Adapted and adjusted from Ellyn's original code and Enderlin et al. (2013)
    
function [beta,U,vm,T,dUdx] = betaSolve(x,h,rho_i,g,m,n,A,H,W,dx,c,ice_end,N,i,U0,x0,t,beta,E,rho_sw)
    
    % Load observed speed, use linear trendline
    U = interp1(x0,U0,x); 
        
    % Set up U, H, x, & dUdx on staggered grid using forward differencing
    Um = (U(1:end-1) + U(2:end))./2; 
        % Make sure U is a row vector
        if size(Um)==[ice_end,1]
            Um=Um';
        end 
        Um = [Um,0];
    Hm0 = (H(1:c-1) + H(2:c))./2; 
        Hm0 = [Hm0,0];
        Hm0(Hm0<0) = 0; %cannot have negative values 
        Hm=Hm0; 
    xm = (x(1:end-1)+x(2:end))./2;
        xm = [xm,(xm(end)+(xm(end)-xm(end-1)))];
        
    % Calculate dhdx, dUdx, dUmdx using forward differencing
    dhdx = zeros(1,length(x));
    dUdx = zeros(1,length(x));
    dUmdx = zeros(1,length(x));
    for k=2:length(x)
        dhdx(k) = (h(k)-h(k-1))/(x(k)-x(k-1));
        dUdx(k) = (U(k)-U(k-1))/(x(k)-x(k-1));
        dUmdx(k) = (Um(k)-Um(k-1))/(xm(k)-xm(k-1));
    end 
        
    % calculate effective viscosity required for 
    % inversion of the stress coefficient matrix
    if n == 3
        gamma = ones(1,length(x));
        for k=1:length(x)
            gamma(k) = U(k).^((1-n)/n); %linearization term for lateral resistance
        end
        gamma(1) = gamma(2); %set linearization term at the divide 
        gamma(gamma>1e+06) = 1e+06; %set the limit so gamma does not approach infinity (minimum U = 1e-09 m s^-1)
        
        % get the rate factor effective viscosity on the staggered grid for the longitudinal stress calculation
        A = A(1:c); % resize so length(A) == length(x) - end at calving front
        Am = (A(1:end-1) + A(2:end))./2; %rate factor on the staggered grid from forward differencing
        Am = [Am 0];
    
        vm = ones(1,length(x));
        for k=1:length(x)
            vm(k) = ((E*Am(k))^(-1/n))*(abs(dUmdx(k)))^((1-n)/n); %effective viscosity (Pa s)
        end
        vm(1) = vm(2);
        vm(end) = 0; % Avoids infinite value at boundary (Am = 0 at terminus)
        vm(vm>8e+16) = 8e+16; %set a maximum value for very low strain rates
        %vm(1:2) = vm(1:2)*10; % Try increasing viscosity at higher altitude
        %vm = 2e14.*ones(1,length(x));  % viscosity of ice at -10 deg C

        if m > 1
            eta = ones(1,length(x));
            for k=1:length(x)
                eta(k) = U(k)^((1-m)/m); %linearization term for basal resistance
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
    
    % Set up the speed and stress terms to solve for the linearized
    % velocity coefficient terms
        % [G_minus(k)*U(k-1)+G(k)*U(k)+G_plus(k)*U(k+1)=T] 
            
    % Solve for the speed coefficients G_minus and G_plus 
    %   and the driving stress T
        G_minus = zeros(1,c); 
        G_plus = zeros(1,c);
        T = zeros(1,c);
        
    % set the upper boundary coefficients
            G_minus(1) = 0;
            G(1) = 1;
            G_plus(1) = 0;
            T(1) = 0;  
    % solve for the coefficients
        for k=2:c
            G_minus(k) = (2/(dx^2))*Hm(k-1)*vm(k-1); 
            G_plus(k) = (2/(dx^2))*Hm(k)*vm(k); 
            T(k) = rho_i*g*H(k)*dhdx(k);
        end 
    % set the lower boundary coefficients by applying the 
    % hydrostatic equilibrium boundary condition 
        G_minus(c) = -1;
        G(c) = 1;
        G_plus(c) = 0;
        T(c) = (E*A(c).*(((rho_i.*g./4).*(H(c).*(1-(rho_i./rho_sw)))).^n)).*dx;
        
    %remove any NaNs from the coefficient vectors
    G_minus(isnan(G_minus)) = 0;
    G(isnan(G)) = 0;
    G_plus(isnan(G_plus)) = 0;
          
      % Solve for beta using the G term (solved from the equation below)
          % [G_minus(k)*U(k-1)+G(k)*U(k)+G_plus(k)*U(k+1)=T]
          % G(k) = (T(k) - G_minus(k)*U(k-1) - G_plus(k)*U(k+1)))/U(k)
          %   where G(k) = -2./(dx(k).^2)*(Hm(k)*vm(k)+Hm(k-1)*vm(k-1)) 
          %             - (beta(k)*N(k)_eta(k))
          %             - (gamma(k)*H(k)/W(k))*(5/(2*A(k)*W(k))^(1/3)); 
          %
          %        beta(k) = [-G(k) - 2./(dx(k).^2)*(Hm(k)*vm(k)+Hm(k-1)*vm(k-1)) 
          %         - (gamma(k)*H(k)/W(k))*(5/(2*A(k)*W(k))^(1/3))]/(N(k)*eta(k));
              
        % Solve for G (for U(k))
        G = []; 
        for k = 2:length(T)
            G(k) = (T(k)-G_minus(k)*Um(k-1)-G_plus(k)*Um(k))/U(k);
        end 
        G(1) = G(2);
        
    %apply the hydrostatic equilibrium boundary condition from the calving
    %front (c) to the end of the ice-covered domain (terminus)
        G_minus(c) = -1;
        G(c) = 1;
        G_plus(c) = 0;
        
        % remove any NaNs from the coefficient vectors
        G_minus(isnan(G_minus)) = 0;
        G_plus(isnan(G_plus)) = 0; 
        
    % Solve for beta (basal roughness factor)
        beta(i,1)=0; % 
        for k=2:c            
            beta(i,k) = (-G(k) - (2/(dx^2))*(Hm(k)*vm(k)+Hm(k-1)*vm(k-1))...
                - (2*gamma(k)*H(k)/W(k))*((5/(A(k)*W(k)))^(1/n)))/(N(k)*eta(k));            
        end 
        beta(i,1) = beta(i,2); %beta(i,c) = 0; % floating
        %beta(beta<=0)=0;
   
    % Run the forward U_convergence with the resulting beta to check success
        %set-up coefficient vectors for the linearized stress terms over the calving front
        %[C(k)*U(k-1)+E(k)*U(k)+G(k)*U(k+1)=Td]
        for k=2:c
            G_minus_rev(k) = (2/(dx^2))*Hm(k-1)*vm(k-1); 
            G_rev(k) = (-2/(dx^2))*(Hm(k-1)*vm(k-1)+Hm(k)*vm(k))-...
                (beta(i,k)*N(k)*eta(k))-...
                ((2*gamma(k)*H(k)/W(k))*((5/(E*A(k)*W(k)))^(1/n)));
            G_plus_rev(k) = (2/(dx^2))*Hm(k)*vm(k);
            T_rev(k) = (rho_i*g*H(k)*dhdx(k));                
        end    

        % Set boundary conditions
        G_minus_rev(1)=G_minus_rev(2); G_rev(1)=G_rev(2); 
        G_plus_rev(1)=G_plus_rev(2);T_rev(1)=T_rev(2); T_rev(c) = 0;
        
        %remove any NaNs from the coefficient vectors
            G_minus_rev(isnan(G_minus_rev)) = 0;
            G_rev(isnan(G_rev)) = 0;
            G_plus_rev(isnan(G_plus_rev)) = 0;

            %create a sparse tri-diagonal matrix for the velocity coefficient vectors
            M = diag(G_minus_rev(2:ice_end),-1) + diag(G_rev(1:ice_end)) + diag(G_plus_rev(1:ice_end-1),1);
            M(1,1) = 1; M(c,c-1) = -1; M(c,c) = 1;
            M = sparse(M);
            
            %make sure Td is a column vector for the inversion 
            if size(T_rev)==[1,c]
                T_rev=T_rev';
            end

            %use the backslash operator to perform the matrix inversion to solve for ice velocities
            Un_rev = M\T_rev; %velocity (m s^-1)

            %remove NaNs
            Un_rev(isnan(Un_rev)) = 0;

            %set velocity at the terminus
            Un_rev(c) = 0;

            %calculate new strain rates (equal to the velocity gradient)
            dUndx_rev = gradient(Un_rev,x); %strain rate (s^-1)

            %make sure Un is a row vector so it can be compared with U
            if size(Un_rev) == [length(x),1]
                Un_rev=Un_rev';
            end

            %make sure dUdx is a row vector
            if size(dUndx_rev) == [length(x),1]
                dUndx_rev=dUndx_rev';
            end

            %check if the difference in speed between iteratons (U vs. Un_rev) meets a set tolerance
            if abs(sum(U)-sum(Un_rev))<0.1*abs(sum(U)) %determine if U has converged sufficiently
                %use sufficiently converged values for speeds & strain rates
                   disp(['t = ',num2str(t(i)/3.1536e7),' yrs: U did not converge sufficiently.']);
            end
            
end 