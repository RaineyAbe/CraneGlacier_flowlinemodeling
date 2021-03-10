%% solve the stress balance equations to obtain the basal roughness factor (beta)
% Adapted by Rainey Aberle from Ellyn Enderlin's flowline model (Enderlin et al., 2013)
    
function [beta,dUdx,Un_rev] = betaSolve(H,c,x,U,n,A,E,m,dx,rho_i,g,h,rho_sw,sigma_b,W,N)
        
    % Set up H & dHdx on staggered grid
    Hm(1:c-1) = (H(2:c) + H(1:c-1))./2; % forward difference
    Hm(c) = (H(c)+H(c-1))./2; % backward difference at c
    dHmdx(1:c-1) = (H(2:c)-H(1:c-1))./(x(2:c)-x(1:c-1)); % forward difference
    dHmdx(c) = (H(c)-H(c-1))./(x(c)-U(c-1)); % backward difference at c    
        
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
    % coefficients up to calving front
    G_minus(2:c-1) = (2./(dx(2:c-1).^2)).*Hm(1:c-2).*vm(1:c-2); %for U(k-1)
    G_plus(2:c-1) = (2./(dx(2:c-1).^2)).*Hm(2:c-1).*vm(2:c-1); %for U(k+1)
    T(2:c-1) = (rho_i.*g.*H(2:c-1).*(h(1:c-2)-h(3:c))./(x(1:c-2)-x(3:c))); %gravitational driving stress
    % upper boundary condition
    T(1) = (rho_i.*g.*H(1).*(h(1)-h(2))./(x(1)-x(2)));      
    % calving front condition
    G_minus(c) = -1;
    G_plus(c) = 0;
    T(c) = (E(c)*A(c).*(((rho_i.*g./4).*((H(c).*(1-(rho_i./rho_sw))-sigma_b./(rho_i.*g)))).^n)).*dx(c); %SOMETHING WRONG HERE
    %remove any NaNs from the coefficient vectors
    G_minus(isnan(G_minus)) = 0;
    G_plus(isnan(G_plus)) = 0;
    T(isnan(T)) = 0;
          
    % Solve for beta using the G term (solved from the equation below)
    % [G_minus(k)*U(k-1)+G(k)*U(k)+G_plus(k)*U(k+1)=T]
    % G(k) = (T(k) - G_minus(k)*U(k-1) - G_plus(k)*U(k+1)))/U(k)
    %   where G(k) = -2./(dx(k).^2)*(Hm(k)*vm(k)+Hm(k-1)*vm(k-1)) 
    %             - (beta(k)*N(k)*eta(k))
    %             - (gamma(k)*H(k)/W(k))*(5/(2*A(k)*W(k))^(1/3)); 
    %
    %        beta(k) = [-G(k) - 2./(dx(k).^2)*(Hm(k)*vm(k)+Hm(k-1)*vm(k-1)) 
    %         - (gamma(k)*H(k)/W(k))*(5/(2*A(k)*W(k))^(1/3))]/(N(k)*eta(k));
    % Solve for G (for U(k))
    G(2:c) = (T(2:c)-G_minus(2:c).*Um(1:c-1)-G_plus(2:c).*Um(2:c))./U(2:c);
    G(1)=G(2);
        
    % Solve the basal roughness factor, beta
    beta(2:c) = (-G(2:c)-(2./(dx(2:c).^2)).*(Hm(2:c).*vm(2:c)+Hm(1:c-1).*vm(1:c-1))...
        -(2.*gamma(2:c).*H(2:c)./W(2:c)).*((5./(A(2:c).*W(2:c))).^(1/n)))./(N(2:c).*eta(2:c));            
    beta(1) = beta(2); beta(c)=0;
    beta(beta<0)=0; % beta cannot be less than 0
   
    % Run the forward U_convergence with the resulting beta to check success
        %set-up coefficient vectors for the linearized stress terms over the calving front
        %[C(k)*U(k-1)+E(k)*U(k)+G(k)*U(k+1)=Td]  
        % upper boundary condition
        G_minus_rev(1) = 0;
        G_rev(1) = -(beta(1).*(N(1)/(rho_i*g)).*eta(1))-...
                (((2*gamma(1).*H(1))./W(1)).*((5/(E(1)*A(1).*W(1))).^(1/n))); 
        G_plus_rev(1) = 0;
        T_rev(1) = (rho_i.*g.*H(1).*(h(1)-h(2))./(x(1)-x(2)));    
        % coefficients up to calving front
        G_minus_rev(2:c-1) = (2./(dx(2:c-1).^2)).*Hm(1:c-2).*vm(1:c-2); %for U(k-1)
        G_rev(2:c-1) = (-2./(dx(2:c-1).^2)).*(Hm(1:c-2).*vm(1:c-2)+Hm(1:c-2).*vm(1:c-2))-...
            (beta(2:c-1).*(N(2:c-1)).*eta(2:c-1))-...
            (((2*gamma(2:c-1).*H(2:c-1))./W(2:c-1)).*((5./(E(2:c-1).*A(2:c-1).*W(2:c-1))).^(1/n))); %for U(k)        
        G_plus_rev(2:c-1) = (2./(dx(2:c-1).^2)).*Hm(2:c-1).*vm(2:c-1); %for U(k+1)
        T_rev(2:c-1) = (rho_i.*g.*H(2:c-1).*(h(3:c)-h(1:c-2))./(x(3:c)-x(1:c-2))); %gravitational driving stress  
        % calving front condition
        G_minus_rev(c) = -1;
        G_rev(c) = 1;
        G_plus_rev(c) = 0;
        T_rev(c) = (E(c)*A(c).*(((rho_i.*g./4).*((H(c).*(1-(rho_i./rho_sw))-sigma_b./(rho_i.*g)))).^n)).*dx(c);
        %remove any NaNs from the coefficient vectors
        G_minus_rev(isnan(G_minus_rev)) = 0;
        G_rev(isnan(G_rev)) = 0;
        G_plus_rev(isnan(G_plus_rev)) = 0;
        T_rev(isnan(T_rev)) = 0;

        %create a sparse tri-diagonal matrix for the velocity coefficient vectors
        M = diag(G_minus_rev(2:c),-1) + diag(G_rev(1:c)) + diag(G_plus_rev(1:c-1),1);
        M = sparse(M);
            
        %make sure Td is a column vector for the inversion 
        if size(T_rev)==[1,c]
            T_rev=T_rev';
        end

        %use the backslash operator to perform the matrix inversion to solve for ice velocities
        Un_rev = M\T_rev; %velocity (m s^-1)
        %remove NaNs
        Un_rev(isnan(Un_rev)) = 0;
        Un_rev(Un_rev<0)=0; % U cannot be less than zero
        
        %make sure Un is a row vector so it can be compared with U
        if size(Un_rev) == [c,1]
            Un_rev=Un_rev';
        end
        
        %calculate new strain rates (equal to the velocity gradient)
        dUndx_rev(1) = (Un_rev(2)-Un_rev(1))/(x(2)-x(1)); % backward difference
        dUndx_rev(2:c-1) = (Un_rev(3:c)-Un_rev(1:c-2))./(x(3:c)-x(1:c-2)); % central difference
        dUndx_rev(c) = (Un_rev(c)-Un_rev(c-1))./(x(c)-x(c-1)); % forward differe ce

        %make sure dUdx is a row vector
        if size(dUndx_rev) == [c,1]
            dUndx_rev=dUndx_rev';
        end

        %check if the difference in speed between iteratons (U vs. Un_rev) meets a set tolerance
        if abs(sum(U)-sum(Un_rev))<0.1*abs(sum(U))==0 %determine if U has converged sufficiently
            %use sufficiently converged values for speeds & strain rates
               disp('U did not converge sufficiently.');
        end
        
        figure(3); clf
        set(gcf,'Position',[1   350   560   420]);
        subplot(2,1,1);
            set(gca,'fontsize',12,'fontweight','bold');
            hold on; legend('Location','southeast'); grid on; title('G terms');
            plot(x,G,'b','displayname','fwd','linewidth',2); 
            plot(x,G_plus,'b','handlevisibility','off','linewidth',2); 
            plot(x,G_minus,'b','handlevisibility','off','linewidth',2);
            plot(x,G_rev,'m','displayname','rev','linewidth',2); 
            plot(x,G_plus_rev,'m','handlevisibility','off','linewidth',2); 
            plot(x,G_minus_rev,'m','handlevisibility','off','linewidth',2);
        subplot(2,1,2);
            set(gca,'fontsize',12,'fontweight','bold');
            hold on; legend('Location','southeast'); grid on; title('T');
            plot(x,T_rev,'m','displayname','rev','linewidth',2); 
            plot(x,T,'b','displayname','fwd','linewidth',2);
            
end 