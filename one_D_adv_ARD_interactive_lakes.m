function [T_np1_k_j,phi_np1_k_j,S_np1_k_j]=one_D_adv_ARD_interactive_lakes(...
    k_i,k_br,k_s,dt,dz,H,c_i,c_br,L,STol,TTol,PhiTol,phi_initial,T_initial,...
    S_initial,S_bottom,rho_i,rho_br,Tm,m,Ttop,Tbottom,g,rho_sw,phi_c,BCoc,BCdn,dT,Lake)

%% Initial Condition Vectors for Phi, T, S, w
phi_initial=phi_initial';
T_initial=T_initial';
phi_np1_km2_j=phi_initial;
T_np1_km1_j=T_initial;
S_initial=S_initial';
S_np1_km1_j=S_initial;
%SinI=SinI_initial;

S_evolve=[]';
T_evolve=[]';
phi_evolve=[]';
SErr=STol+1;
TErr=TTol+1;
PhiErr=PhiTol+1;
counter=0;
matrix_dimension=H/dz+1;


%% Itterattions over k
while SErr>STol || TErr>TTol || PhiErr>PhiTol;
    counter=counter+1;
    if counter>25
        disp(counter)
    else
    end

    %% Flux from gravity drainage parameterization
    [T_source,S_source]=one_D_advective_flux_lakes(phi_np1_km2_j,S_np1_km1_j...
        ,T_np1_km1_j,H,dz,dt,S_bottom,g,rho_sw);
    
    %% Freezing point depression calculations
    %   delta_T=1.853*(S_np1_km1_j/28);        %% Old Freezing point depression
    %   due to salt (linear approx)
    
    if Lake==1
        delta_T=Tm-(-(9.1969758*(10^-5)*S_np1_km1_j.^2)-0.03942059*S_np1_km1_j+...
       272.63617665);                        %% Liquidus point using FREZCHEM seawater
    elseif Lake==2
        delta_T=Tm-(-(1.333489497*(10^-5)*S_np1_km1_j.^2)-0.01612951864*S_np1_km1_j+...
            273.055175687);                       %% Liquidus point using FREZCHEM Europa ocean

% BC LAKES LIQUIDUS CURVES
    elseif Lake==3
    delta_T=Tm-(-(2.1302*(10^-5)*S_np1_km1_j.^2)-0.01307*S_np1_km1_j+...
       273.0403);                        %% Liquidus point using FREZCHEM Salt Lake
    elseif Lake==4
    delta_T=Tm-(-(2.9486*(10^-5)*S_np1_km1_j.^2)-0.012162*S_np1_km1_j+...
       273.00);                        %% Liquidus point using FREZCHEM Basque Lake 2.2 
    elseif Lake==5
    delta_T=Tm-(-(3.2684*(10^-5)*S_np1_km1_j.^2)-0.010903*S_np1_km1_j+...
       273.0587);                        %% Liquidus point using FREZCHEM Basque Lake 2.3
    else
    delta_T=Tm-(-(1.0579*(10^-5)*S_np1_km1_j.^2)-0.033321*S_np1_km1_j+...
       272.8703);                        %% Liquidus point using FREZCHEM Last Chance Lake 
    end

    Hs=c_i*(Tm-delta_T);                     %% Enthalpy Calculation
    Hsmelt=c_i*Tm;                           %% Enthalpy to melt ice

    %% Calculating solid fraction at k-1 iteration for use in
    %% temperature, salinity, and brine velocity at k iteration
    
    %% Value of Phi(n+1,k-1,j)
    phi_np1_km1_j=phi_np1_km2_j;
    for i=1:matrix_dimension;
    %% Test to see if melting ice or freezing brine
    En=c_i*T_np1_km1_j(i)+L*phi_np1_km2_j(i);
        if En<Hs(i)+phi_c*L
            phi_np1_km1_j(i)=phi_c;
        elseif En>Hs(i)+L;
            phi_np1_km1_j(i)=1;
        else
            phi_np1_km1_j(i)=((c_i*T_np1_km1_j(i)+phi_np1_km2_j(i)*L)-Hs(i))/...
            L;
        end;
    end;

    %% Reassigning Phi(n+1,k-2) as Phi(n+1,k-1) for next
    %% iteration
    phi_np1_km2_j=phi_np1_km1_j;
    
    
    %% Create K matrix and solve for T
    D=((phi_np1_km1_j*k_br+(1-phi_np1_km1_j)*k_i)*dt)./...
        ((dz^2)*(phi_np1_km1_j*(rho_br*c_br)+(1-phi_np1_km1_j)...
        *(rho_i*c_i)));
    for i=1:matrix_dimension;
            a(i)=1+2*D(i);
            b(i)=-D(i);
            c(i)=-D(i);
    end        

    b(end)=[];
    c(1)=[];
    y=T_initial-((rho_i*L./((phi_np1_km1_j*(rho_br*c_br)+...
        (1-phi_np1_km1_j)*(rho_i*c_i)))).*...
        (phi_np1_km1_j-phi_initial))+(1/((phi_np1_km1_j*(rho_br*c_br)+...
        (1-phi_np1_km1_j)*(rho_i*c_i))))'.*T_source;

    y(1)=y(1)+D(1)*Ttop;
    
    if BCdn==1
        y(end)=y(end)+D(end)*Tbottom;
    else
        c(end)=c(end)-D(end-1);
        y(end)=y(end)+2*D(end)*dz*dT;
    end

    % Solve for T
    T_np1_k_j=Thomas_Trid(a,b,c,y');
    
    T_np1_km1_j=T_np1_k_j;
    
    %% Create K_s matrix and solve for S
    one_over_phi1(1:matrix_dimension)=1;
    one_over_phi=rdivide(one_over_phi1',phi_np1_km1_j);
    %one_over_one_minus_phi=rdivide(one_over_phi1',1-phi_np1_km1_j);
    for i=1:matrix_dimension;
        if phi_np1_km1_j(i)<=phi_c;
            one_over_phi(i)=0;
        else
        end;
    end;
%     for i=1:matrix_dimension;
%         if one_over_one_minus_phi(i)==Inf;
%             one_over_one_minus_phi(i)=0;
%         end;
%     end;
    k_s_mod1(1:matrix_dimension)=k_s;
    k_s_mod=times(k_s_mod1',phi_np1_km1_j.^m);
    D_s_mod=((k_s_mod*dt)/dz^2).*one_over_phi;
    
    % Spatially averaged diffusion tensors
%     D_s_plus=(D_s_mod+circshift(D_s_mod,1))/2;
%     D_s_plus(1)=D_s_mod(1);
%     D_s_minus=(D_s_mod+circshift(D_s_mod,-1))/2;
%     D_s_minus(end)=D_s_mod(end);
% 
%     a2(1)=1+D_s_minus(1);
%     b2(1)=-D_s_minus(1);
%     for i=2:matrix_dimension;
%             a2(i)=1+D_s_plus(i)+D_s_minus(i);
%             b2(i)=-D_s_minus(i);
%             c2(i)=-D_s_plus(i);
%     end
    
    a2(1)=1+D_s_mod(1);
    b2(1)=-D_s_mod(1);
    for i=2:matrix_dimension;
            a2(i)=1+D_s_mod(i)+D_s_mod(i);
            b2(i)=-D_s_mod(i);
            c2(i)=-D_s_mod(i);
    end

    b2(end)=[];
    c2(1)=[];
    y2=S_initial-(rho_i/rho_br)*one_over_phi.*S_initial...
        .*(phi_np1_km1_j-phi_initial)+one_over_phi.*S_source;
%     for i=1:length(phi_np1_km1_j)
%         if phi_np1_km1_j(i)==0
%             y2(i)=0;
%             SinI(i)=one_over_one_minus_phi(i)*(SinI(i)*(1-phi_initial(i))-S_initial(i)*(phi_np1_km1_j(i)-phi_initial(i)));
%         elseif phi_np1_km1_j(i)-phi_initial(i)<1
%             y2(i)=S_initial(i)-(rho_i/rho_br)*one_over_phi(i)*(1-P_s)*S_initial(i)...
%                 *(phi_np1_km1_j(i)-phi_initial(i))+one_over_phi(i)*S_source(i);
%             SinI(i)=one_over_one_minus_phi(i)*(SinI(i)*(1-phi_initial(i))-S_initial(i)*P_s*(phi_np1_km1_j(i)-phi_initial(i)));
%         else
%             y2(i)=S_initial(i)-(rho_i/rho_br)*one_over_phi(i)*S_initial(i)...
%                 *(phi_np1_km1_j(i)-phi_initial(i))+one_over_phi(i)*SinI(i)*(phi_np1_km1_j(i)-phi_initial(i))...
%                 +one_over_phi(i)*S_source(i);
%         end
%     end

    if BCoc==1
        y2(end)=y2(end)+D_s_mod(end)*S_bottom;
    else
        a2(end)=1+D_s_mod(end);
    end


    %Solve for S
    S_np1_k_j=Thomas_Trid(a2,b2,c2,y2');
    
    S_np1_km1_j=S_np1_k_j;
    
    %% Appending value to matrix to check for convergence
    S_evolve=[S_evolve S_np1_k_j];
    T_evolve=[T_evolve T_np1_k_j];
    phi_evolve=[phi_evolve phi_np1_km1_j];
    
    %% Reassign T vector to use in next k iteration
    T_np1_km1_j=T_np1_k_j;
    phi_np1_k_j=phi_np1_km1_j;
    % tolerance tests
    if counter==1;
        SErr=SErr;
        TErr=TErr;
        PhiErr=PhiErr;
    else
        TErr=max(abs(T_evolve(:,counter)-T_evolve(:,counter-1)));
        PhiErr=max(abs(phi_evolve(:,counter)-phi_evolve(:,counter-1)));
        SErr=max(abs(S_evolve(:,counter)-S_evolve(:,counter-1)));
    end;
end



