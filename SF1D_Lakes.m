%% This is the main control file for simulating the evolution of shallow lake solidification
%  and will be used with 'one_D_Adv_ARD_interactive_lakes.m',
%  'one_D_advective_flux_lakes.m', and 'Thomas_Trid.m'. The majority of
%  modifications will be made to the input parameters and the output will
%  be Temperature, Porosity, Salinity, and Bulk Salinity profiles.
%  (variables - Temperature, Liquid_Fraction, Salinity,
%  Porosity'.*Salinity)

%% Input Parameters 
%% Change these as desired - these govern the majority of the results
Lake=6;                  %% Select lake chemistry (1-Seawater, 2-MgSO4 Eruopa Ocean, 3-Salt Lake, 4-Basque Lake 2.2,
                         %  5-Basque Lake 2.3, 6-Last Chance Lake
                         
S_bottom=70;            %% Ocean Salinity (ppt)
S_sat0=70;               %% Initial ocean Saturation (ppt) NaCl - 357
Sat_slope=5;          %% Slope of temperature dependent saturation (ppt/K)
Ttop=255;                %% Ice top temp (K) (use T_surf below if you want evolving surface temp)

[Brine_T,Brine_rho]=Ocean_Param_Calc(S_bottom,Lake); %% Calculates ocean temp (0.01K above liquidus) and density from salnity
                                              %  and ocean composition

Tbottom=Brine_T;         %% Ocean Temp (K) [MgSO4 12.3ppt-273.0K, 100ppt-271.31K, 282ppt-267.5K NaCl 34ppt-271.2K
                         %                  NaCl 100ppt - 267.78, NaCl 200ppt - 261.08] *always put above liquidus*
Tm=273.15;               %% Melt temperature of pure ice (K)
k_i=2;                   %% Thermal conductivity (ice) W/m*K
k_br=.6;                 %% '                            ' (brine)
k_s=2*10^-9;             %% Diffusivity for Salt
dt=100;                   %% Time step (s)
tf=5000000;              %% Final time (s)
dz=.001;                 %% Spatial discretization (m)
H=0.45;                  %% Domain size (m)
c_br=3985;               %% Specific heat of seawater (J/K)
c_i=2000;                %% '              ' ice
L=334774;                %% Latent heat of fusion ice<->water (J/Kg)
L_salt=14265734;         %% Latent heat of formation for precipitated salts (epsomite - 13744870) (J/Kg)
                         %  (Natron - 14265734)
rho_i=917;               %% Density of Ice (Kg/m^3)
rho_br=Brine_rho;        %% Density of Ocean (used in volume averaging - 1D grav. drainage uses delta S) 34ppt NaCl-1027
                         %  100ppt NaCl-1081 200ppt NaCl-1168 12.3ppt MgSO4-1012, 100pppt-1103, 282ppt-1323
STol=0.1;                %% Error Tolerance for Salinity
PhiTol=0.005;            %% Error Tolerance for Porosity
TTol=0.01;               %% Error Tolerance for Temperature
m=2;                     %% Cementation exponent
g=9.8;                   %% Earth Gravity
mu=1.88*10^-3;           %% Viscosity of Ocean
phi_c=0.05;              %% Critical Porosity

%% Sinusoidal surface temp, BE SURE TO CHANGE 'Ttop' in the calling of 'one_D_adv_ARD_interactive_lakes.m'
%  to 'T_surf(n)' if using this

T_surf=273.15-15*sin(2*pi*[1:dt:31557600]/31557600);  %% Sinusoidal surface temp


%% Bottom Boundary Condition options Open=1 or Closed=2 bottom boundary (Salt)
%BCoc=1;
BCoc=2;
S_0=H*S_bottom;
%% Dirichlet=1 or Neumann=2 bottom boundary condition (Temp)
%BCdn=1;
BCdn=2;

dT=0.025;      %% Temperature gradient at bottom boundary (K/m) - geothermal flux

%% Plotting Vector
Depth=[0:-dz:-H];

%% Initial Condition Vectors for Phi and T  and S
    T_initial=0*[0:dz:H]+Tbottom;
    S_initial=0*[0:dz:H]+S_bottom;
    phi_initial=0*[0:dz:H]+1;

%% Vectors for function iterations
Temperature=[];
Liquid_Fraction=[];
Salinity=[];

T_hold=[];
Phi_hold=[];
Sal_hold=[];

Salt_Precip=0;
Salt_Layer=0;

%% Looping Over Time
for n=1:tf/dt
    if n==1;

        [NewTemp,NewPhi,NewS]=one_D_adv_ARD_interactive_lakes(...
        k_i,k_br,k_s,dt,dz,H,c_i,c_br,L,STol,TTol,PhiTol,phi_initial,...
        T_initial,S_initial,S_bottom,rho_i,...
        rho_br,Tm,m,T_surf(n),Tbottom,g,rho_br,phi_c,BCoc,BCdn,dT,Lake);

        Temperature=NewTemp;
        Liquid_Fraction=NewPhi;
        Salinity=NewS;
        
        % Marking the ice-ocean interface
        IO_int=0;
        if BCoc==2
            for i=1:length(NewPhi)
                if NewPhi(i)<0.99
                    IO_int=IO_int+1;
                else
                end
            end
            
            % Calculating salinity of remaining brine due to gravity
            % drainage
            S_ice=dz*sum(Liquid_Fraction(1:IO_int).*Salinity(1:IO_int));
            S_br=S_0-S_ice;
            S_ppt=S_br/((length(NewPhi)-IO_int)*dz);
            Salinity(IO_int+1:end)=S_ppt;
            
            % Salt precipitation and accompanying heat of formation
            S_sat=S_sat0+Sat_slope*(Temperature(end)-Tbottom);
            if S_ppt>S_sat
                S_hold=(S_ppt-S_sat)*(length(NewPhi)-IO_int);
                Salt_Precip=Salt_Precip+S_hold;
                Salt_Layer=floor(Salt_Precip/1000);
                Salinity(IO_int+1:end)=S_sat;
                % Temperature mod
                del_T=((L_salt/c_br)*(S_hold/1000))/(length(NewPhi)-IO_int);
                Temperature(IO_int+1:end)=Temperature(IO_int+1:end)+del_T;
                if Salt_Layer==length(NewPhi)-IO_int
                    break
                else
                end
            else
            end
        else
        end

    else
        [NewTemp,NewPhi,NewS]=one_D_adv_ARD_interactive_lakes(...
        k_i,k_br,k_s,dt,dz,H,c_i,c_br,L,STol,TTol,PhiTol,Liquid_Fraction',...
        Temperature',Salinity',S_ppt,rho_i,...
        rho_br,Tm,m,T_surf(n),Tbottom,g,rho_br,phi_c,BCoc,BCdn,dT,Lake);
        
        
        Temperature=NewTemp;
        Liquid_Fraction=NewPhi;
        Salinity=NewS;
        
        if ismember(n*dt,[0:100000:10000000])==1
            T_hold=[T_hold Temperature];
            Phi_hold=[Phi_hold Liquid_Fraction];
            Sal_hold=[Sal_hold Salinity];
        else
        end
        
        IO_int=0;
        if BCoc==2
            for i=1:length(NewPhi)
                if NewPhi(i)<0.99
                    IO_int=IO_int+1;
                else
                end
            end
            
            % Calculating salinity of remaining brine due to gravity
            % drainage
            S_ice=dz*sum(Liquid_Fraction(1:IO_int).*Salinity(1:IO_int));
            S_br=S_0-S_ice;
            S_ppt=S_br/((length(NewPhi)-IO_int)*dz);
            Salinity(IO_int+1:end)=S_ppt;
            
            % Salt precipitation and accompanying heat of formation
            S_sat=S_sat0+Sat_slope*(Temperature(end)-Tbottom);
            if S_ppt>S_sat
                S_hold=(S_ppt-S_sat)*(length(NewPhi)-IO_int);
                Salt_Precip=Salt_Precip+S_hold;
                Salt_Layer=floor(Salt_Precip/1000);
                Salinity(IO_int+1:end)=S_sat;
                % Temperature mod
                del_T=((L_salt/c_br)*(S_hold/1000))/(length(NewPhi)-IO_int);
                Temperature(IO_int+1:end)=Temperature(IO_int+1:end)+del_T;
                if Salt_Layer==length(NewPhi)-IO_int
                    break
                else
                end
            else
            end
        else
        end
%         if BCoc==2
%             for i=1:length(NewPhi)
%                 if NewPhi(i)<1
%                     IO_int=IO_int+1;
%                 else
%                 end
%             end
%             S_ice=dz*sum(Liquid_Fraction(1:IO_int).*Salinity(1:IO_int));
%             S_br=S_0-S_ice;
%             S_ppt=S_br/((length(NewPhi)-IO_int)*dz);
%             Salinity(IO_int+1:end)=S_ppt;
%         else
%         end
        
    end
end

% for i=1:length(Liquid_Fraction)
%     if Liquid_Fraction(i)<=phi_c
%         Porosity(i)=phi_c;
%     else
%         Porosity(i)=Liquid_Fraction(i);
%     end
% end

%% Bulk salinity vs depth
figure
subplot(1,4,1)
plot(T_hold,Depth)
title('Temperature vs. Depth')
xlabel('Temperature (K)')
ylabel('Depth (m)')
subplot(1,4,2)
plot(Phi_hold,Depth)
title('Porosity vs. Depth')
xlabel('Porosity')
ylabel('Depth (m)')
subplot(1,4,3)
plot(Sal_hold,Depth)
title('Salinity vs. Depth')
xlabel('Salinity (ppt)')
ylabel('Depth (m)')
xlim([0 max(Salinity)+5])
% subplot(1,5,4)
% plot(SinI,Depth)
% title('Ice Salinity vs. Depth')
% xlabel('Salinity (ppt)')
% ylabel('Depth (m)')
subplot(1,4,4)
plot(Phi_hold.*Sal_hold,Depth)
title('Bulk Salinity vs. Depth')
xlabel('Bulk Salinity (ppt)')
ylabel('Depth (m)')
xlim([0 max(Salinity.*Liquid_Fraction)+5])

figure
subplot(1,4,1)
plot(Temperature,Depth)
title('Temperature vs. Depth')
xlabel('Temperature (K)')
ylabel('Depth (m)')
subplot(1,4,2)
plot(Liquid_Fraction,Depth)
title('Porosity vs. Depth')
xlabel('Porosity')
ylabel('Depth (m)')
subplot(1,4,3)
plot(Salinity,Depth)
title('Salinity vs. Depth')
xlabel('Salinity (ppt)')
ylabel('Depth (m)')
xlim([0 max(Salinity)+5])
% subplot(1,5,4)
% plot(SinI,Depth)
% title('Ice Salinity vs. Depth')
% xlabel('Salinity (ppt)')
% ylabel('Depth (m)')
subplot(1,4,4)
plot(Liquid_Fraction.*Salinity,Depth)
title('Bulk Salinity vs. Depth')
xlabel('Bulk Salinity (ppt)')
ylabel('Depth (m)')
xlim([0 max(Salinity.*Liquid_Fraction)+5])

% velocity=[];
% for i=2:length(vel_hold)
%     velocity(i-1)=(vel_hold(i,1)-vel_hold(i-1,1))/(vel_hold(i,2)-vel_hold(i-1,2));
% end
% 
% figure
% plot(vel_hold(2:end,2),velocity)
% title('Growth Velocity vs. Time')
% xlabel('Growth Velocity (m/s)')
% ylabel('Bulk Salinity (ppt)')

% k_eff=(1/S_bottom)*salt_hold(2:end);
% log_k_eff=log((1./k_eff)-1);
% figure
% plot(100*velocity,log_k_eff)
% title('Growth Velocity vs. Log(Keff)')
% xlabel('Growth Velocity (cm/s)')
% ylabel('ln((1/Keff)-1)')

% figure
% plot([1:31557600],-15*sin(2*pi*[1:31557600]/31557600)+5*sin(2*pi*[1:31557600]/86400))