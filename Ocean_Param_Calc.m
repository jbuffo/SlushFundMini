%% function calculating ocean temperature and density using liquidus curves and given salinity
function [T,rho]=Ocean_Param_Calc(S,Lake)

if Lake==1
        T=(-(9.1969758*(10^-5)*S^2)-0.03942059*S+272.63617665)+0.01;   %% Liquidus point using FREZCHEM seawater
        rho=1000+0.81*S;
    elseif Lake==2
        T=(-(1.333489497*(10^-5)*S^2)-0.01612951864*S+273.055175687)+0.01; %% Liquidus point using FREZCHEM Europa ocean
        rho=1000+S;
    elseif Lake==3
        T=(-(2.1302*(10^-5)*S^2)-0.01307*S+273.0403)+0.01; %% Liquidus point using FREZCHEM Salt Lake
        rho=1000+S;
    elseif Lake==4
        T=(-(2.9486*(10^-5)*S^2)-0.012162*S+273.00)+0.01;   %% Liquidus point using FREZCHEM Basque Lake 2.2 
        rho=1000+S;
    elseif Lake==5
        T=(-(3.2684*(10^-5)*S^2)-0.010903*S+273.0587)+0.01; %% Liquidus point using FREZCHEM Basque Lake 2.3
        rho=1000+S;
    else
        T=(-(1.0579*(10^-5)*S^2)-0.033321*S+272.8703+0.01); %% Liquidus point using FREZCHEM Last Chance Lake 
        rho=1000+0.71*S;
end


    
% SL_liq=((-2.1302*10^-5)*salt_sweep.^2)-0.01307*salt_sweep+273.0403;
% BL2p2_liq=((-2.9486*10^-5)*salt_sweep.^2)-0.012162*salt_sweep+273;
% BL2p3_liq=((-3.2684*10^-5)*salt_sweep.^2)-0.010903*salt_sweep+273.0587;
% LCL_liq_all_Na=((-1.0579*10^-5)*salt_sweep.^2)-0.033321*salt_sweep+272.8703;

