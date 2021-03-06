function P2=liquidus(T_HIGH, T_LOW, C_HIGH, C_LOW, INC, SALINITY, Na, K, Ca,...
    Mg, Cl, Sulfate, Nitrate)
%sw=[0 1 1 0.469 0.0102 0.0103 0.0528 0.546 0.0282 0.0 0.0 0.0 273.15 272.15 1.0];
%europa=[0 1 1 0.0491 0.001964 0.009637 0.06271 0.02087 0.08744 0.0 0.0 0.0 273.15 272.15 1.0];
precip=[];
temps=[];
conc=[];
T=272.15:-1:261.15;
%Thigh=273.15;
Thigh=T_HIGH;
%Tlow=243.15;
Tlow=T_LOW;
%chigh=9;
chigh=C_HIGH;
%clow=1;
clow=C_LOW;
%inc=2.00;
inc=INC;
%salinity=35;
salinity=SALINITY;
%for j=1:50
%    temps(j,:)=T;
%end
%for j=1:50
%    conc(j,:)=1:10;
%end
hydr=[];
mirab=[];
caSO4=[];
ice=[];
hydr2=[];
mirab2=[];
caSO42=[];
ice2=zeros((chigh-clow)/inc+1,int8((Thigh-Tlow)/inc+1));
T2=zeros((chigh-clow)/inc+1,int8((Thigh-Tlow)/inc+1));
conc2=zeros((chigh-clow)/inc+1,int8((Thigh-Tlow)/inc+1));
k=0;
for j=Thigh:-inc:Tlow
    k=k+1;
    m=0;
    for i=clow:inc:chigh
        m=m+1;
    %seawater
%     [conditions,solution,solid]=FREZCHEM([0 1 1 i*0.469 i*0.0102 i*0.0103 i*0.0528 i*0.546...
%         i*0.0282 0.0 0.0 0.0 274.15 j 274.15-j]);
    [conditions,solution,solid]=FREZCHEM([0 1 1 i*Na i*K i*Ca i*Mg i*Cl...
        i*Sulfate i*Nitrate 0.0 0.0 274.15 j 274.15-j]);
    %Europa
    %[conditions,solution,solid]=FREZCHEM([0 1 1 i*0.0491 i*0.001964 i*0.009637 i*0.06271 i*0.02087...
    %    i*0.08744 0.0 0.0 0.0 274.15 j 274.15-j]);
    hydr=[hydr;i,j,solid(2,1)];
    mirab=[mirab;i,j,solid(11,1)];
    caSO4=[caSO4;i,j,solid(18,1)];
    ice=[ice;i,j,solid(1,1)];
    ice2(m,k)=solid(1,1);
    T2(m,k)=j;
    conc2(m,k)=i;
    
    hydr2=[hydr2;(sum(solution(:,2))/sum(solution(:,1))),j,solid(2,1)];
    mirab2=[mirab2;(sum(solution(:,2))/sum(solution(:,1))),j,solid(11,1)];
    caSO42=[caSO42;(sum(solution(:,2))/sum(solution(:,1))),j,solid(18,1)];
    end
end

count=[];
Tj=[];
Tjm1=[];
concj=[];
concjm1=[];
for i=1:(chigh-clow)/inc+1;
    for j=1:(Thigh-Tlow)/inc+1;
        if ice2(i,j)>0
            Tj=[Tj T2(i,j)];
            Tjm1=[Tjm1 T2(i,j-1)];
            concj=[concj conc2(i,j)];
            concjm1=[concjm1 conc2(i,j-1)];
            break
        else
        end
    end
end


%precip=precip(2:end,:);
% figure
% subplot(2,3,1)
% scatter3(hydr(:,1),hydr(:,2),hydr(:,3),3,hydr(:,3))
% xlabel('Concentration (xSeawater)')
% ylabel('Temperature (K)')
% zlabel('Precipitate (mol)')
% title('Hydrohalite Precipitation')
% subplot(2,3,2)
% scatter3(mirab(:,1),mirab(:,2),mirab(:,3),3,mirab(:,3))
% xlabel('Concentration (xSeawater)')
% ylabel('Temperature (K)')
% zlabel('Precipitate (mol)')
% title('Mirabilite Precipitation')
% subplot(2,3,3)
% scatter3(caSO4(:,1),caSO4(:,2),caSO4(:,3),3,caSO4(:,3))
% xlabel('Concentration (xSeawater)')
% ylabel('Temperature (K)')
% zlabel('Precipitate (mol)')
% title('CaSO4 Hydrate Precipitation')
% subplot(2,3,4)
% scatter3(hydr2(:,1),hydr2(:,2),hydr2(:,3),3,hydr2(:,3))
% xlabel('Concentration (xSeawater)')
% ylabel('Temperature (K)')
% zlabel('Precipitate (mol)')
% title('Hydrohalite Precipitation')
% subplot(2,3,5)
% scatter3(mirab2(:,1),mirab2(:,2),mirab2(:,3),3,mirab2(:,3))
% xlabel('Concentration (xSeawater)')
% ylabel('Temperature (K)')
% zlabel('Precipitate (mol)')
% title('Mirabilite Precipitation')
% subplot(2,3,6)
% scatter3(caSO42(:,1),caSO42(:,2),caSO42(:,3),3,caSO42(:,3))
% xlabel('Concentration (xSeawater)')
% ylabel('Temperature (K)')
% zlabel('Precipitate (mol)')
% title('CaSO4 Hydrate Precipitation')
% 
% figure
% scatter3(ice(:,1),ice(:,2),ice(:,3),3,ice(:,3))
% xlabel('Concentration (xSeawater)')
% ylabel('Temperature (K)')
% zlabel('Precipitate (mol)')
% title('Ice (mol)')

P=polyfit(salinity*(concj+concjm1)/2,(Tj+Tjm1)/2,1);
P2=polyfit(salinity*(concj+concjm1)/2,(Tj+Tjm1)/2,2);
fitline=P(1)*salinity*[clow:inc:chigh]+P(2);
fitline2=P2(1)*(salinity*[clow:inc:chigh]).^2+P2(2)*salinity*[clow:inc:chigh]+P2(3);
% figure
% plot(salinity*concj,Tj,salinity*concjm1,Tjm1,salinity*[clow:inc:chigh],fitline,salinity*[clow:inc:chigh],fitline2)
% xlabel('Concentration (ppt)')
% ylabel('Temperature (K)')
% quad_fit=fittype('(a*S^2)+(b*S)+c','dependent',{'Tmelt'},'independent',{'S'},'coefficients',{'a','b','c'});
% quad_line=fit(salinity*(concj'+concjm1')/2,(Tj'+Tjm1')/2,quad_fit);

figure
scatter3(salinity*ice(:,1),ice(:,2),ice(:,3),3,ice(:,3))
hold on
plot3(salinity*[clow:inc/10:chigh],...
    P2(1)*(salinity*[clow:inc/10:chigh]).^2+P2(2)*salinity*[clow:inc/10:chigh]+P2(3),...
    0*[clow:inc/10:chigh],'k')
scatter3(salinity*[clow:inc:chigh],fitline2,0*[clow:inc:chigh],30,'red','filled')
scatter3(salinity*[clow:inc/10:chigh],...
    P2(1)*(salinity*[clow:inc/10:chigh]).^2+P2(2)*salinity*[clow:inc/10:chigh]+P2(3),...
    0*[clow:inc/10:chigh],10,'red','filled')
surf(salinity*[clow:inc:chigh],Thigh:-inc:Tlow,transpose(ice2),'FaceAlpha',0.5)
text(max(salinity*ice(:,1)),max(ice(:,2)),max(ice(:,3))-(0.1*(max(ice(:,3))-min(ice(:,3)))),...
    strcat('Tmelt=(',num2str(P2(1)),')S^2+(',num2str(P2(2)),')S+',num2str(P2(3))),...
    'FontSize',12)
xlabel('Concentration (ppt)')
ylabel('Temperature (K)')
zlabel('Precipitated Ice (mol)')
title('Ice (mol)')

save('Coefficients.mat','P2')
