%THE ONE TO USE (UPDATE) TU VALIDATE.only change line 196 and 197 (S and
%Tair)
clc
clear all
close all

%Given operational data
%PS255 ECSOLAR Polycristallin
%
C = 1;
Voc_r=37.6 * C;  %[V]
Isc_r=8.79 * C;  %[A]
Vmp_r=30.4 * C;  %[V]
Imp_r=8.39 * C;  %[A]
Sref=1000;   %[W/m2]
beta= -0.33*Voc_r/100;         % [1/K] voltage temp coefficient
alpha= 0.055*Isc_r/100;        % [1/K] current T coefficient
Tnoct=43;                      % [°C] 45+/-2
 


%  Voc_r=37.1;  %[V]
% Isc_r=8.32;  %[A]
% Vmp_r=29.2;  %[V]
% Imp_r=7.88;  %[A]
% Sref=1000;   %[W/m2]
% beta= -0.33*Voc_r/100;         % [1/K] voltage temp coefficient
% alpha= 0.055*Isc_r/100; 

%BP365
% Voc_r=22.1;  %[V]
% Isc_r=3.99;  %[A]
% Vmp_r=17.6;  %[V]
% Imp_r=3.69;  %[A]
% Sref=1000;   %[W/m2]
% beta= -0.08;        % [1/K] voltage temp coefficient
% alpha= 0.065*Isc_r;        % [1/K] current T coefficient


Tref=25;     %[°C]
Tref=25+273.15;  %[K]
Eg_r=1.125;   %[eV] material band gap at 25°C for polycristalline silicon


Ns=60 * C; %number cell series
Np=1; %number cell parallel
%number of cell series, series to increase V while parallel to increase I.


q=1.602e-19;    % [C] electronic charge
k=1.3806503e-23;     % [J/K] Boltzmann's constant  r

% Some constant are introduced to avoid cluttering

c1=k*Tref/q;
c2=Voc_r/(Ns*c1);
c3=Vmp_r/(Ns*c1);
c4=Imp_r/(Np*c1);
c5=Isc_r/(Np*c1);

% Temperature and Irradiance, different conditions

T=24+273.15;
S=1000;

%Bandgap Energy (eV)
Eg_r=1.17-0.000473*((Tref)^2/(Tref+636));
Eg=1.17-0.000473*(T^2/(T+636));

%Temperature-dependent factor for shunt conductance
KT=(T/Tref)^3*exp(Eg_r*q/(k*Tref)-Eg*q/(k*T));


%Solving non linear system

%Rmax, n_r=2
Rs_max=(2/c4*(1+lambertw(-exp((c2-2-2*c3)/2)))+c3/c4);

%Initialize
x0 = [Rs_max/2 1];
lb=  [0 0.5];
ub=  [Rs_max 2];

foff = @(x) TwoParameters(x,Ns,k,Tref,q, beta, T, alpha,c1,c2,c3,c4,c5, KT)

options = optimset('TolFun', 1e-18, 'Display','iter');


%y = fsolve(foff, x0, options);
y = lsqnonlin(foff,x0,lb,ub, options);
%y=fmincon(foff,x0,[],[],[],[],lb,ub,[], options);


[res,Rs_r,n_r]= TwoParameters(y,Ns,k,Tref,q, beta, T, alpha,c1,c2,c3,c4,c5, KT);


if Rs_r==(c2-c3)/c4;
    
    error('Rs not acceptable')
end

%Il_r,Iz_r and Rp_r are now calculated

Il_r = c1*c2*c4 / (c3-c4*Rs_r) + ((c1*c4*(2*c3-c2)) / (c3-c4*Rs_r) * ...
    (exp(c2/n_r)-1-c2/n_r * exp((c3+c4*Rs_r)/n_r)))/(exp(c2/n_r)+...
    ((c4*Rs_r+c3-c2)/n_r-1)*exp((c3+c4*Rs_r)/n_r));

Iz_r = c1*c4 / (c3-c4*Rs_r) * (2*c3-c2) / ((exp(c2/n_r)+...
    ((c4*Rs_r+c3-c2)/n_r-1)*exp((c3+c4*Rs_r)/n_r)));

%Gp_r = (c4*(((1+(c3-c4*Rs_r)/n_r))*exp((c3+c4*Rs_r-c2)/n_r)-1))/((1+((c4*Rs_r+c3-c2)/n_r-1)*exp((c3+c4*Rs_r-c2)/n_r))*(c4*Rs_r-c3));

Gp_r = abs((Il_r-Iz_r*(exp(c5*Rs_r/n_r)-1)-c1*c5)/(c1*c5*Rs_r));

Rp_r=1/Gp_r;

%I-V curve at STC
 a_ref=k*Tref*Ns*n_r/q;
 
 V_r=linspace(0,Voc_r,100);
 Idc_r=zeros(1,100);   
    index = 1;
    
    for V_r=linspace(0,Voc_r,100)
        Idc_r(index) = fsolve( @(I_r) -I_r + Np*Il_r-Np*Iz_r*(exp((V_r+I_r*Ns/Np*Rs_r)/a_ref)-1)-(V_r+I_r*Rs_r*Ns/Np)/(Ns/Np*Rp_r), 0);
        index = index + 1;
    end
  
  V_r=linspace(0,Voc_r,100);
% 
% I_r=Np*Il_r-Np*Iz_r*(exp((V_r+I_r*Ns/Np*Rs_r)/a_ref)-1)-(V_r+I_r*Rs_r*Ns/Np)/(Ns/Np*Rp_r);
% 
% figure(1)
% plot(V_r,Idc_r)

%I-V curve different T
%Rs_r=0.02;
deltaT = [5 8 10 12];
IDC = zeros(100,length(deltaT));
T=Tref;

S=780*ones(1,length(deltaT));

for i=1:length(deltaT)
    
    T=deltaT(i)+273.15;

    
   % T=Tref;
   
    
    a(i)=a_ref*T/Tref;
    
    Eg(i)=1.17-0.000473*(T^2/(T+636));
    
    Iz(i) = Iz_r*(T/Tref)^3*exp(q/(n_r*k)*(Eg_r/Tref-Eg(i)/T));
    
    Il(i)=S(i)/Sref*(Il_r + alpha*(T-Tref));
    
    Rp=Sref/S(i)*Rp_r;
    
    Rs=Rs_r;
    
    Voc(i)= Voc_r + beta*(T-Tref);
    
    
%     I=zeros(1,100);
%     V=linspace(0,Voc,100);
%     I = Np*Il-Np*Iz*(exp((V+I*Ns/Np*Rs)/a)-1)-(V+I*Rs*Ns/Np)/(Ns/Np*Rp);
    

    V(:,i)=linspace(0,Voc(i),100);
    [row,~]=size(V);
    
    for kk = 1:row
        IDC(kk,i) = fsolve( @(I) -I + Np*Il(i)-Np*Iz(i)*(exp((V(kk,i)+I*Ns/Np*Rs)/a(i))-1)-(V(kk,i)+I*Rs*Ns/Np)/(Ns/Np*Rp), 0);
        P(kk,i)=IDC(kk,i)*V(kk,i);        
    end
    
end
figure;
for i = 1 : size(V,2)
  plot(V(:,i),IDC(:,i))
  hold on;
end
figure;
for i = 1 : size(V,2)
  plot(V(:,i),P(:,i))
  hold on;
  PM(i) = max(P(:,i));
end
% S_star=repmat(S,100,1);
% figure(2)
% surf(S_star,V,IDC);%,S_star);
% ylabel('voltage')
% zlabel('current')
% xlabel('Irradiance')
% %legend('25°C','35°C','45°C','55°C')
% axis([0 inf 0 inf])
% 
% figure(3)
% surf(S_star,IDC,V);%,S_star);
% zlabel('voltage')
% ylabel('current')
% xlabel('Irradiance')
% %legend('25°C','35°C','45°C','55°C')
% axis([0 inf 0 inf])
% 
% figure(4)
% plot(S_star,V)
% %legend('25°C','35°C','45°C','55°C')
% axis([0 inf 0 inf])

nop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Maximum power from the module

Tair=23
S=155;
Tc=Tair+(Tnoct-20)*S/800;

T=Tc+273.15;
a=a_ref*T/Tref;

Eg=1.17-0.000473*(T^2/(T+636));

Iz = Iz_r*(T/Tref)^3*exp(q/(n_r*k)*(Eg_r/Tref-Eg/T));


Il=S/Sref*(Il_r + alpha*(T-Tref));

Rp=Sref/S*Rp_r;

Rs=Rs_r;

Voc= Voc_r+beta*(T-Tref);
Isc=Isc_r+alpha*(T-Tref);

%Iz=(Isc-(Voc-Isc*Rs)/Rp)*exp(-Voc/a);

V=linspace(0,Voc,100);
IDC = zeros(1,100);
index = 1;

for V=linspace(0,Voc,100)
    IDC(index) = fsolve( @(I) -I + Np*Il-Np*Iz*(exp((V+I*Ns/Np*Rs)/a)-1)-(V+I*Rs*Ns/Np)/(Ns/Np*Rp), 0);
    P(index)=IDC(index)*V;
     index = index + 1;
    
end
V=linspace(0,Voc,100);
res=-IDC + Np*Il-Np*Iz*(exp((V+IDC*Ns/Np*Rs)/a)-1)-(V+IDC*Rs*Ns/Np)/(Ns/Np*Rp);

%Maximum power from converter

  max=find(P==max(P));
  
  Vmp=V(max)
  Imp=IDC(max)
  
  Pdc1=(Vmp*14)*(2*Imp)

  Pdc2=(Vmp*11)*Imp;
  
  
  Pdc=Pdc1+Pdc2;
  PDC=39*P(max)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

