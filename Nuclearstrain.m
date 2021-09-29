%Hongyuan Zhu
%Xi'an Jiaotong University
%9 September 2021
% The code is calculating the change of nuclear strain when cells are transferred from TCP to soft gels (SR/HR)

F0=150; %pN, traction force on TCP, obtained from motor-clutch model
F1=125; %pN, traction force on SR, obtained from motor-clutch model
F2=56; %pN, traction force on HR, obtained from motor-clutch model

S=1e4; %nm2, the nuclear area drew by actin filament

EN=102; %pN, effective elasticity of nucleus
etaN=4500; %pN*day, effetive viscocity of nucleus

T=10;%day
dt=0.01; %day
N=(T+1)/dt;
i=0;
epsilonN1=zeros(N+1,1); %Nuclear strain TCP-SR
epsilonN2=zeros(N+1,1); %Nuclear strain TCP-HR
Time=zeros(N+1,1);
for t=0:dt:T %day
    i=i+1;
    epsilonN1(i,1)=(1/EN+t/etaN)*F0;
    epsilonN2(i,1)=epsilonN1(i,1);
    Time(i,1)=t;
end

for t=0:dt:3
    i=i+1;
    epsilonN1(i,1)=T/etaN*F0+(1/EN+t/etaN)*F1;
    epsilonN2(i,1)=T/etaN*F0+(1/EN+t/etaN)*F2;
    Time(i,1)=t+T;
end
