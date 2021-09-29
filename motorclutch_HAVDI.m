%Hongyuan Zhu
%Xi'an Jiaotong University
%9 September 2021

%Reduplicate the Chan & Bangasser Motor-Clutch model 
%Including the effect of HAVDI on kon-stiffness relationship into the model
% 0 - disengaged clutch
% 1 - engaged clutch

clear
clc

tic

nm=75; %number of myosin motors
Fm=-2; %motor stall force in pN
vu=-110; %unloaded motor velocity in nm/s
nc=75; %number of molecular clutches
% kon=0.3; %On rate constant in 1/s
koff=0.1; %Off rate constant in 1/s
Fb=-6; %Bond rupture force in pN

% koff1=1e-3; %Off rate constant in 1/s
% koff2=2.5; %Off rate constant in 1/s
% Fb=-9.4; %Bond rupture force in pN
% Fb2=-10.3; % fitting constant in pN

Kc=0.8; %Clutch spring constant in pN/nm
gain=0; %gain of feedback loop
ra=550; %radius of adhesion in nm

events=1e5; %number of events to simulate
stiffness=logspace(0,6,301); %define stiffness vector in kPa
% stiffness=1e6;
retro=zeros(length(stiffness),1); %initialize retro flow vector
failfreq=zeros(length(stiffness),1); %initialize failure frequency vector
avnumcon=zeros(length(stiffness),1); %initialize average number of engaged clutches vector
avcough=zeros(length(stiffness),1); %initialize average koff vector
avtime=zeros(length(stiffness),1); %initialize average timestep vector
avtrac=zeros(length(stiffness),1); %initialize average traction force vector
avsubpos=zeros(length(stiffness),1); %initialize average substrate position vector
avFc=zeros(length(stiffness),1); %initialize average clutch force vector
binds=zeros(length(stiffness),1); %initialize number of binding events
failures=zeros(length(stiffness),1); %initialize number of failure events
cyct=zeros(length(stiffness),1); %initialize cycle time
velcv=zeros(length(stiffness),1); %initialize velocity cv
rebcyc=zeros(length(stiffness),1); %initialize rebinds per cycle
kon=zeros(length(stiffness),1); 
for jj=1:length(stiffness)
    
    Ksub=4*pi*ra/9/1000*stiffness(jj) %Substrate spring constant in pN/nm

       kon(jj)=0.57*exp(-72/(40+stiffness(jj)));%SR
%     kon(jj)=0.42*exp(-72/(40+stiffness(jj)));%HR
  
    %initialize vectors for plots
    cstate=zeros(1,nc); %clutch state vector
    cunbind=zeros(1,nc); %clutch unbind state vector
    crebind=zeros(1,nc); %clutch rebind state vector
    xc=zeros(1,nc); %clutch position vector
    t=zeros(1,events+1); %time vector
    subpos=zeros(1,events+1); %substrate position vector
    numcon=zeros(1,events+1); %number of engaged clutches vector
    numcoff=zeros(1,events+1); %number of disengaged clutches vector
    vel=zeros(1,events+1); %velocity vector
    timestep=zeros(1,events+1); %vector of dt's
    cough=zeros(1,events+1); %koff vector
    Ft=zeros(1,events+1); %traction force vector
    totFc=zeros(1,events+1); %mean engaged clutch tension
    
    % Set inital state
    i=1;
    ceng=find(cstate==1); %find indices of engaged clutches
    cdisen=find(cstate==0); %find indices of disengaged clutches
    vf=vu; %claculate actin filament velocity
    xc(ceng)=xc(ceng)+vf*0.005; %claculate positions of engaged clutches(dt=0.005)
    xsub=(Kc*sum(xc(ceng))/(Ksub+length(ceng)*Kc)); %calculate substrate posisiton
    %vf=vu*(1-((Ksub*xsub)/(nm*Fm))); %claculate actin filament velocity
    %xc(ceng)=xc(ceng)+vf*0.005; %claculate positions of engaged clutches (dt=0.005)
    xc(cdisen)=xsub; %calculate posisiton of disengaged clutches
    Fc=Kc*(xc-xsub); %calculate force on each clutch
    t(i)=0;
    subpos(i)=xsub;
    numcon(i)=length(ceng);
    numcoff(i)=length(cdisen);
    vel(i)=-vf;
    timestep(i)=0;
    
    %Event stepping
    while i<=events
        i=i+1;
        
        %calculate clutch binding times
        if isempty(cdisen)
            tbind=inf;
        else
            tbind=-log(rand(1,length(cdisen)))/kon(jj);
        end
        
        %calculate clutch unbinding times
        if isempty(ceng)
            tunbind=inf;
%             cough(i)=koff1+koff2;
            cough(i)=koff;
            totFc(i)=0;
        else
            tunbind=-log(rand(1,length(ceng)))./(koff*exp(Fc(ceng)./(Fb+gain*Fc(ceng))));
            cough(i)=mean(koff*exp(Fc(ceng)./(Fb+gain*Fc(ceng))));
%             tunbind=-log(rand(1,length(ceng)))./(koff1*exp(Fc(ceng)./Fb)+koff2*exp(-Fc(ceng)./Fb2));
%             cough(i)=mean(koff1*exp(Fc(ceng)./Fb)+koff2*exp(-Fc(ceng)./Fb2));
            totFc(i)=mean(Fc(ceng));
        end
        
        %find minimum time and execute that event
        [dtbind indbind]=min(tbind);
        [dtunbind indunbind]=min(tunbind);
        if dtbind<dtunbind %free clutch bind to actin
            cstate(cdisen(indbind))=1;
            dt=dtbind;
            binds(jj)=binds(jj)+1;
            if cunbind(cdisen(indbind))==1 %if clutch has already unbound during the cycle
                crebind(cdisen(indbind))=crebind(cdisen(indbind))+1;
            end
        else %engaged clutch disengages from actin
            cstate(ceng(indunbind))=0;
            dt=dtunbind;
            cunbind(ceng(indunbind))=1;
        end
        
        ceng=find(cstate==1); %find indices of engaged clutches
        cdisen=find(cstate==0); %find indices of disengaged clutches
        Ftrac=Ksub*xsub; %claculate traction force
        vf=vu*(1-((Ksub*xsub)/(nm*Fm))); %claculate actin filament velocity
        %vf=vu;
        xc(ceng)=xc(ceng)+vf*dt; %claculate positions of engaged clutches
        xsub=(Kc*sum(xc(ceng))/(Ksub+length(ceng)*Kc)); %calculate substrate posisiton
        %Ftrac=Ksub*xsub; %claculate traction force - old
        %vf=vu*(1-((Ksub*xsub)/(nm*Fm))); %claculate actin filament - old
        %xc(ceng)=xc(ceng)+vf*dt; %claculate positions of engaged clutches - old
        xc(cdisen)=xsub; %calculate posisiton of disengaged clutches
        Fc=Kc*(xc-xsub); %calculate force on each clutch
        
        if xsub==0 %reset unbind vector at failure event
            cunbind=zeros(1,nc);
        end
        
        t(i)=t(i-1)+dt;
        subpos(i)=xsub;
        numcon(i)=length(ceng);
        numcoff(i)=length(cdisen);
        vel(i)=-vf;
        timestep(i)=dt;
        Ft(i)=Ftrac;
    end
    cyctime=diff(t(subpos==0)); %cycle time
    
    retro(jj,1)=sum((vel.*timestep)/t(events+1)); %weighted average retrograde flowrate
    failfreq(jj,1)=(length(numcon(numcon==0))-1)/(t(events+1)); %failures/s
    avnumcon(jj,1)=sum((numcon.*timestep)/t(events+1)); %average number of engaged clutches
    avcough(jj,1)=sum((cough.*timestep)/t(events+1)); %average koff
    avtime(jj,1)=mean(timestep); %average timestep
    avtrac(jj,1)=abs(sum((Ft.*timestep)/t(events+1))); %average traction force
    avsubpos(jj,1)=sum((subpos.*timestep)/t(events+1)); %average substrate position
    avFc(jj,1)=sum((totFc.*timestep)/t(events+1)); %average force
    failures(jj,1)=(length(numcon(numcon==0))-1);
    cyct(jj,1)=mean(cyctime); %mean cycle time
    velcv(jj,1)=((sum(timestep.*((vel-retro(jj)).^2)/t(events+1)))^0.5)/retro(jj); %weighted actin flow coeffieicent of variation
    rebcyc(jj,1)=sum(crebind)./failures(jj);
end
figure;
subplot(1,2,1)
semilogx(stiffness,kon)
xlabel('Substrate stiffness (kPa)')
ylabel('kon (s-1)')

subplot(1,2,2)
semilogx(stiffness,avtrac)
xlabel('Substrate stiffness (kPa)')
ylabel('Mean traction force (pN)')

save('SR.mat')

toc