% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Date: 09/01/2024
% Author: Bruno Canal
%
% System level switching energy of SAR ADC employing window switching and 
% entire CDAC splitting strategy
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all
%close all
clc

Cun=1; % 1 for enrgy in CV^2, choose a value for real energy calculation
Vdd=1; % 1 for enrgy in CV^2, choose a value for real energy calculation

Vref=Vdd;
Vcm=Vref/2;
Nbits=10;
LSB=Vref/(2^Nbits);

%%%%%%%%%%%%%%%%%%%%%% CDAC weight distribution %%%%%%%%%%%%%%%%%%%%%%%%%%%
CDAC_P1=[96 64 32 32 16 7 4 2 1 1 1];
CDAC_P2=[96 64 32 32 16 7 4 2 1 1 1];
CDAC_N1=[96 64 32 32 16 7 4 2 1 1 1];
CDAC_N2=[96 64 32 32 16 7 4 2 1 1 1];

Ctot_P=sum(CDAC_P1)+sum(CDAC_P2);
Ctot_N=sum(CDAC_N1)+sum(CDAC_N2);

numCap=length(CDAC_P1);

VDAC_P=zeros(1,numCap-2);
VDAC_N=zeros(1,numCap-2);

%%%%%%%%%%%%%%%%%%%%%% Define the switching value %%%%%%%%%%%%%%%%
   
VindArray=(-Vref:Vref/2^(Nbits-1):Vref);
E_switch_thesis=zeros(1,length(VindArray));

for i=1:length(VindArray)
    
    Ecode=0;
    Ep=zeros(1,numCap);
    En=zeros(1,numCap);  
    
    InputValue=VindArray(i);
    
    P1switch=zeros(1,length(CDAC_P1));
    P2switch=ones(1,length(CDAC_P2));
    N1switch=zeros(1,length(CDAC_N1));
    N2switch=ones(1,length(CDAC_N2));
    
    DvP1=zeros(1, numCap);
    DvN1=zeros(1, numCap);
    DvP2=zeros(1, numCap);
    DvN2=zeros(1, numCap);
        
    %%% Sampling %%%
    
    Vp=Vcm;
    Vn=Vcm;
    Vbp0_DACP1=ones(1,length(CDAC_P1))*(Vcm+InputValue/2);
    Vbp0_DACP2=ones(1,length(CDAC_P1))*(Vcm+InputValue/2);
    Vbp0_DACN1=ones(1,length(CDAC_P1))*(Vcm-InputValue/2);
    Vbp0_DACN2=ones(1,length(CDAC_P1))*(Vcm-InputValue/2); 

    Vbpf_DACP1=P1switch*Vref;
    Vbpf_DACP2=P2switch*Vref;
    Vbpf_DACN1=N1switch*Vref;
    Vbpf_DACN2=N2switch*Vref;

    VDAC_P(1)=Vp;
    VDAC_N(1)=Vn;
    
    Vp=VDAC_P(1)+((sum(CDAC_P1.*(Vbpf_DACP1-Vbp0_DACP1))+sum(CDAC_P2.*(Vbpf_DACP2-Vbp0_DACP2)))/Ctot_P);
    Vn=VDAC_N(1)+((sum(CDAC_N1.*(Vbpf_DACN1-Vbp0_DACN1))+sum(CDAC_N2.*(Vbpf_DACN2-Vbp0_DACN2)))/Ctot_N);
   
    Vd=Vp-Vn;

    %%% Window Swithcing %%%  
    
    %%% First SAR Cycle
    if (InputValue>(5*Vref/8)) || (InputValue<(-5*Vref/8))
       TDC=0;
       if Vd<0
           binVal(1)=1;
           P1switch(1:3)=1;
           N2switch(1:3)=0;
       else
           binVal(1)=0;
           N1switch(1:3)=1;
           P2switch(1:3)=0;
       end
       Vbp0_DACP1=Vbpf_DACP1;
       Vbp0_DACP2=Vbpf_DACP2;
       Vbp0_DACN1=Vbpf_DACN1;
       Vbp0_DACN2=Vbpf_DACN2;
       
       Vbpf_DACP1=P1switch*Vref;
       Vbpf_DACP2=P2switch*Vref;
       Vbpf_DACN1=N1switch*Vref;
       Vbpf_DACN2=N2switch*Vref;
       
       VDAC_P(1)=Vp;
       VDAC_N(1)=Vn;
       
       Vp=VDAC_P(1)+((sum(CDAC_P1.*(Vbpf_DACP1-Vbp0_DACP1))+sum(CDAC_P2.*(Vbpf_DACP2-Vbp0_DACP2)))/Ctot_P);
       Vn=VDAC_N(1)+((sum(CDAC_N1.*(Vbpf_DACN1-Vbp0_DACN1))+sum(CDAC_N2.*(Vbpf_DACN2-Vbp0_DACN2)))/Ctot_N);
     
    elseif (InputValue>(3*Vref/8)) || (InputValue<(-3*Vref/8))
        TDC=1;
        if Vd<0
           binVal(1)=1;
           P1switch(2:3)=1;
           N2switch(2:3)=0;
        else
           binVal(1)=0;
           N1switch(2:3)=1;
           P2switch(2:3)=0;
        end
        Vbp0_DACP1=Vbpf_DACP1;
        Vbp0_DACP2=Vbpf_DACP2;
        Vbp0_DACN1=Vbpf_DACN1;
        Vbp0_DACN2=Vbpf_DACN2;
        
        Vbpf_DACP1=P1switch*Vref;
        Vbpf_DACP2=P2switch*Vref;
        Vbpf_DACN1=N1switch*Vref;
        Vbpf_DACN2=N2switch*Vref;
        
        VDAC_P(1)=Vp;
        VDAC_N(1)=Vn;
        
        Vp=VDAC_P(1)+((sum(CDAC_P1.*(Vbpf_DACP1-Vbp0_DACP1))+sum(CDAC_P2.*(Vbpf_DACP2-Vbp0_DACP2)))/Ctot_P);
        Vn=VDAC_N(1)+((sum(CDAC_N1.*(Vbpf_DACN1-Vbp0_DACN1))+sum(CDAC_N2.*(Vbpf_DACN2-Vbp0_DACN2)))/Ctot_N);
       
    else    
        TDC=2;
        if Vd<0
           binVal(1)=1;
           P1switch(3)=1;
           N2switch(3)=0;
        else
           binVal(1)=0;
           N1switch(3)=1;
           P2switch(3)=0; 
        end
        
        Vbp0_DACP1=Vbpf_DACP1;
        Vbp0_DACP2=Vbpf_DACP2;
        Vbp0_DACN1=Vbpf_DACN1;
        Vbp0_DACN2=Vbpf_DACN2;
        
        Vbpf_DACP1=P1switch*Vref;
        Vbpf_DACP2=P2switch*Vref;
        Vbpf_DACN1=N1switch*Vref;
        Vbpf_DACN2=N2switch*Vref;
        
        VDAC_P(1)=Vp;
        VDAC_N(1)=Vn;
        
        Vp=VDAC_P(1)+((sum(CDAC_P1.*(Vbpf_DACP1-Vbp0_DACP1))+sum(CDAC_P2.*(Vbpf_DACP2-Vbp0_DACP2)))/Ctot_P);
        Vn=VDAC_N(1)+((sum(CDAC_N1.*(Vbpf_DACN1-Vbp0_DACN1))+sum(CDAC_N2.*(Vbpf_DACN2-Vbp0_DACN2)))/Ctot_N);
        
    end
    
    DvP1=(Vp-Vbpf_DACP1)-(VDAC_P(1)-Vbp0_DACP1);
    DvN1=(Vn-Vbpf_DACN1)-(VDAC_N(1)-Vbp0_DACN1);
    DvP2=(Vp-Vbpf_DACP2)-(VDAC_P(1)-Vbp0_DACP2);
    DvN2=(Vn-Vbpf_DACN2)-(VDAC_N(1)-Vbp0_DACN2);  
    
    Ep=-(Vbpf_DACP1.*CDAC_P1).*DvP1-(Vbpf_DACP2.*CDAC_P2).*DvP2;
    En=-(Vbpf_DACN1.*CDAC_N1).*DvN1-(Vbpf_DACN2.*CDAC_N2).*DvN2;  
    
    Vd=Vp-Vn;
    Ecode=Ecode+sum(Ep)+sum(En);


    for capIndex=4:length(CDAC_P1)
        bitIndex=capIndex-2;
        
        DvP1=zeros(1, numCap);
        DvN1=zeros(1, numCap);
        DvP2=zeros(1, numCap);
        DvN2=zeros(1, numCap);
        
        if Vd<0
            binVal(bitIndex)=1;   
            P1switch(capIndex)=1;
            if capIndex ~= numCap
                N2switch(capIndex)=0;
            end
        else
            binVal(bitIndex)=0;
            N1switch(capIndex)=1;   
            if capIndex ~= numCap
                P2switch(capIndex)=0;
            end            
        end
        %%%%%%% CRS %%%%%%%%%%%
        if (bitIndex==2 && (binVal(1) ~= binVal(2)))
            P1switch(3:4)=zeros(1,2);
            P2switch(3:4)=ones(1,2);
            N1switch(3:4)=zeros(1,2);
            N2switch(3:4)=ones(1,2);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%

        Vbp0_DACP1=Vbpf_DACP1;
        Vbp0_DACP2=Vbpf_DACP2;
        Vbp0_DACN1=Vbpf_DACN1;
        Vbp0_DACN2=Vbpf_DACN2;

        Vbpf_DACP1=P1switch*Vref;
        Vbpf_DACP2=P2switch*Vref;
        Vbpf_DACN1=N1switch*Vref;
        Vbpf_DACN2=N2switch*Vref;

        VDAC_P(bitIndex)=Vp;
        VDAC_N(bitIndex)=Vn;

        Vp=VDAC_P(bitIndex)+((sum(CDAC_P1.*(Vbpf_DACP1-Vbp0_DACP1))+sum(CDAC_P2.*(Vbpf_DACP2-Vbp0_DACP2)))/Ctot_P);
        Vn=VDAC_N(bitIndex)+((sum(CDAC_N1.*(Vbpf_DACN1-Vbp0_DACN1))+sum(CDAC_N2.*(Vbpf_DACN2-Vbp0_DACN2)))/Ctot_N);
        
        DvP1=(Vp-Vbpf_DACP1)-(VDAC_P(bitIndex)-Vbp0_DACP1);
        DvN1=(Vn-Vbpf_DACN1)-(VDAC_N(bitIndex)-Vbp0_DACN1);
        DvP2=(Vp-Vbpf_DACP2)-(VDAC_P(bitIndex)-Vbp0_DACP2);
        DvN2=(Vn-Vbpf_DACN2)-(VDAC_N(bitIndex)-Vbp0_DACN2);
        
        Ep=-(Vbpf_DACP1.*CDAC_P1).*DvP1-(Vbpf_DACP2.*CDAC_P2).*DvP2;
        En=-(Vbpf_DACN1.*CDAC_N1).*DvN1-(Vbpf_DACN2.*CDAC_N2).*DvN2;

        Vd=Vp-Vn;
        Ecode=Ecode+sum(Ep)+sum(En);

    end

     E_switch_thesis(i)=Ecode;
end

Avg_energy_thesis=mean(E_switch_thesis)
   
    
figure(1)
hold on
stairs(VDAC_P,'r')
stairs(VDAC_N)
axis([1 10  0 1])
hold off

% %%%%%%%%%%%%%%% Array Simulation %%%%%%%%%%%%%%%%%%%%%%%%%    

figure(2)
hold on
plot(E_switch_thesis,'g')
axis([0 1024 0 150])
hold off

    
    