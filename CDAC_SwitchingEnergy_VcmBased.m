% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Date: 09/01/2024
% Author: Bruno Canal
%
% System level switching energy of SAR ADC employing Vcm-Based switching
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc

Cun=1; % 1 for enrgy in CV^2, choose a value for real energy calculation
Vdd=1; % 1 for enrgy in CV^2, choose a value for real energy calculation
Vref=Vdd;
Vcm=Vref/2;
Nbits=10;
LSB=Vref/(2^Nbits);

CDAC_P=[192 128 64 64 32 14 8 4 2 2 1];
CDAC_N=[192 128 64 64 32 14 8 4 2 2 1];

Ctot_P=sum(CDAC_P);
Ctot_N=sum(CDAC_N);

numCap=length(CDAC_P);

VDAC_P=zeros(1,numCap-2);
VDAC_N=zeros(1,numCap-2);

%%%%%%%%%%%%%%%%%%%%%% Define the switching value %%%%%%%%%%%%%%%%
   
VindArray=(-Vref:Vref/2^(Nbits-1):Vref);
E_switch_thesis=zeros(1,length(VindArray));

for i=1:length(VindArray)
    
    VDAC_P=zeros(1,numCap);
    VDAC_N=zeros(1,numCap);
    
    Ecode=0;
    Ep=zeros(1,numCap);
    En=zeros(1,numCap); 
    InputValue=VindArray(i);
    
    Pswitch=zeros(1,length(CDAC_P));
    Nswitch=zeros(1,length(CDAC_P));
    
    DvP=zeros(1, numCap);
    DvN=zeros(1, numCap);
        
    %%% Sampling %%%
    
    Vp=Vcm;
    Vn=Vcm;
    Vbp0_DACP=ones(1,length(CDAC_P))*(Vcm+InputValue/2);
    Vbp0_DACN=ones(1,length(CDAC_P))*(Vcm-InputValue/2); 

    Vbpf_DACP=Vcm+Pswitch*Vcm;
    Vbpf_DACN=Vcm+Nswitch*Vcm;

    VDAC_P(1)=Vp;
    VDAC_N(1)=Vn;

    Vp=VDAC_P(1)+((sum(CDAC_P.*(Vbpf_DACP-Vbp0_DACP)))/Ctot_P);
    Vn=VDAC_N(1)+((sum(CDAC_N.*(Vbpf_DACN-Vbp0_DACN)))/Ctot_N);
    
    Vd=Vp-Vn;

    %%% Remaining SAR Cycles

    for capIndex=1:length(CDAC_P)
        bitIndex=capIndex;

        if Vd<0
            binVal(bitIndex)=1;
            Pswitch(capIndex)=1;
            Nswitch(capIndex)=-1;
        else
            binVal(bitIndex)=0;
            Nswitch(capIndex)=1;   
            Pswitch(capIndex)=-1;
        end

        Vbp0_DACP=Vbpf_DACP;
        Vbp0_DACN=Vbpf_DACN;

        Vbpf_DACP=Vcm+Pswitch*Vcm;
        Vbpf_DACN=Vcm+Nswitch*Vcm;

        VDAC_P(bitIndex)=Vp;
        VDAC_N(bitIndex)=Vn;

        Vp=VDAC_P(bitIndex)+((sum(CDAC_P.*(Vbpf_DACP-Vbp0_DACP)))/Ctot_P);
        Vn=VDAC_N(bitIndex)+((sum(CDAC_N.*(Vbpf_DACN-Vbp0_DACN)))/Ctot_N);
        
        DvP=(Vp-Vbpf_DACP)-(VDAC_P(bitIndex)-Vbp0_DACP);
        DvN=(Vn-Vbpf_DACN)-(VDAC_N(bitIndex)-Vbp0_DACN);
        
        Ep=(-Vbpf_DACP.*CDAC_P).*DvP;
        En=(-Vbpf_DACN.*CDAC_N).*DvN;
      
        Vd=Vp-Vn;

        Ecode=Ecode+sum(Ep)+sum(En);
    end

     E_switch_thesis(i)=Ecode;
end

Avg_energy_thesis=mean(E_switch_thesis)
   
%%%%%%%%%%%%%%% VDAC Sequence %%%%%%%%%%%%%%%%%%%%%%%%%     
figure(1)
hold on
stairs(VDAC_P,'r')
stairs(VDAC_N)
axis([1 10  0 1])
hold off

%%%%%%%%%%%%%%% Switching Energy %%%%%%%%%%%%%%%%%%%%%%%%%    
figure(2)
hold on
plot(E_switch_thesis,'r')
axis([0 1024 0 250])
hold off

    
    