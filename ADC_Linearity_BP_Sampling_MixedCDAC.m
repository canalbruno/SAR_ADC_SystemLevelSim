% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Date: 27/04/2022
% Author: Bruno Canal
% 
% Linearity extraction of a SAR ADC with redundant weights on the DAC array
% The extraction considers a sweep input signal of "StepNum" steps from -Vref
% to +Vref.
% The analysis considers the Gaussian Normal Distribution of capacitor 
% variabilities whith the variance of "VarCun" and "MCSamples"
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear all
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Enviromental and constant
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k = 1.3806504e-23;
T=293;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script Options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
KtCNoise='false'; % 'true' to simulate with KT/C Noise effect

StepNum=10000;
MCsample=10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ADC Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Cun=1;
VarCun=Cun*0.01;

numBits=10;
capWeight=[192 128 64 0 0 0 0 0 0 0 0];
cap1Weight=[0 0 0 32 16 7 4 2 1 1 1];
cap2Weight=[0 0 0 32 16 7 4 2 1 1 1];

CtotIdeal=sum(capWeight)+sum(cap1Weight)+sum(cap2Weight);
numCap=length(cap1Weight);
Vref=0.6;
Vcm=Vref/2;
LSB=2*Vref/2^(numBits);

INL_v6=zeros(numBits,MCsample);
DNL_v6=zeros(numBits,MCsample);

compOffset=0; % Comparator Offset

decVal=zeros(MCsample,StepNum);


for MC=1:MCsample
    
    CDAC_P=zeros(1,length(capWeight));
    CDAC_N=zeros(1,length(capWeight));
    CDAC_P1=zeros(1,length(cap1Weight));
    CDAC_P2=zeros(1,length(cap1Weight));
    CDAC_N1=zeros(1,length(cap2Weight));
    CDAC_N2=zeros(1,length(cap2Weight));
    
    %%%%%%%%%%%%%%%%%%%%%%% Generating the Cap Array %%%%%%%%%%%%%%%%%
    for capIndex=1:length(cap1Weight)
        for capMult=1:cap1Weight(capIndex)
            CDAC_P1(capIndex)=CDAC_P1(capIndex)+normrnd(Cun,VarCun);
            CDAC_N1(capIndex)=CDAC_N1(capIndex)+normrnd(Cun,VarCun);
        end
    end
    
    for capIndex=1:length(cap2Weight)
        for capMult=1:cap2Weight(capIndex)
            CDAC_P2(capIndex)=CDAC_P2(capIndex)+normrnd(Cun,VarCun);
            CDAC_N2(capIndex)=CDAC_N2(capIndex)+normrnd(Cun,VarCun);
        end
    end
    
    for capIndex=1:length(capWeight)
        for capMult=1:capWeight(capIndex)
            CDAC_P(capIndex)=CDAC_P(capIndex)+normrnd(Cun,VarCun);
            CDAC_N(capIndex)=CDAC_N(capIndex)+normrnd(Cun,VarCun);
        end
    end
    
    Ctot_P=sum(CDAC_P)+sum(CDAC_P1)+sum(CDAC_P2);
    Ctot_N=sum(CDAC_N)+sum(CDAC_N1)+sum(CDAC_N2);
       
    if KtCNoise
        SampNoiseP=sqrt(k*T/Ctot_P);
        SampNoiseN=sqrt(k*T/Ctot_N);
    else
        SampNoiseP=0;
        SampNoiseN=0;
    end

    %%%%%%%%%%%%%%%%%%%%%% Define the switching value %%%%%%%%%%%%%%%%

    
    binVal=zeros(StepNum,numCap-1);
    digVal=zeros(StepNum,numBits);
    Index2VFactor=((1200/(StepNum-1))-(600+1200/(StepNum-1)))/1000;
    
    
    Pswitch=zeros(StepNum,numCap);
    Nswitch=zeros(StepNum,numCap);
    P1switch=zeros(StepNum,numCap);
    P2switch=ones(StepNum,numCap);
    N1switch=zeros(StepNum,numCap);
    N2switch=ones(StepNum,numCap);
        
    for i=1:StepNum
        
        RSD_val=zeros(1,4);
        InputValue=(i*(1200/(StepNum-1))-(600+1200/(StepNum-1)))/1000;
        
        Vp=Vcm;
        Vn=Vcm;
        Vbp0_DACP=ones(1,length(capWeight))*(Vcm+InputValue/2);
        Vbp0_DACN=ones(1,length(capWeight))*(Vcm-InputValue/2); 
        Vbp0_DACP1=ones(1,length(cap1Weight))*(Vcm+InputValue/2);
        Vbp0_DACP2=ones(1,length(cap1Weight))*(Vcm+InputValue/2);
        Vbp0_DACN1=ones(1,length(cap1Weight))*(Vcm-InputValue/2);
        Vbp0_DACN2=ones(1,length(cap1Weight))*(Vcm-InputValue/2);
        
        Vbpf_DACP=Vcm+Pswitch(i,:)*Vcm;
        Vbpf_DACN=Vcm+Nswitch(i,:)*Vcm;
        Vbpf_DACP1=P1switch(i,:)*Vref;
        Vbpf_DACP2=P2switch(i,:)*Vref;
        Vbpf_DACN1=N1switch(i,:)*Vref;
        Vbpf_DACN2=N2switch(i,:)*Vref;
        Vp=Vp+normrnd(0,SampNoiseP)+((sum(CDAC_P.*(Vbpf_DACP-Vbp0_DACP))+sum(CDAC_P1.*(Vbpf_DACP1-Vbp0_DACP1))+sum(CDAC_P2.*(Vbpf_DACP2-Vbp0_DACP2)))/Ctot_P);
        Vn=Vn+normrnd(0,SampNoiseN)+((sum(CDAC_N.*(Vbpf_DACN-Vbp0_DACN))+sum(CDAC_N1.*(Vbpf_DACN1-Vbp0_DACN1))+sum(CDAC_N2.*(Vbpf_DACN2-Vbp0_DACN2)))/Ctot_N);
        Vd=Vp-Vn;
        
        %%%%% Window Swithcing %%%%%%%%        
        %%%% First SAR Cycle
        if (InputValue>(5*Vref/8)) || (InputValue<(-5*Vref/8))
%              display('Big')
           TDC=0;
           if Vd<compOffset
               binVal(i,1)=1;
               Pswitch(i,1:3)=1; 
               Nswitch(i,1:3)=-1;
           else
               binVal(i,1)=0;
               Pswitch(i,1:3)=-1;
               Nswitch(i,1:3)=1;  
           end
        elseif (InputValue>(3*Vref/8)) || (InputValue<(-3*Vref/8))
%              display('Mid')
            TDC=1;
            if Vd<compOffset
               binVal(i,1)=1;
               Pswitch(i,2:3)=1; 
               Nswitch(i,2:3)=-1;
           else
               binVal(i,1)=0;
               Pswitch(i,2:3)=-1;
               Nswitch(i,2:3)=1;  
            end
        else    
%              display('Small')
            TDC=2;
            if Vd<compOffset
               binVal(i,1)=1;
               Pswitch(i,3)=1; 
               Nswitch(i,3)=-1;
           else
               binVal(i,1)=0;
               Pswitch(i,3)=-1;
               Nswitch(i,3)=1;  
            end
        end
                    
        Vbp0_DACP=Vbpf_DACP;
        Vbp0_DACN=Vbpf_DACN;
        Vbp0_DACP1=Vbpf_DACP1;
        Vbp0_DACP2=Vbpf_DACP2;
        Vbp0_DACN1=Vbpf_DACN1;
        Vbp0_DACN2=Vbpf_DACN2;

        Vbpf_DACP=Vcm+Pswitch(i,:)*Vcm;
        Vbpf_DACN=Vcm+Nswitch(i,:)*Vcm;
        Vbpf_DACP1=P1switch(i,:)*Vref;
        Vbpf_DACP2=P2switch(i,:)*Vref;
        Vbpf_DACN1=N1switch(i,:)*Vref;
        Vbpf_DACN2=N2switch(i,:)*Vref;

        
        for capIndex=4:numCap
            bitIndex=capIndex-2;
            Vp=Vp+((sum(CDAC_P.*(Vbpf_DACP-Vbp0_DACP))+sum(CDAC_P1.*(Vbpf_DACP1-Vbp0_DACP1))+sum(CDAC_P2.*(Vbpf_DACP2-Vbp0_DACP2)))/Ctot_P);
            Vn=Vn+((sum(CDAC_N.*(Vbpf_DACN-Vbp0_DACN))+sum(CDAC_N1.*(Vbpf_DACN1-Vbp0_DACN1))+sum(CDAC_N2.*(Vbpf_DACN2-Vbp0_DACN2)))/Ctot_N);
            Vd=Vp-Vn;

         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            if Vd<compOffset
                binVal(i,bitIndex)=1;
                P1switch(i,capIndex)=1;
                if capIndex ~= numCap
                    N2switch(i,capIndex)=0;
                end
            else
                binVal(i,bitIndex)=0;
                if capIndex ~= numCap
                    P2switch(i,capIndex)=0;
                end
                N1switch(i,capIndex)=1;               
            end
            
            %%%%%%%% CRS %%%%%%%%%%%
            if (bitIndex==2 && (binVal(i,1) ~= binVal(i,2)))
                
                Pswitch(i,3)=0;
                Nswitch(i,3)=0;
                P1switch(i,4)=0;
                P2switch(i,4)=1;
                N1switch(i,4)=0; 
                N2switch(i,4)=1; 
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%
            
            Vbp0_DACP=Vbpf_DACP;
            Vbp0_DACN=Vbpf_DACN;
            Vbp0_DACP1=Vbpf_DACP1;
            Vbp0_DACP2=Vbpf_DACP2;
            Vbp0_DACN1=Vbpf_DACN1;
            Vbp0_DACN2=Vbpf_DACN2;
            
            Vbpf_DACP=Vcm+Pswitch(i,:)*Vcm;
            Vbpf_DACN=Vcm+Nswitch(i,:)*Vcm;
            Vbpf_DACP1=P1switch(i,:)*Vref;
            Vbpf_DACP2=P2switch(i,:)*Vref;
            Vbpf_DACN1=N1switch(i,:)*Vref;
            Vbpf_DACN2=N2switch(i,:)*Vref;
            
        end
        
        %%%% Last SAR comparison %%%%
        Vp=Vp+((sum(CDAC_P.*(Vbpf_DACP-Vbp0_DACP))+sum(CDAC_P1.*(Vbpf_DACP1-Vbp0_DACP1))+sum(CDAC_P2.*(Vbpf_DACP2-Vbp0_DACP2)))/Ctot_P);
        Vn=Vn+((sum(CDAC_N.*(Vbpf_DACN-Vbp0_DACN))+sum(CDAC_N1.*(Vbpf_DACN1-Vbp0_DACN1))+sum(CDAC_N2.*(Vbpf_DACN2-Vbp0_DACN2)))/Ctot_N);
        Vd=Vp-Vn;
        Residue(MC,i)=Vd;
        
        if Vd<compOffset
            binVal(i,capIndex-1)=1;
        else
            binVal(i,capIndex-1)=0;
        end
        
        %%%% Digital Error Correction (DEC) - For redundancy recovery %%%%
        %%%% TDC
        switch TDC
            case 0 % LARGE AMPLITUDE
                RSD_val(1)=binVal(i,1); %D<9>
                RSD_val(2)=binVal(i,1); %D<8>
                RSD_val(3)=binVal(i,2); %D<7>
                RSD_val(4)=binVal(i,3); %D<6>
            case 1 % MEDIUM AMPLITUDE
                RSD_val(1)=binVal(i,1); %D<9>
                RSD_val(2)=xor(xor(binVal(i,1),binVal(i,2)),binVal(i,3)); %D<8>
                RSD_val(3)=~RSD_val(2); %D<7>
                RSD_val(4)=~binVal(i,3); %D<6>
            case 2 % SMALL AMPLITUDE
                Sum3Bit=binVal(i,1)+binVal(i,2)+binVal(i,3);
                if Sum3Bit>1 %D<9>
                    RSD_val(1)=1;
                else
                    RSD_val(1)=0;
                end
                RSD_val(2)=~RSD_val(1); %D<8>
                RSD_val(3)=xor(xor(binVal(i,1),binVal(i,2)),binVal(i,3)); %D<7>
                RSD_val(4)=~binVal(i,3); %D<6>
            otherwise
                RSD_val(1:5)=1;
                display('Error')
        end
   
        digVal(i,numBits)=binVal(i,numBits); %D<0>
        digVal(i,numBits-1)=mod(binVal(i,numBits-1)+binVal(i,numBits-2)+binVal(i,numBits-2),2); %D<1>
        if (binVal(i,numBits-1)+binVal(i,numBits-2)+binVal(i,numBits-2))>1; carry=1; else carry=0; end
            
        digVal(i,numBits-2)=mod(binVal(i,numBits-6)+binVal(i,numBits-3)+carry,2);%D<2>
        if (binVal(i,numBits-6)+binVal(i,numBits-3)+carry)>1; carry=1; else carry=0; end
                        
        digVal(i,numBits-3)=mod(binVal(i,numBits-4)+binVal(i,numBits-6)+carry,2);%D<3>        
        if (binVal(i,numBits-4)+binVal(i,numBits-6)+carry)>1; carry=1; else carry=0; end
            
        digVal(i,numBits-4)=mod(binVal(i,numBits-5)+binVal(i,numBits-6)+carry,2); %D<4>       
        if (binVal(i,numBits-5)+binVal(i,numBits-6)+carry)>1; carry=1; else carry=0; end
            
        digVal(i,numBits-5)=carry; %D<5>
        
        digVal(i,numBits-6)=RSD_val(4); %D<6>
        digVal(i,numBits-7)=RSD_val(3); %D<7>            
        digVal(i,numBits-8)=RSD_val(2); %D<8>        
        digVal(i,numBits-9)=RSD_val(1); %D<9>
            
        decVal(MC,i)=bi2de(digVal(i,:),'left-msb');
    end
    
    %%%% Find the transistions of Output ADC Value %%%%
    trValue=zeros(1,2^(numBits));
    auxVal=decVal(MC,1);
    trValue(1)=-0.6;
    trValueIndex=2;
    
    for indexSearch=1:StepNum
       if auxVal<decVal(MC,indexSearch)
           trValue(trValueIndex)=(indexSearch*(1200/(StepNum-1))-(600+1200/(StepNum-1)))/1000;
           save(trValueIndex)=indexSearch;
           auxVal=decVal(MC,indexSearch);
           trValueIndex=trValueIndex+1;
       end        
    end
    
    for i=3:2^(numBits)-1
        DNL_v6(i-2,MC)=(((trValue(i)-trValue(i-1))-LSB)/(LSB));
        INL_v6(i-2,MC)=sum(DNL_v6(:,MC));
    end
    
end



for i=1:length(INL_v6(:,1))
    RMS_INL_v6(i)=rms(INL_v6(i,:));
    RMS_DNL_v6(i)=rms(DNL_v6(i,:));
end
Max_RMSBGWININL=max(RMS_INL_v6)
Max_RMSBGWINDNL=max(RMS_DNL_v6)

figure (6)
    subplot(4,1,1)
    plot(DNL_v6)
    ylabel('DNL')
    xlim([0 2^(numBits)])
    subplot(4,1,2)
    plot(INL_v6)
    ylabel('INL')
    xlim([0 2^(numBits)])
    subplot(4,1,3)
    plot(RMS_INL_v6)
    ylabel('RMS INL')
    xlim([0 2^(numBits)])
    %ylim([0 0.4])
    subplot(4,1,4)
    plot(RMS_DNL_v6)
    ylabel('RMS DNL')
    xlim([0 2^(numBits)])
    %ylim([0 0.4])
    
figure (61)
    hold on
    plot(decVal')
    hold off
    
       

    
    