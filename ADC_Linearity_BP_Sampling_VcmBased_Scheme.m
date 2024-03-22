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
% Script Options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
StepNum=10000;
MCsample=100;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ADC Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Cun=1;
VarCun=Cun*0.01;

numBits=10;
capWeight=[192 128 64 64 32 14 8 4 2 2 1 1]; % the last capacitor works just in the sampling phase
CtotIdeal=sum(capWeight)+1;
numCap=length(capWeight);
Vref=0.6;
Vcm=Vref/2;
LSB=2*Vref/2^(numBits);


INL_v44=zeros(numBits,MCsample);
DNL_v44=zeros(numBits,MCsample);


for MC=1:MCsample
    
    CDAC_P=zeros(1,numCap);
    CDAC_N=zeros(1,numCap);
    
    %%%%%%%%%%%%%%%%%%%%%%% Generating the Cap Array %%%%%%%%%%%%%%%%%
    for capIndex=1:numCap
        for capMult=1:capWeight(capIndex)
            CDAC_P(capIndex)=CDAC_P(capIndex)+normrnd(Cun,VarCun);
            CDAC_N(capIndex)=CDAC_N(capIndex)+normrnd(Cun,VarCun);
        end
    end
    
    Ctot_P=sum(CDAC_P);%+normrnd(Cun,VarCun);
    Ctot_N=sum(CDAC_N);%+normrnd(Cun,VarCun);

    %%%%%%%%%%%%%%%%%%%%%% Define the switching value %%%%%%%%%%%%%%%%
    
    binVal=zeros(StepNum,numCap-1);
    digVal=zeros(StepNum,numBits);
    decVal=zeros(1,StepNum);
    Index2VFactor=((1200/(StepNum-1))-(600+1200/(StepNum-1)))/1000;
    
    
    for i=1:StepNum
        Pswitch=zeros(1,numCap);
        Nswitch=zeros(1,numCap);
        RSD_val=zeros(1,5);
        InputValue=(i*(1200/(StepNum-1))-(600+1200/(StepNum-1)))/1000;
        
        % At sampling cycle Top plate receive Vcm
        Vp=Vcm;
        Vn=Vcm;
        
        % At sampling cycle Bottom-plate receive Input signal
        Vbp0_DACP=ones(1,length(capWeight))*(Vcm+InputValue/2);
        Vbp0_DACN=ones(1,length(capWeight))*(Vcm-InputValue/2); 
        
        % At the end of sampling cycle bottom-plate receive Vcm (P and N switch are 0)
        Vbpf_DACP=Vcm+Pswitch*Vcm;
        Vbpf_DACN=Vcm+Nswitch*Vcm;
    
        for capIndex=1:numCap-1
            Vp=Vp+((sum(CDAC_P.*(Vbpf_DACP-Vbp0_DACP)))/Ctot_P);
            Vn=Vn+((sum(CDAC_N.*(Vbpf_DACN-Vbp0_DACN)))/Ctot_N);
            Vd=Vp-Vn;

         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            if Vd<0
                binVal(i,capIndex)=1;
                Pswitch(capIndex)=1; 
                Nswitch(capIndex)=-1;              
            else
                binVal(i,capIndex)=0;
                Pswitch(capIndex)=-1;
                Nswitch(capIndex)=1;                
            end
            
            Vbp0_DACP=Vbpf_DACP;
            Vbp0_DACN=Vbpf_DACN;
            Vbpf_DACP=Vcm+Pswitch*Vcm;
            Vbpf_DACN=Vcm+Nswitch*Vcm;
        end
        
        %%%% Last SAR comparison %%%%            
        
        Vp=Vp+((sum(CDAC_P.*(Vbpf_DACP-Vbp0_DACP)))/Ctot_P);
        Vn=Vn+((sum(CDAC_N.*(Vbpf_DACN-Vbp0_DACN)))/Ctot_N);
        Vd=Vp-Vn;
        
        if Vd<0
            binVal(i,capIndex+1)=1;
        else
            binVal(i,capIndex+1)=0;
        end
        

        %%%% Digital Error Correction (DEC) - For redundancy recovery %%%%
        
        digVal(i,numBits)=binVal(i,numCap); %D<0>
        digVal(i,numBits-1)=binVal(i,numCap-1); %D<1>             
        digVal(i,numBits-2)=mod(binVal(i,numCap-6)+binVal(i,numCap-3)+binVal(i,numCap-2),2);%D<2>
        if (binVal(i,numCap-6)+binVal(i,numCap-3)+binVal(i,numCap-2))>1; carry=1; else carry=0; end
                       
        digVal(i,numBits-3)=mod(binVal(i,numCap-4)+binVal(i,numCap-6)+carry,2);%D<3>        
        if (binVal(i,numCap-4)+binVal(i,numCap-6)+carry)>1; carry=1; else carry=0; end
            
        digVal(i,numBits-4)=mod(binVal(i,numCap-5)+binVal(i,numCap-6)+carry,2); %D<4>       
        if (binVal(i,numCap-5)+binVal(i,numCap-6)+carry)>1; carry=1; else carry=0; end
            
        digVal(i,numBits-5)=carry; %D<5>
        
        digVal(i,numBits-6)=binVal(i,numCap-7); %D<6>
        
        digVal(i,numBits-7)=mod(binVal(i,numCap-8)+binVal(i,numCap-9)+binVal(i,numCap-11),2);%D<7>
        if (binVal(i,numCap-8)+binVal(i,numCap-9)+binVal(i,numCap-11))>1; carry=1; else carry=0; end
        
        digVal(i,numBits-8)=mod(binVal(i,numCap-10)+binVal(i,numCap-11)+carry,2);%D<8>
        if (binVal(i,numCap-10)+binVal(i,numCap-11)+carry)>1; carry=1; else carry=0; end
              
        digVal(i,numBits-9)=carry; %D<9>
            
        decVal(i)=bi2de(digVal(i,:),'left-msb');
                
    end
    
    %%%% Find the transistions of Output ADC Value %%%%
    trValue=zeros(1,2^(numBits));
    auxVal=decVal(1);
    trValue(1)=-0.6;
    trValueIndex=2;
    
    for indexSearch=1:StepNum
       if auxVal<decVal(indexSearch)
           trValue(trValueIndex)=(indexSearch*(1200/(StepNum-1))-(600+1200/(StepNum-1)))/1000;
           save(trValueIndex)=indexSearch;
           auxVal=decVal(indexSearch);
           trValueIndex=trValueIndex+1;
       end        
    end
    
    for i=3:2^(numBits)-1
        DNL_v44(i-2,MC)=(((trValue(i)-trValue(i-1))-LSB)/(LSB));
        INL_v44(i-2,MC)=sum(DNL_v44(:,MC));
    end
    
end



for i=1:length(INL_v44(:,1))
    RMS_INL_v44(i)=rms(INL_v44(i,:));
    RMS_DNL_v44(i)=rms(DNL_v44(i,:));
end
Max_RMSBGWININL=max(RMS_INL_v44)
Max_RMSBGWINDNL=max(RMS_DNL_v44)

figure (44)
    subplot(4,1,1)
    plot(DNL_v44)
    ylabel('DNL')
    xlim([0 2^(numBits)])
    subplot(4,1,2)
    plot(INL_v44)
    ylabel('INL')
    xlim([0 2^(numBits)])
    subplot(4,1,3)
    plot(RMS_INL_v44)
    ylabel('RMS INL')
    xlim([0 2^(numBits)])
    subplot(4,1,4)
    plot(RMS_DNL_v44)
    ylabel('RMS DNL')
    xlim([0 2^(numBits)])
    
    

    
    