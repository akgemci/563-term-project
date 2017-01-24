%% EE563 Fall 2016 Term project 
% Aysel Akgemci-1875566
%% initilization
Sampfreq=25600;                             %sampling frequency given
n=0.2*Sampfreq;                             %number of samples
t=(0:(n-1)).*(1/Sampfreq);                  %time vector for given sample rate
%fundamental components
iA1=1750*sin(100*pi.*t)*sqrt(2);
iB1=1261.9*sin(100*pi.*t-0.7439*pi)*sqrt(2);
iC1=1261.9*sin(100*pi.*t+0.7439*pi)*sqrt(2);
%95Hz components
iA95=190*sin(190*pi.*t)*sqrt(2);
iB95=95.4*sin(190*pi.*t-0.9711*pi)*sqrt(2);
iC95=95.4*sin(190*pi.*t+0.9711*pi)*sqrt(2);
%100Hz components
iA2=154*sin(200*pi.*t)*sqrt(2);
iB2=77.95*sin(200*pi.*t-0.9503*pi)*sqrt(2);
iC2=77.95*sin(200*pi.*t+0.9503*pi)*sqrt(2);
%105Hz components
iA105=180*sin(210*pi.*t)*sqrt(2);
iB105=96.4*sin(210*pi.*t-0.883*pi)*sqrt(2);
iC105=96.4*sin(210*pi.*t+0.883*pi)*sqrt(2);
%total currents
iA=iA1+iA95+iA2+iA105;
iB=iB1+iB95+iB2+iB105;
iC=iC1+iC95+iC2+iC105;

 %plot(t,iA,t,iB,t,iC);
 %% part A
 %find alpha beta currents for fundamental component
 ialpha_50Hz = sqrt(2/3)*(iA-(iB/2)-(iC/2));
 ibeta_50Hz = sqrt(2/3)*(sqrt(3)*(iB/2)-sqrt(3)*(iC/2));
 izero_50Hz = sqrt(2/3)*(iA/(sqrt(2))+iB/(sqrt(2))+iC/(sqrt(2)));
 % find d q currents for both positive and negative directions
 idpos_50Hz = ialpha_50Hz.*cos(2*pi*50*t)+ibeta_50Hz.*sin(2*pi*50*t);
 iqpos_50Hz = -ialpha_50Hz.*sin(2*pi*50*t)+ibeta_50Hz.*cos(2*pi*50*t);
 
 idneg_50Hz = ialpha_50Hz.*cos(2*pi*50*t)-ibeta_50Hz.*sin(2*pi*50*t);
 iqneg_50Hz = ialpha_50Hz.*sin(2*pi*50*t)+ibeta_50Hz.*cos(2*pi*50*t);
 %find the mean instead of using LPF
 idposmean_50Hz = mean(idpos_50Hz);
 iqposmean_50Hz = mean(iqpos_50Hz);
 idnegmean_50Hz = mean(idneg_50Hz);
 iqnegmean_50Hz = mean(iqneg_50Hz);
 % perform reverse transformation
 ialpha_newpos50Hz = (idposmean_50Hz * cos(2*pi*50*t)) - (iqposmean_50Hz * sin(2*pi*50*t));
 ibeta_newpos50Hz = (idposmean_50Hz * sin(2*pi*50*t)) + (iqposmean_50Hz * cos(2*pi*50*t));
 
 ialpha_newneg50Hz = (idnegmean_50Hz * cos(2*pi*50*t)) - (iqnegmean_50Hz * sin((-2)*pi*50*t));
 ibeta_newneg50Hz = (idnegmean_50Hz * sin((-2)*pi*50*t)) + (iqnegmean_50Hz * cos(2*pi*50*t));
 % positive and negative sequence components of fundamental current
 ia_pos50Hz = sqrt(2/3) * ialpha_newpos50Hz;
 ia_neg50Hz = sqrt(2/3) * ialpha_newneg50Hz;
 
 ib_pos50Hz = sqrt(2/3) * ((-0.5) * ialpha_newpos50Hz + (sqrt(3)/2) * ibeta_newpos50Hz);
 ib_neg50Hz = sqrt(2/3) * ((-0.5) * ialpha_newneg50Hz + (sqrt(3)/2) * ibeta_newneg50Hz);
 
 ic_pos50Hz = sqrt(2/3) * ((-0.5) * ialpha_newpos50Hz - (sqrt(3)/2) * ibeta_newpos50Hz);
 ic_neg50Hz = sqrt(2/3) * ((-0.5) * ialpha_newneg50Hz - (sqrt(3)/2) * ibeta_newneg50Hz);

 %plot(t,ia_pos50Hz,t,ib_pos50Hz,t,ic_pos50Hz);
 %plot(t,ia_neg50Hz,t,ib_neg50Hz,t,ic_neg50Hz);
ia_50Hz = ia_pos50Hz+ia_neg50Hz;
ib_50Hz = ib_pos50Hz+ib_neg50Hz;
ic_50Hz = ic_pos50Hz+ic_neg50Hz;

iA_H = iA - ia_50Hz;
iB_H = iB - ib_50Hz;
iC_H = iC - ic_50Hz;


%% Perform the same procedure for 100Hz component
 ialpha_100Hz = sqrt(2/3)*(iA_H-(iB_H/2)-(iC_H/2));
 ibeta_100Hz = sqrt(2/3)*(sqrt(3)*(iB_H/2)-sqrt(3)*(iC_H/2));
 izero_100Hz = sqrt(2/3)*(iA_H/(sqrt(2))+iB_H/(sqrt(2))+iC_H/(sqrt(2)));
 
 idpos_100Hz = ialpha_100Hz.*cos(2*pi*100*t)+ibeta_100Hz.*sin(2*pi*100*t);
 iqpos_100Hz = -ialpha_100Hz.*sin(2*pi*100*t)+ibeta_100Hz.*cos(2*pi*100*t);
 
 idneg_100Hz = ialpha_100Hz.*cos(2*pi*100*t)-ibeta_100Hz.*sin(2*pi*100*t);
 iqneg_100Hz = ialpha_100Hz.*sin(2*pi*100*t)+ibeta_100Hz.*cos(2*pi*100*t);
 
 idposmean_100Hz = mean(idpos_100Hz);
 iqposmean_100Hz = mean(iqpos_100Hz);
 idnegmean_100Hz = mean(idneg_100Hz);
 iqnegmean_100Hz = mean(iqneg_100Hz);
 
 ialpha_newpos100Hz = (idposmean_100Hz * cos(2*pi*100*t)) - (iqposmean_100Hz * sin(2*pi*100*t));
 ibeta_newpos100Hz = (idposmean_100Hz * sin(2*pi*100*t)) + (iqposmean_100Hz * cos(2*pi*100*t));
 
 ialpha_newneg100Hz = (idnegmean_100Hz * cos(2*pi*100*t)) - (iqnegmean_100Hz * sin((-2)*pi*100*t));
 ibeta_newneg100Hz = (idnegmean_100Hz * sin((-2)*pi*100*t)) + (iqnegmean_100Hz * cos(2*pi*100*t));
 
 ia_pos100Hz = sqrt(2/3) * ialpha_newpos100Hz;
 ia_neg100Hz = sqrt(2/3) * ialpha_newneg100Hz;
 
 ib_pos100Hz = sqrt(2/3) * ((-0.5) * ialpha_newpos100Hz + (sqrt(3)/2) * ibeta_newpos100Hz);
 ib_neg100Hz = sqrt(2/3) * ((-0.5) * ialpha_newneg100Hz + (sqrt(3)/2) * ibeta_newneg100Hz);
 
 ic_pos100Hz = sqrt(2/3) * ((-0.5) * ialpha_newpos100Hz - (sqrt(3)/2) * ibeta_newpos100Hz);
 ic_neg100Hz = sqrt(2/3) * ((-0.5) * ialpha_newneg100Hz - (sqrt(3)/2) * ibeta_newneg100Hz);

 %plot(t,ia_pos100Hz,t,ib_pos100Hz,t,ic_pos100Hz);
 %plot(t,ia_neg100Hz,t,ib_neg100Hz,t,ic_neg100Hz);
 
ia_100Hz = ia_pos100Hz+ia_neg100Hz;
ib_100Hz = ib_pos100Hz+ib_neg100Hz;
ic_100Hz = ic_pos100Hz+ic_neg100Hz;

iA_iH = iA_H - ia_100Hz;
iB_iH = iB_H - ib_100Hz;
iC_iH = iC_H - ic_100Hz;
%plot(t,iA_iH,t,iB_iH,t,iC_iH)


%% Perform the same procedure for 105Hz component
 ialpha_105Hz = sqrt(2/3)*(iA_iH-(iB_iH/2)-(iC_iH/2));
 ibeta_105Hz = sqrt(2/3)*(sqrt(3)*(iB_iH/2)-sqrt(3)*(iC_iH/2));
 izero_105Hz = sqrt(2/3)*(iA_iH/(sqrt(2))+iB_iH/(sqrt(2))+iC_iH/(sqrt(2)));
 
 idpos_105Hz = ialpha_105Hz.*cos(2*pi*105*t)+ibeta_105Hz.*sin(2*pi*105*t);
 iqpos_105Hz = -ialpha_105Hz.*sin(2*pi*105*t)+ibeta_105Hz.*cos(2*pi*105*t);
 
 idneg_105Hz = ialpha_105Hz.*cos(2*pi*105*t)-ibeta_105Hz.*sin(2*pi*105*t);
 iqneg_105Hz = ialpha_105Hz.*sin(2*pi*105*t)+ibeta_105Hz.*cos(2*pi*105*t);
 
 idposmean_105Hz = mean(idpos_105Hz);
 iqposmean_105Hz = mean(iqpos_105Hz);
 idnegmean_105Hz = mean(idneg_105Hz);
 iqnegmean_105Hz = mean(iqneg_105Hz);
 
 ialpha_newpos105Hz = (idposmean_105Hz * cos(2*pi*105*t)) - (iqposmean_105Hz * sin(2*pi*105*t));
 ibeta_newpos105Hz = (idposmean_105Hz * sin(2*pi*105*t)) + (iqposmean_105Hz * cos(2*pi*105*t));
 
 ialpha_newneg105Hz = (idnegmean_105Hz * cos(2*pi*105*t)) - (iqnegmean_105Hz * sin((-2)*pi*105*t));
 ibeta_newneg105Hz = (idnegmean_105Hz * sin((-2)*pi*105*t)) + (iqnegmean_105Hz * cos(2*pi*105*t));
 
 ia_pos105Hz = sqrt(2/3) * ialpha_newpos105Hz;
 ia_neg105Hz = sqrt(2/3) * ialpha_newneg105Hz;
 
 ib_pos105Hz = sqrt(2/3) * ((-0.5) * ialpha_newpos105Hz + (sqrt(3)/2) * ibeta_newpos105Hz);
 ib_neg105Hz = sqrt(2/3) * ((-0.5) * ialpha_newneg105Hz + (sqrt(3)/2) * ibeta_newneg105Hz);
 
 ic_pos105Hz = sqrt(2/3) * ((-0.5) * ialpha_newpos105Hz - (sqrt(3)/2) * ibeta_newpos105Hz);
 ic_neg105Hz = sqrt(2/3) * ((-0.5) * ialpha_newneg105Hz - (sqrt(3)/2) * ibeta_newneg105Hz);

 %plot(t,ia_pos105Hz,t,ib_pos105Hz,t,ic_pos105Hz);
 %plot(t,ia_neg105Hz,t,ib_neg105Hz,t,ic_neg105Hz);
 
ia_105Hz = ia_pos105Hz+ia_neg105Hz;
ib_105Hz = ib_pos105Hz+ib_neg105Hz;
ic_105Hz = ic_pos105Hz+ic_neg105Hz;

iA_iiH = iA_iH - ia_105Hz;
iB_iiH = iB_iH - ib_105Hz;
iC_iiH = iC_iH - ic_105Hz;


%% Perform the same procedure for 95Hz component
 ialpha_95Hz = sqrt(2/3)*(iA_iiH-(iB_iiH/2)-(iC_iiH/2));
 ibeta_95Hz = sqrt(2/3)*(sqrt(3)*(iB_iiH/2)-sqrt(3)*(iC_iiH/2));
 izero_95Hz = sqrt(2/3)*(iA_iiH/(sqrt(2))+iB_iiH/(sqrt(2))+iC_iiH/(sqrt(2)));
 
 idpos_95Hz = ialpha_95Hz.*cos(2*pi*95*t)+ibeta_95Hz.*sin(2*pi*95*t);
 iqpos_95Hz = -ialpha_95Hz.*sin(2*pi*95*t)+ibeta_95Hz.*cos(2*pi*95*t);
 
 idneg_95Hz = ialpha_95Hz.*cos(2*pi*95*t)-ibeta_95Hz.*sin(2*pi*95*t);
 iqneg_95Hz = ialpha_95Hz.*sin(2*pi*95*t)+ibeta_95Hz.*cos(2*pi*95*t);
 
 idposmean_95Hz = mean(idpos_95Hz);
 iqposmean_95Hz = mean(iqpos_95Hz);
 idnegmean_95Hz = mean(idneg_95Hz);
 iqnegmean_95Hz = mean(iqneg_95Hz);
 
 ialpha_newpos95Hz = (idposmean_95Hz * cos(2*pi*95*t)) - (iqposmean_95Hz * sin(2*pi*95*t));
 ibeta_newpos95Hz = (idposmean_95Hz * sin(2*pi*95*t)) + (iqposmean_95Hz * cos(2*pi*95*t));
 
 ialpha_newneg95Hz = (idnegmean_95Hz * cos(2*pi*95*t)) - (iqnegmean_95Hz * sin((-2)*pi*95*t));
 ibeta_newneg95Hz = (idnegmean_95Hz * sin((-2)*pi*95*t)) + (iqnegmean_95Hz * cos(2*pi*95*t));
 
 ia_pos95Hz = sqrt(2/3) * ialpha_newpos95Hz;
 ia_neg95Hz = sqrt(2/3) * ialpha_newneg95Hz;
 
 ib_pos95Hz = sqrt(2/3) * ((-0.5) * ialpha_newpos95Hz + (sqrt(3)/2) * ibeta_newpos95Hz);
 ib_neg95Hz = sqrt(2/3) * ((-0.5) * ialpha_newneg95Hz + (sqrt(3)/2) * ibeta_newneg95Hz);
 
 ic_pos95Hz = sqrt(2/3) * ((-0.5) * ialpha_newpos95Hz - (sqrt(3)/2) * ibeta_newpos95Hz);
 ic_neg95Hz = sqrt(2/3) * ((-0.5) * ialpha_newneg95Hz - (sqrt(3)/2) * ibeta_newneg95Hz);

 %plot(t,ia_pos95Hz,t,ib_pos95Hz,t,ic_pos95Hz);
 %plot(t,ia_neg95Hz,t,ib_neg95Hz,t,ic_neg95Hz);
 
ia_95Hz = ia_pos95Hz+ia_neg95Hz;
ib_95Hz = ib_pos95Hz+ib_neg95Hz;
ic_95Hz = ic_pos95Hz+ic_neg95Hz;

iA_0H = iA_iiH - ia_95Hz;
iB_0H = iB_iiH - ib_95Hz;
iC_0H = iC_iiH - ic_95Hz;
%plot(t,iA_0H,t,iB_0H,t,iC_0H);

 %% part B analytical solution
 
 % phasor representations
 ia1=1750*sqrt(2);
 ib1=1261.9*sqrt(2)*exp(-i*0.7439*pi);
 ic1=1261.9*sqrt(2)*exp(i*0.7439*pi);
 
 ia95=190*sqrt(2);
 ib95=95.4*sqrt(2)*exp(-i*0.9711*pi);
 ic95=95.4*sqrt(2)*exp(i*0.9711*pi);
 
 ia2=154*sqrt(2);
 ib2=77.95*sqrt(2)*exp(-i*0.9503*pi);
 ic2=77.95*sqrt(2)*exp(i*0.9503*pi);
 
 ia105=180*sqrt(2);
 ib105=96.4*sqrt(2)*exp(-i*0.883*pi);
 ic105=96.4*sqrt(2)*exp(i*0.883*pi);
 
 a=-(1/2)+(sqrt(3)/2)*i;
 C=(1/3)*[1 1 1; 1 a a^2;1 a^2 a];
 I3ph1=[ia1 ;ib1 ;ic1];
 I3ph95=[ia95 ;ib95 ;ic95];
 I3ph2=[ia2 ;ib2 ;ic2];
 I3ph105=[ia105; ib105 ;ic105];
 
 Iseq1=C*I3ph1
 Iseq95=C*I3ph95
 Iseq2=C*I3ph2
 Iseq105=C*I3ph105
 Ia_analytical_pos=(Iseq1(2))*sin(2*pi*50*t);
 %plot(t,Ia_analytical_pos,t,ia_pos50Hz);
 
 
 %% part C
 
 d_ft_50Hz=fft(ia_pos50Hz);
 plot(d_ft_50Hz)
 
 

 
 
 


