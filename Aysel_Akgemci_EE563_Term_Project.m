%% EE563 Fall 2016 Term project 
% Aysel Akgemci-1875566
% uncomment the corresponding plot segment to reobtain the graphs
%% initilization
Sampfreq=25600;            %sampling frequency given
n=0.2*Sampfreq;            %number of samples
t=(0:(n-1)).*(1/Sampfreq); %time vector for given sample rate
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
 %sum of all positive sequence components for each phases
 ia_POS = ia_pos50Hz + ia_pos100Hz + ia_pos95Hz + ia_pos105Hz;
 ib_POS = ib_pos50Hz + ib_pos100Hz + ib_pos95Hz + ib_pos105Hz;
 ic_POS = ic_pos50Hz + ic_pos100Hz + ic_pos95Hz + ic_pos105Hz;
 
 %sum of all negative sequence components for each phases
 ia_NEG = ia_neg50Hz + ia_neg100Hz + ia_neg95Hz + ia_neg105Hz;
 ib_NEG = ib_neg50Hz + ib_neg100Hz + ib_neg95Hz + ib_neg105Hz;
 ic_NEG = ic_neg50Hz + ic_neg100Hz + ic_neg95Hz + ic_neg105Hz;

 
%ia positive sequence
ia_pos_FFT = fft(ia_POS);
iaposfreq = (0:length(ia_pos_FFT)-1)*25600/length(ia_pos_FFT);
% generate empty template for filtered components
filter_50Hz_iapos = zeros(1,length(ia_pos_FFT));
filter_95Hz_iapos = zeros(1,length(ia_pos_FFT));
filter_100Hz_iapos = zeros(1,length(ia_pos_FFT));
filter_105Hz_iapos = zeros(1,length(ia_pos_FFT));
% assign the desired part of array to the empty template
% FFT is symmetric, two elemets need to exist in the filtered current FFTs
filter_50Hz_iapos(11) = ia_pos_FFT(11);
filter_50Hz_iapos(5120-11+2) = ia_pos_FFT(5120-11+2);
filter_95Hz_iapos(20) = ia_pos_FFT(20);
filter_95Hz_iapos(5120-20+2) = ia_pos_FFT(5120-20+2);
filter_100Hz_iapos(21) = ia_pos_FFT(21);
filter_100Hz_iapos(5120-21+2) = ia_pos_FFT(5120-21+2);
filter_105Hz_iapos(22) = ia_pos_FFT(22);
filter_105Hz_iapos(5120-22+2) = ia_pos_FFT(5120-22+2);

% figure;
% plot(t,ifft(filter_50Hz_iapos));
% hold on;
% plot(t,ifft(filter_95Hz_iapos));
% plot(t,ifft(filter_100Hz_iapos));
% plot(t,ifft(filter_105Hz_iapos));
% legend('50Hz comp','95Hz comp','100Hz comp','105Hz comp')
% xlabel('time');
% ylabel('I(A)');

%ia negative sequence
ia_neg_FFF = fft(ia_NEG);
ianegfreq = (0:length(ia_neg_FFF)-1)*25600/length(ia_neg_FFF);
% generate empty template for filtered components
filter_50Hz_ianeg = zeros(1,length(ia_neg_FFF));
filter_95Hz_ianeg = zeros(1,length(ia_neg_FFF));
filter_100Hz_ianeg = zeros(1,length(ia_neg_FFF));
filter_105Hz_ianeg = zeros(1,length(ia_neg_FFF));
% assign the desired part of array to the empty template
% FFT is symmetric, two elemets need to exist in the filtered current FFTs
filter_50Hz_ianeg(11) = ia_neg_FFF(11);
filter_50Hz_ianeg(5120-11+2) = ia_neg_FFF(5120-11+2);
filter_95Hz_ianeg(20) = ia_neg_FFF(20);
filter_95Hz_ianeg(5120-20+2) = ia_neg_FFF(5120-20+2);
filter_100Hz_ianeg(21) = ia_neg_FFF(21);
filter_100Hz_ianeg(5120-21+2) = ia_neg_FFF(5120-21+2);
filter_105Hz_ianeg(22) = ia_neg_FFF(22);
filter_105Hz_ianeg(5120-22+2) = ia_neg_FFF(5120-22+2);

% figure;
% plot(t,ifft(filter_50Hz_ianeg));
% hold on;
% plot(t,ifft(filter_95Hz_ianeg));
% plot(t,ifft(filter_100Hz_ianeg));
% plot(t,ifft(filter_105Hz_ianeg));
% legend('50Hz comp','95Hz comp','100Hz comp','105Hz comp')
% xlabel('time');
% ylabel('I(A)');


%ib positive sequence
ib_pos_FFT = fft(ib_POS);
ibposfreq = (0:length(ib_pos_FFT)-1)*25600/length(ib_pos_FFT);
% generate empty template for filtered components
filter_50Hz_ibpos = zeros(1,length(ib_pos_FFT));
filter_95Hz_ibpos = zeros(1,length(ib_pos_FFT));
filter_100Hz_ibpos = zeros(1,length(ib_pos_FFT));
filter_105Hz_ibpos = zeros(1,length(ib_pos_FFT));
% assign the desired part of array to the empty template
% FFT is symmetric, two elemets need to exist in the filtered current FFTs
filter_50Hz_ibpos(11) = ib_pos_FFT(11);
filter_50Hz_ibpos(5120-11+2) = ib_pos_FFT(5120-11+2);
filter_95Hz_ibpos(20) = ib_pos_FFT(20);
filter_95Hz_ibpos(5120-20+2) = ib_pos_FFT(5120-20+2);
filter_100Hz_ibpos(21) = ib_pos_FFT(21);
filter_100Hz_ibpos(5120-21+2) = ib_pos_FFT(5120-21+2);
filter_105Hz_ibpos(22) = ib_pos_FFT(22);
filter_105Hz_ibpos(5120-22+2) = ib_pos_FFT(5120-22+2);

% figure;
% plot(t,ifft(filter_50Hz_ibpos));
% hold on;
% plot(t,ifft(filter_95Hz_ibpos));
% plot(t,ifft(filter_100Hz_ibpos));
% plot(t,ifft(filter_105Hz_ibpos));
% legend('50Hz comp','95Hz comp','100Hz comp','105Hz comp')
% xlabel('time');
% ylabel('I(A)');

%ib negative sequence
ib_neg_FFT = fft(ib_NEG);
ibnegfreq = (0:length(ib_neg_FFT)-1)*25600/length(ib_neg_FFT);
% generate empty template for filtered components
filter_50Hz_ibneg = zeros(1,length(ib_neg_FFT));
filter_95Hz_ibneg = zeros(1,length(ib_neg_FFT));
filter_100Hz_ibneg = zeros(1,length(ib_neg_FFT));
filter_105Hz_ibneg = zeros(1,length(ib_neg_FFT));
% assign the desired part of array to the empty template
% FFT is symmetric, two elemets need to exist in the filtered current FFTs
filter_50Hz_ibneg(11) = ib_neg_FFT(11);
filter_50Hz_ibneg(5120-11+2) = ib_neg_FFT(5120-11+2);
filter_95Hz_ibneg(20) = ib_neg_FFT(20);
filter_95Hz_ibneg(5120-20+2) = ib_neg_FFT(5120-20+2);
filter_100Hz_ibneg(21) = ib_neg_FFT(21);
filter_100Hz_ibneg(5120-21+2) = ib_neg_FFT(5120-21+2);
filter_105Hz_ibneg(22) = ib_neg_FFT(22);
filter_105Hz_ibneg(5120-22+2) = ib_neg_FFT(5120-22+2);

% figure;
% plot(t,ifft(filter_50Hz_ibneg));
% hold on;
% plot(t,ifft(filter_95Hz_ibneg));
% plot(t,ifft(filter_100Hz_ibneg));
% plot(t,ifft(filter_105Hz_ibneg));
% legend('50Hz comp','95Hz comp','100Hz comp','105Hz comp')
% xlabel('time');
% ylabel('I(A)');


%ic positive sequence
ic_pos_FFT = fft(ic_POS);
icposfreq = (0:length(ic_pos_FFT)-1)*25600/length(ic_pos_FFT);
% generate empty template for filtered components
filter_50Hz_icpos = zeros(1,length(ic_pos_FFT));
filter_95Hz_icpos = zeros(1,length(ic_pos_FFT));
filter_100Hz_icpos = zeros(1,length(ic_pos_FFT));
filter_105Hz_icpos = zeros(1,length(ic_pos_FFT));
% assign the desired part of array to the empty template
% FFT is symmetric, two elemets need to exist in the filtered current FFTs
filter_50Hz_icpos(11) = ic_pos_FFT(11);
filter_50Hz_icpos(5120-11+2) = ic_pos_FFT(5120-11+2);
filter_95Hz_icpos(20) = ic_pos_FFT(20);
filter_95Hz_icpos(5120-20+2) = ic_pos_FFT(5120-20+2);
filter_100Hz_icpos(21) = ic_pos_FFT(21);
filter_100Hz_icpos(5120-21+2) = ic_pos_FFT(5120-21+2);
filter_105Hz_icpos(22) = ic_pos_FFT(22);
filter_105Hz_icpos(5120-22+2) = ic_pos_FFT(5120-22+2);

% figure;
% plot(t,ifft(filter_50Hz_icpos));
% hold on;
% plot(t,ifft(filter_95Hz_icpos));
% plot(t,ifft(filter_100Hz_icpos));
% plot(t,ifft(filter_105Hz_icpos));
% legend('50Hz comp','95Hz comp','100Hz comp','105Hz comp')
% xlabel('time');
% ylabel('I(A)');


%ic negative sequence
ic_neg_FFT = fft(ic_NEG);
icnegfreq = (0:length(ic_neg_FFT)-1)*25600/length(ic_neg_FFT);
% generate empty template for filtered components
filter_50Hz_icneg = zeros(1,length(ic_neg_FFT));
filter_95Hz_icneg = zeros(1,length(ic_neg_FFT));
filter_100Hz_icneg = zeros(1,length(ic_neg_FFT));
filter_105Hz_icneg = zeros(1,length(ic_neg_FFT));
% assign the desired part of array to the empty template
% FFT is symmetric, two elemets need to exist in the filtered current FFTs
filter_50Hz_icneg(11) = ic_neg_FFT(11);
filter_50Hz_icneg(5120-11+2) = ic_neg_FFT(5120-11+2);
filter_95Hz_icneg(20) = ic_neg_FFT(20);
filter_95Hz_icneg(5120-20+2) = ic_neg_FFT(5120-20+2);
filter_100Hz_icneg(21) = ic_neg_FFT(21);
filter_100Hz_icneg(5120-21+2) = ic_neg_FFT(5120-21+2);
filter_105Hz_icneg(22) = ic_neg_FFT(22);
filter_105Hz_icneg(5120-22+2) = ic_neg_FFT(5120-22+2);
% 
% figure;
% plot(t,ifft(filter_50Hz_icneg));
% hold on;
% plot(t,ifft(filter_95Hz_icneg));
% plot(t,ifft(filter_100Hz_icneg));
% plot(t,ifft(filter_105Hz_icneg));
% legend('50Hz comp','95Hz comp','100Hz comp','105Hz comp')
% xlabel('time');
% ylabel('I(A)');

% calculation pf RMS values
RMS_50Hz_pos=rms(ifft(filter_50Hz_iapos))
RMS_50Hz_neg=rms(ifft(filter_50Hz_ianeg))
RMS_95Hz_pos= rms(ifft(filter_95Hz_iapos))
RMS_95Hz_neg=rms(ifft(filter_95Hz_ianeg))
RMS_100Hz_pos=rms(ifft(filter_100Hz_iapos))
RMS_100Hz_neg=rms(ifft(filter_100Hz_ianeg))
RMS_105Hz_pos=rms(ifft(filter_105Hz_iapos))
RMS_105Hz_neg=rms(ifft(filter_105Hz_ianeg))


