%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                       Speech prject I
%% Start/End Voiced/Unvoiced Speech Detection on test wave form H.1.wav
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Intialization and setting of training data
clc;
clear all;
train_signal = cell(7,1);
log_energy = cell(7,1);
zero_crossing = cell(7,1);
window_time =30*1e-3; % time 30 msec
overlap_time =10*1e-3;     % time 10 msec
% Read and store training data set Samples
%[train_signal{2,1},~] = wavread('H.3.wav');
%[train_signal{3,1},~] =  wavread('H.4.wav');
%[train_signal{4,1},~] = wavread('D.6.2.wav');
%[train_signal{5,1},~] = wavread('D.6.3.wav');
%[train_signal{6,1},~] = wavread('D.6.4.wav');
%[train_signal{7,1},~] = wavread('D.6.5.wav');

window_length =window_time*Fs;
overlap_length = overlap_time*Fs;

% for all audio files filter for noise removal and calculate log energy and ZCR
for i_t =1:7
        clear b;
        clear a;
        clear y_filt;
        clear y_signal;
        clear y_frames;
        clear Result_energy;
        clear square_energy;
        clear y_sig_zc;clear result_zero;
        
        y_signal = train_signal{i_t,1};
        [b,a] =fir1(101,(60*(2/Fs)),'high');  % 60Hz cutoff high pass filter
        speech_filter = filter(b,a,y_signal);
        speech_frames = buffer(speech_filter,window_length,(window_length-overlap_length),'nodelay'); % overlap of 60 percent
        win= hamming(window_length);

  
          % Short time log energy
                for i=1:size(speech_frames,2)
                 square_energy(1,i) = sum((win(:,1).*speech_frames(:,i)).^2,1);
                end
                log_energy(i_t,1) = {(10*log10(square_energy))- max(10*log10(square_energy))};
               
          % Zero crossing
                speech_zc = sign(speech_frames);
                [r,c] = size(speech_zc);
                result_zero = zeros(c,1);
                for col =1:c
                        for row=1:r-1
                            c = speech_zc(row,col) - speech_zc(row+1,col);
                                    result_zero(col,1) = result_zero(col,1) +abs(c);
                        end
                end
               zero_crossing(i_t,1) = {(overlap_length)/(2*window_length) .*result_zero};
end

X_train_log = cell(7,1);
X_train_ZCR = cell(7,1);

%% Manual collection of time locations of voiced speech present in training data samples...
%% for determining thresholds to determine voiced and unvoiced regions.
X_train_log (1,1)= {[0.08,0.14;0.51,0.66;0.79,0.93;1.36,1.51;1.98,2.11;2.30,2.44;2.48,2.53]};
X_train_log (2,1) = {[0.08,0.12; 0.20,0.33;0.50,0.54;0.75,0.85;1.24,1.28;1.41,1.56;2.14,2.19;2.22,2.3]};
X_train_log(3,1) = {[0.01,0.10;0.78,0.91;0.98,1.06;1.55,1.58;2.04,2.07;2.27,2.55]};
X_train_log(4,1) = {[0.09,0.12;0.43,0.56;0.65,0.87;1.0,1.48;1.95,2.14;2.32,2.41;2.61,2.63;2.68,2.94]};
X_train_log(5,1) = {[0.09,0.17;0.70,0.76;1.23,1.30;1.51,1.70;1.88,2.02;2.13,2.29;2.56,2.83]};
X_train_log(6,1) = {[0.10,0.16;0.62,0.75;1.48,1.59;1.80,1.93;2.24,2.32;2.53,2.77]};
X_train_log(7,1) = {[0.02,0.1;0.41,0.51;1.03,1.27;1.79,1.93;2.37,2.67]};

for i =1:7
    X_train_ZCR(i,1)= {floor(X_train_log{i,1}.*100)+1};
end

i=0;i1=0;ii=10;

for i2 =1:7
    
clear v_zcr;
clear v_en;
clear v;
v = X_train_ZCR{i2,1};
v_zcr = zero_crossing{i2,1};
v_en = log_energy{i2,1};
v_en_l = length(v_en);

[r,c] = size(v);

%% Extracting silenced /voiced/unvoiced frames of energy and ZCr for histogram plotting
for r1=1:r
    
    %%%Extracting values of silence from the start
    if r1==1
         f_1=length(1:v(r1,1));
         f_1 = f_1-1;
       if (f_1>0)
        x_silence_en(1,ii+1:f_1+ii) = v_en(1,1:f_1);
        x_silence_zcr(1,ii+1:f_1+ii) = v_zcr(1:f_1,1);
         ii=ii+f_1;
       end
    end
    
    %%% Extracting values of silence from the end
    if r1==r
        
         f_1=length(v(r1,2):v_en_l);
         f_1 = f_1-1;
         disp(f_1)
       if (f_1>0)
        x_silence_en(1,ii+1:f_1+ii) = v_en(1,(v(r1,2)+1):v_en_l);
        x_silence_zcr(1,ii+1:f_1+ii) = v_zcr(v(r1,2)+1:v_en_l,1);
         ii=ii+f_1;
       end
    end
         f = length(v(r1,1):v(r1,2));
        x_unvoiced(1,i+1:f+i)=v_en(1,v(r1,1):v(r1,2));
        x_unvoiced_zcr(1,i+1:f+i) = v_zcr(v(r1,1):v(r1,2),1);
        i=i+f;
end

for r1=1:r-1 
         f1 =length(v(r1,2)+1:v(r1+1,1)-1);
         x_voiced(1,i1+1:f1+i1)=v_en(1,v(r1,2)+1:v(r1+1,1)-1);
         x_voiced_zcr(1,i1+1:f1+i1) = v_zcr(v(r1,2)+1:v(r1+1,1)-1,1);
        i1=i1+f1;
end



end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting of histograms and scatter plot based on the training data sample set
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
subplot(3,1,1);
hist(x_unvoiced);
set(gca,'xlim',[-60,0]);
xlabel('dB Energy');
title('Histogram of log energy of unoviced');

subplot(3,1,2);
hist(x_voiced);
xlabel('dB Energy');
set(gca,'xlim',[-60,0]);
title('Histogram of log energy of voiced ');
subplot(3,1,3);

hist(x_silence_en);
set(gca,'xlim',[-60,0]);
xlabel('dB Energy');
title('Histogram of log energy of silence ');


figure(2)
subplot(3,1,1)
hist(x_unvoiced_zcr);
xlabel('ZCR');
title('Histogram of zero crossing of unvoiced');
subplot(3,1,2)
hist(x_voiced_zcr);
xlabel('ZCR');
title('Histogram of zero crossing of voiced');
subplot(3,1,3);
hist(x_silence_zcr);
xlabel('ZCR');
title('Histogram of zero crossing of silence ');

figure(3);
X_trial = [x_silence_zcr,x_unvoiced_zcr,x_voiced_zcr];
Y_trial = [x_silence_en,x_unvoiced,x_voiced];
Y_label=zeros(1,1913);
Y_label(1,1:97)=3;
Y_label(1,98:751)=2;
Y_label(1,752:1913)=1;
group =Y_label;
clr =['r','g','b'];
sym=['o','+','*'];

gscatter(X_trial,Y_trial,group,clr,sym);
xlabel('ZCR');
ylabel('log energy');
title('log-energy vs zero crossing')
legend('voiced','unvoiced','silence','Location','southeast');



%% Calculating the mean of log energy of silence and voiced and unvoiced

voiced_mean = mean(x_voiced);
uv_mean = mean(x_unvoiced);
 s_mean = mean(x_silence_en);
eavg = voiced_mean;
esig = std(x_voiced);
zcavg = mean(x_silence_zcr);
zcsig = std(x_silence_zcr);
%% Setting Threshods based on training data set
                IZCT = max(35,(zcavg+3*zcsig));
                ITU = -20;
                ITR = -30;
               
               
%% Testing of test sample to get voiced and unvoiced  and silence frames

    %[y_test,Fs_test] = wavread('H.1.wav');
    [b1,a1] =fir1(101,(60*2/1000),'high');
    y_filt1 = filter(b1,a1,y_test);
    % to check signal with noise uncomment the two lines below and run code
%       y_snr = awgn(y_filt1,20);
%       y_filt1 = y_snr;
    speech_frames1 = buffer(y_filt1,window_length,(window_length-overlap_length),'nodelay');
    win1= hamming(window_length);

   
    for j=1:330
        Test_energy(j,:) = (win1(j,1)*speech_frames1(j,:)).^2;
    end
 square_energy_test = sum(Test_energy,1);
 log_energy_test = (10*log10(square_energy_test))-max(10*log10(square_energy_test)); 
 
 
 test_zc1 = sign(speech_frames1);
 [r1,c1] = size(test_zc1);
 test_zero1 = zeros(c1,1);
 
 for col =1:c1
     for row=1:r1-1
         c1 = test_zc1(row,col) - test_zc1(row+1,col);
         test_zero1(col,1) = test_zero1(col,1) +abs(c1);
     end
 end
 zero_crossing1 = (110)/(2*330) .*test_zero1;
 
%% Start and End point Detection
l_limit = length(log_energy_test);
 B1=0;E1=0;
 j=1;
 for kk =1:l_limit
     if (log_energy_test(1,kk)>=ITR)
         v =kk;
         
         while(v>0)
             if(zero_crossing1(v,1)>IZCT)
                B1(j,1)=v-1;
                j=j+1;
                break;
             else
                 v=v-1;
             end
         end
         
         break;
     else
         
     end
 end
 j=1;
  for kk =0:l_limit
      if (log_energy_test(1,(l_limit-kk))>=ITR)
         v =(l_limit-kk);
         
         while(v<l_limit)
             if(zero_crossing1(v,1)>IZCT)
                E1(j,1)=v;
                j=j+1;
                break;
             else
                 v=v+1;
             end
         end
         
         break;
     else
        % disp(l_limit-kk)
     end
  end
 %% conditions for when noise is added to the original signal
%   if E1 ==0
%       E1=l_limit;
%   end
%   if B1==0
%       B1=1;
%   end
 
  
  
  
  %% Checking of voiced and un voiced signals 
  [r_1,c_1] =size(log_energy_test);
 i=1;j=1;m=1;
 X_test = zeros(l_limit,1);
 X_test(:,1) =2;
 X_test(1:B1,1)=3;
 X_test(E1:l_limit,1)=3;

for cc = B1:E1 
 
if (log_energy_test(1,cc)>=ITR && zero_crossing1(cc,1)<IZCT)
         X_test(cc,1)=1;
         
elseif(log_energy_test(1,cc)<-40 && zero_crossing1(cc,1)<30)
X_test(cc,1)=3;
else
     X_test(cc,1)=2;


     end
 end
         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  This section compares the result with algorithm based result with 
%  manually extracted results from praat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

X_time = [0.14,0.32; 0.35, 0.51; 0.71, 1.03; 1.21,1.70; 1.85,2.10; 2.27,2.57];
X_voiced_frame = (floor(X_time.*11000)./110)+1;
X_utime = [0.12,0.14; 0.32,0.35; 0.51, 0.71; 1.03, 1.21; 1.70, 1.85; 2.10, 2.27; 2.57,2.81];
X_unvo_frame = (floor(X_utime.*11000)./110)+1;

X_train=zeros(l_limit,1);
X_train(:,1) = 3;

for j = 1:6
          X_train(X_voiced_frame(j,1):X_voiced_frame(j,2),1)=1;
          
     end     
for j = 1:7
          X_train(X_unvo_frame(j,1):X_unvo_frame(j,2),1)=2;
end
     
     
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  X_final_answer = [X_train X_test];
  ccr =confusionmat(X_train,X_test)
  ccr_percent = sum(diag(ccr))./sum(sum(ccr))*100
  %%%%%%%%%%%%%%%% THE END%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% display waveform with start and end points and voiced and unvoiced frames
  % for test data
  
 s_f =length(log_energy_test);
B1_f(1,1:s_f) = B1;
E1_f(1,1:s_f) = E1;

figure(4)
l=1:s_f;
hold on;
plot(l,log_energy_test,'color','r');
hold on;
plot(l,(10*log10(X_test)),'color','b');
hold on;
plot(l,(10*log10(X_train.*10)),'color','g');
title('Log Energy of Test Signal H.1.wav');
legend('Energy','Classified Result','Actual Result','Location','southeast')
BB =(B1-1)*110/Fs;
EE = (E1-1)*110/Fs;
disp('Start and End time for speech signal are')
disp(BB)
disp(EE)

figure(5)
N = length(y_test);
T = (0:1:N-1)/Fs;
plot(T,y_test);
 line([BB,BB],get(gca,'ylim'),'color','r');
 line([EE,EE],get(gca,'ylim'),'color','b'); 
 legend('END','START');
  title('Speech Signal with Start and End points marked');

  
  
