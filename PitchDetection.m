%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Speech prject I
%% Start/End Voiced/Unvoiced Speech Detection on test wave form H.1.wav
%% Pitch detection on test waveform H.1.wav
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Intialization and setting of training data
clc;
clear all;
%% Thresholds got from VUS training set
ITR = -30;
ITU =-20;
IZCT = 35;
window_time =30*1e-3; % time 30 msec
overlap_time =10*1e-3;     % time 10 msec

% Testing on Test sample H.1.wav  
[y_test,Fs] = wavread('H.1.wav');
window_length =window_time*Fs;
overlap_length = overlap_time*Fs;

% Filter speech to eliminate noise

[b1,a1] =fir1(301,(60*2/1000),'high');
y_filt1 = filter(b1,a1,y_test);
y_frames1 = buffer(y_filt1,window_length,(window_length-overlap_length),'nodelay');
win1= hamming(window_length);

% Calculate Ennergy in Signal
for j=1:window_length
    Test_energy(j,:) = (win1(j,1)*y_frames1(j,:)).^2;
end
square_energy_test = sum(Test_energy,1);
log_energy_test = (10*log10(square_energy_test))-max(10*log10(square_energy_test)); 
 
% Calculate ZCR 
 test_zc1 = sign(y_frames1);
 [r1,c1] = size(test_zc1);
 test_zero1 = zeros(c1,1);
 
 for col =1:c1
     for row=1:r1-1
         c1 = test_zc1(row,col) - test_zc1(row+1,col);
         test_zero1(col,1) = test_zero1(col,1) +abs(c1);
     end
 end
 zero_crossing1 = (overlap_length)/(2*window_length) .*test_zero1;

 
 
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
         disp(v)
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
         
     end
  end
  
  
  
  
  
  %% Checking of voiced and unvoiced signals
 
  [r_1,c_1] =size(log_energy_test);
 i=1;j=1;m=1;
 X_test = zeros(l_limit,1);
 X_test(:,1) =2;
 X_test(1:B1,1)=3;
 X_test(E1:l_limit,1)=3;

for cc = B1:E1 
 
if(log_energy_test(1,cc)>=-30 || zero_crossing1(cc,1)<35)
         X_test(cc,1)=1;
         
elseif(log_energy_test(1,cc)<-40 &&zero_crossing1(cc,1)<35)
X_test(cc,1)=3;
else
     X_test(cc,1)=2;
          
     end
 end
         
  %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%data got from praat for voiced and unvoiced frames
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

X_time = [0.14,0.32; 0.36, 0.51; 0.71, 1.03; 1.21,1.70; 1.85,2.10; 2.27,2.57];
X_voiced_frame = (floor(X_time.*11000)./110)+1;
X_utime = [0.12,0.14; 0.32,0.35; 0.51, 0.70; 1.03, 1.21; 1.70, 1.85; 2.10, 2.27; 2.57,2.81];
X_unvo_frame = (floor(X_utime.*11000)./110)+1;

X_train=zeros(l_limit,1);
X_train(:,1) = 3;

for j = 1:6
          X_train(X_voiced_frame(j,1):X_voiced_frame(j,2),1)=1;
          
     end     
for j = 1:7
          X_train(X_unvo_frame(j,1):X_unvo_frame(j,2),1)=2;
end
     
 X_final_answer = [X_train X_test];
          ccr =confusionmat(X_train,X_test);

         ccr_percent = sum(diag(ccr))./sum(sum(ccr))*100     
 
  





  
  %%%%%%%%%%%%% THE PITCH
  %%%%%%%%%%%%% DETECTION%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  i_p=0;t=1;i_p1=0;i_p2=0;
  p=1;len=0;
  X_test1 =X_test;
  while(  p<length(X_test1))
       % Pull out all speech segments that are predicted as voice by thresholding      
      if (X_test1(p,1)==1 )
          len = (p-1)*110;
          X_pitch_voiced(i_p+1:overlap_length+i_p,1)=y_test(len:((len+overlap_length)-1),1);
          i_p = i_p+110;
      end
      
      % Pull out speech segments that are got from praat
      if(X_train(p,1)==1)
         X_pitch_voiced1(i_p1+1:110+i_p1,1)=y_test(len:((len+overlap_length)-1),1);
         i_p1 = i_p1+110;
      end
      
      %Pull only those speech segments that are correctly classifies as
      %voiced
      if(X_test1(p,1)==1 && X_train(p,1)==1)
         X_pitch_voiced2(i_p2+1:110+i_p2,1)=y_test(len:((len+110)-1),1);
          i_p2 = i_p2+110;
      end
      p=p+1;
  end
  
  % Calculate pitch for all the three speech segments
   y_pitch_test=cell(3,1);
   pitch_out_final = cell(3,1);
   y_pitch_test(1,1) ={X_pitch_voiced};
   y_pitch_test(2,1) ={wavread('H_1_voiced.wav')};%{X_pitch_voiced1};
   y_pitch_test(3,1) ={X_pitch_voiced2};%{ wavread('H_1_voiced.wav')};
  
   % looping to calculate pitch for all three segments
   for pi =1:3
      clear X_pitch_filt;
      clear v_pitch;
      clear x_temp;
      clear final_pitch;
         [b2,a2] = fir1(301,[60 900]*2/Fs,'bandpass'); 
         X_pitch_filt = filter(b2,a2,y_pitch_test{pi,1});
         v_pitch =buffer(X_pitch_filt,(window_length),(window_length-overlap_length),'nodelay');
     
      [r,c] = size(v_pitch);
      x_temp=zeros(r,c);
         for d=1:c
                        a_p = max(abs(v_pitch(1:100,d)));
                        b_p = max(abs(v_pitch(230:window_length,d)));
                        clip = 0.68* min(a_p,b_p);
                            for t1=1:330
                               if (v_pitch(t1,d)>clip)
                                  x_temp(t1,d)=v_pitch(t1,d)-clip;
                               elseif(v_pitch(t1,d)<-clip)
                                    x_temp(t1,d) = v_pitch(t1,d)+clip;
                               else
                                     x_temp(t1,d)=0;
                               end
                             end
                         r = xcorr(x_temp(:,d),Fs/50);
                         ms_2=floor(Fs/500);                 % 2ms
                         ms_20=floor(Fs/50);                 % 20ms
                                             
                         r = r(floor(length(r)/2):end);
                         [maxi,idx]=max(r(ms_2:ms_20));
                      
                          pitch_f = Fs/(ms_2+idx-1);
                          final_pitch(d,1)=pitch_f;
                        
 end  
 final_pitch1 = medfilt1(final_pitch,4);
 
 pitch_out_final(pi,1) = {final_pitch1};
   end
   
   figure(1)
    
     p_l1 = length(pitch_out_final{1,1});
    T = (0:1:p_l1-1)/100;
    subplot(3,1,1);
    plot(T,pitch_out_final{1,1},'*');
    title('Pitch of  all frames classified as voice by Thresholding');
    xlabel('Time in sec');
    ylabel('Pitch Frequency');
   
    p_l2 = length(pitch_out_final{2,1});
    T2 = (0:1:p_l2-1)/100;
    subplot(3,1,2)
    plot(T2,pitch_out_final{2,1},'*');
   title('Pitch of voiced segments got from praat');
    
    xlabel('Time in sec');
    ylabel('Pitch Frequency');
    
    subplot(3,1,3)
    
    p_l3 = length(pitch_out_final{3,1});
    T3 = (0:1:p_l3-1)/100;
    plot(T3,pitch_out_final{3,1},'*');
    title('Pitch of frames correctly classified as voice by Thresholding');
     xlabel('Time in sec');
    ylabel('Pitch Frequency');
    
    %%%
    figure(2)
p_l2 = length(pitch_out_final{2,1});
    T2 = (0:1:p_l2-1)/100;
    
    plot(T2,pitch_out_final{2,1},'*','color','b');
   title('Pitch of voiced segments got from praat');
    
    xlabel('Time in sec');
    ylabel('Pitch Frequency');
    hold on;
    
    
    p_l3 = length(pitch_out_final{3,1});
    T3 = (0:1:p_l3-1)/100;
    plot(T3,pitch_out_final{3,1},'*','color','g');
    title('Pitch of frames correctly classified as voice by Thresholding');
     xlabel('Time in sec');
    ylabel('Pitch Frequency');
    legend('Original pitch signal','Pitch of correctly detected voice frames through thresholding','Location','southeast');




v_pitch =pitch_out_final{3,1};
v_pitch_main = pitch_out_final{2,1};
lk = length(v_pitch_main);
dc = zeros(lk,1);
dc(1:length(v_pitch),1) = v_pitch;
diff = (v_pitch_main-dc);
mean(diff);
disp('Average difference in the pitch frequencies is')
disp(mean(diff))
   