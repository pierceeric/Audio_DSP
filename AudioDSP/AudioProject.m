%
% Eric Pierce
% University of North Dakota
% Digital Signal Processing
% 
%
% This program demonstrates LP Filter, HP Filter, Band-Pass Filter,
% Echo, Down-Sampling and Up-Sampling of input audio files.
% 
% After applying each sound effect, the audio file is played,
% then a FFT is performed on the signal and a plot is displayed,
% showing a graphical representation of the audio file just heard.
%

clc;
clear;
close all;

% Set path where files will be created
cd 'C:\Users\epier\Desktop\DSP_MATLAB\AudioProject\'

% Place music file into array of pressure values between -1 and 1
% fs is sampling frequency from file creation, hithouse - 44100, pacman 11025
% Uncomment to select audio track to use
[audioSample,sampleFrequency] = audioread("HithouseSample.wav");
%[audioSample,sampleFrequency] = audioread("pacman_beginning.wav");
%[audioSample,sampleFrequency] = audioread("test.wav");

% Calulate duration of each sample
dt = 1/sampleFrequency;

% Find number of samples
numberOfSamples = length(audioSample);

% Create a time vector from 0 to time of final sample
t = 0:dt:((numberOfSamples-1)*dt);  

% Calculate the duration of the Audio Clip in seconds
audioDuration = numberOfSamples/sampleFrequency; 

% Play audio vector over speakers
sound(audioSample,sampleFrequency,16);
pause(audioDuration) % pause duration of sample to listen

% Plot Signal in Time Domain
t = linspace(0,numberOfSamples*dt,numberOfSamples);
subplot(5,1,1);
plot(t,audioSample);
title('Time-Domain Audio Signal');
xlabel('time (s)');
ylabel('Amplitude');

% Plot Signal in Frequency Domain
fftOfAudioSample = fft(audioSample);      % DFT of signal
subplot(5,1,2)
plot(real(fftOfAudioSample));
xlim([20 20000]);     % limit to human hearing range
ylim([0 1500])        % arbitrary amplitude limit
title('Freq-Domain Audio Signal')
xlabel('freq (Hz)');
ylabel('Amplitude');

% Wait for user to press a button
disp('Press a button to apply a Low-Pass Filter')
w = waitforbuttonpress; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Low-Pass Filter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cutOffFrequency = 100;    % Set the cutoff frequency in Hz
filterOrder = 100;        % Set the filter order -this had to be
% much higher than I anticipated to cut out highs 
% 50 for PacMan and 100 for Peter Slaghuis
lowpassFilter = designfilt('lowpassfir', ...
    'CutoffFrequency', cutOffFrequency, ...
    'FilterOrder', filterOrder, ...
    'SampleRate', sampleFrequency);

% Apply the low-pass filter to the input signal 
outputSignal = filter(lowpassFilter, audioSample);

% Write the filtered signal to an .wav file
fileName = 'audio_lpf.wav';
audiowrite(fileName, outputSignal, sampleFrequency);

% Play audio vector over speakers 
sound(outputSignal,sampleFrequency,16) 
pause(audioDuration) % pause duration of sample to listen

% Perform FFT
x = fft(outputSignal);      % DFT of filtered signal

% Graph Filtered Frequency Domain 
subplot(5,1,3)
plot(real(x));
xlim([20 20000]);     % limit to human hearing range
ylim([0 1500])        % arbitrary amplitude limit
title('Freq-Domain LPF Audio Signal')
xlabel('freq (Hz)');
ylabel('Amplitude');

disp('Press a button to apply a Band-Pass Filter')
w = waitforbuttonpress; % wait for user to press a button

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Band-Pass Filter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bandpassFilter = designfilt('bandpassfir', ...
    'FilterOrder', 200, ...
    'CutoffFrequency1', 300, ...
    'CutoffFrequency2', 1000, ...
    'SampleRate',sampleFrequency);

% Apply the band-pass filter to the input signal 
outputSignal = filter(bandpassFilter, audioSample);

% Write the filtered signal to an .wav file
fileName = 'audio_bp.wav';
audiowrite(fileName, outputSignal, sampleFrequency);

% Play audio vector over speakers 
sound(outputSignal,sampleFrequency,16)
pause(audioDuration) % pause duration of sample to listen

% Perform FFT
x = fft(outputSignal);      % DFT of filtered signal

% Graph Filtered Frequency Domain 
subplot(5,1,4)
plot(real(x));
xlim([20 20000]);     % limit to human hearing range
ylim([0 1500])        % arbitrary amplitude limit
title('Freq-Domain Band-Pass Audio Signal')
xlabel('freq (Hz)');
ylabel('Amplitude');

% Wait for the user to press a button
disp('Press a button to apply a High-Pass Filter')
w = waitforbuttonpress; 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% High-Pass Filter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

highpassFilter = designfilt('highpassfir', ...
    'FilterOrder', 50, ...
    'CutoffFrequency',2000, ...
    'SampleRate',sampleFrequency);

% Apply the band-pass filter to the input signal 
outputSignal = filter(highpassFilter, audioSample);

% Write the filtered signal to an .wav file
fileName = 'audio_hpf.wav';
audiowrite(fileName, outputSignal, sampleFrequency);

% Play audio vector over speakers 
sound(outputSignal,sampleFrequency,16)
pause(audioDuration) % pause duration of sample to listen

% Perform FFT
x = fft(outputSignal);      % DFT of filtered signal

% Graph Filtered Frequency Domain 
subplot(5,1,5)
plot(real(x));
xlim([20 20000]);     % limit to human hearing range
ylim([0 1500])        % arbitrary amplitude limit
title('Freq-Domain HPF Audio Signal')
xlabel('freq (Hz)');
ylabel('Amplitude');

% Wait for the user to press a button
disp('Press a button to apply an Echo Effect')
w = waitforbuttonpress; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Echo
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sounds cool with PacMan

% Parameters to play with
delayInSeconds = 0.15 ; % Set the delay time in seconds
decayFactor = 0.5;    % Reduces volume of the echo

% Calculate number of preceding zeros for delay
numberOfZeros = round(delayInSeconds * sampleFrequency); 

% Create a column array of zeros
columnOfZeros = zeros(numberOfZeros, 1); 

% Create a delayed version of the input signal
delayedSignal = [columnOfZeros; audioSample];

% Add zeros to end of audioSample to equal size of delayed signal
audioSampleWithZeros = audioSample;
audioSampleWithZeros(numel(delayedSignal))= 0;

% Mix the original and decayed delayed signals
outputSignal = audioSampleWithZeros + (decayFactor * delayedSignal);

% Write the Echo Signal to an audio file
fileName = 'audio_echo.wav';
audiowrite(fileName, outputSignal, sampleFrequency);

% Play audio vector over speakers 
sound(outputSignal,sampleFrequency,16)
pause(audioDuration) % pause duration of sample to listen

% Plot Signal in Time Domain
figure;
% Find number of samples
numberOfSamples = length(outputSignal);

t = linspace(0,numberOfSamples*dt,numberOfSamples);
figure;
subplot(2,1,1);
plot(t,outputSignal);
title('Time-Domain Echo Signal');
xlabel('time (s)');
ylabel('Amplitude');

% Perform FFT
x = fft(outputSignal);   % DFT of signal

% Graph Filtered Frequency Domain 
subplot(2,1,2)
plot(real(x));
xlim([20 20000]);     % limit to human hearing range
ylim([0 1500])        % arbitrary amplitude limit
title('Freq-Domain Echo Signal')
xlabel('freq (Hz)');
ylabel('Amplitude');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Down Sampling by 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Wait for the user to press a button
disp('Press a button to apply Down Sampling by 2')
w = waitforbuttonpress; 

% Set the desired output sample rate
outputSampleRate = round(sampleFrequency / 2);

% Downsample the audio signal using the resample function
outputSignal = resample(audioSample, outputSampleRate, sampleFrequency);

% Write the downsampled signal to the output audio file
audiowrite('audio_downsampleby002.wav', outputSignal, outputSampleRate);

% Play audio vector over speakers 
sound(outputSignal,outputSampleRate,16)
pause(audioDuration) % pause duration of sample to listen

% Perform FFT
x = fft(outputSignal);   % DFT of signal

% Graph Filtered Frequency Domain 
figure;

subplot(7,1,1)
plot(real(x));
xlim([20 20000]);     % limit to human hearing range
ylim([0 1500])        % arbitrary amplitude limit
title('Freq-Domain Downsample by 2 Audio Signal')
xlabel('freq (Hz)');
ylabel('Amplitude');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Down Sampling by 4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Wait for the user to press a button
disp('Press a button to apply Down Sampling by 4')
w = waitforbuttonpress; 

% Set the desired output sample rate
outputSampleRate = round(sampleFrequency / 4);

% Downsample the audio signal using the resample function
outputSignal = resample(audioSample, outputSampleRate, sampleFrequency);

% Write the downsampled signal to the output audio file
 audiowrite('audio_downsampleby004.wav', outputSignal, outputSampleRate);

% Play audio vector over speakers 
if outputSampleRate > 1000
    sound(outputSignal,outputSampleRate,16)
    pause(audioDuration) % pause duration of sample to listen
end

% Perform FFT
x = fft(outputSignal);   % DFT of signal

% Graph Filtered Frequency Domain 
subplot(7,1,2)
plot(real(x));
xlim([20 20000]);     % limit to human hearing range
ylim([0 1500])        % arbitrary amplitude limit
title('Freq-Domain Downsample by 4 Audio Signal')
xlabel('freq (Hz)');
ylabel('Amplitude');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Down Sampling by 8
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Wait for the user to press a button
disp('Press a button to apply Down Sampling by 8')
w = waitforbuttonpress; 

% Set the desired output sample rate
outputSampleRate = round(sampleFrequency / 8);

% Downsample the audio signal using the resample function
outputSignal = resample(audioSample, outputSampleRate, sampleFrequency);

% Write the downsampled signal to the output audio file
audiowrite('audio_downsampleby008.wav', outputSignal, outputSampleRate);
 
% Play audio vector over speakers 
if outputSampleRate > 1000
    sound(outputSignal,outputSampleRate,16)
    pause(audioDuration) % pause duration of sample to listen
end

% Perform FFT
x = fft(outputSignal);   % DFT of signal

% Graph Filtered Frequency Domain 
subplot(7,1,3)
plot(real(x));
xlim([20 20000]);     % limit to human hearing range
ylim([0 1500])        % arbitrary amplitude limit
title('Freq-Domain Downsample by 8 Audio Signal')
xlabel('freq (Hz)');
ylabel('Amplitude');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Down Sampling by 16
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Wait for the user to press a button
disp('Press a button to apply Down Sampling by 16')
w = waitforbuttonpress; 

% Set the desired output sample rate
outputSampleRate = round(sampleFrequency / 16);

% Downsample the audio signal using the resample function
outputSignal = resample(audioSample, outputSampleRate, sampleFrequency);

% Write the downsampled signal to the output audio file
audiowrite('audio_downsampleby016.wav', outputSignal, outputSampleRate);
 
% Play audio vector over speakers 
if outputSampleRate > 1000
    sound(outputSignal,outputSampleRate,16)
    pause(audioDuration) % pause duration of sample to listen
end

% Perform FFT
x = fft(outputSignal);   % DFT of signal

% Graph Filtered Frequency Domain 
subplot(7,1,4)
plot(real(x));
xlim([20 20000]);     % limit to human hearing range
ylim([0 1500])        % arbitrary amplitude limit
title('Freq-Domain Downsample by 16 Audio Signal')
xlabel('freq (Hz)');
ylabel('Amplitude');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Down Sampling by 32
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% High frequencies disappear since sample rate is dropping

% Wait for the user to press a button
disp('Press a button to apply Down Sampling by 32')
w = waitforbuttonpress; 

% Set the desired output sample rate
outputSampleRate = round(sampleFrequency / 32);

% Downsample the audio signal using the resample function
outputSignal = resample(audioSample, outputSampleRate, sampleFrequency);

% Write the downsampled signal to the output audio file
audiowrite('audio_downsampleby032.wav', outputSignal, outputSampleRate);
 
% Play audio vector over speakers 
if outputSampleRate > 1000
    sound(outputSignal,outputSampleRate,16)
    pause(audioDuration) % pause duration of sample to listen
end

% Perform FFT
x = fft(outputSignal);   % DFT of signal

% Graph Filtered Frequency Domain 
subplot(7,1,5)
plot(real(x));
xlim([20 20000]);     % limit to human hearing range
ylim([0 1500])        % arbitrary amplitude limit
title('Freq-Domain Downsample by 32 Audio Signal')
xlabel('freq (Hz)');
ylabel('Amplitude');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Down Sampling by 64
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Wait for the user to press a button
disp('Press a button to apply Down Sampling by 64')
w = waitforbuttonpress; 

% Set the desired output sample rate
outputSampleRate = round(sampleFrequency / 64);

% Downsample the audio signal using the resample function
outputSignal = resample(audioSample, outputSampleRate, sampleFrequency);

% Write the downsampled signal to the output audio file
audiowrite('audio_downsampleby064.wav', outputSignal, outputSampleRate);

% minimum Sample Rate for playing sound in MATLAB is 1000 Hz.
% the audio file can still be created and played in another player.
% Play audio vector over speakers 
if outputSampleRate > 1000
    sound(outputSignal,outputSampleRate,16)
    pause(audioDuration) % pause duration of sample to listen
end

% Perform FFT
x = fft(outputSignal);   % DFT of signal

% Graph Filtered Frequency Domain 
subplot(7,1,6)
plot(real(x));
xlim([20 20000]);     % limit to human hearing range
ylim([0 1500])        % arbitrary amplitude limit
title('Freq-Domain Downsample by 64 Audio Signal')
xlabel('freq (Hz)');
ylabel('Amplitude');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Down Sampling by 128
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Wait for the user to press a button
disp('Press a button to apply Down Sampling by 128')
w = waitforbuttonpress; 

% Set the desired output sample rate
outputSampleRate = round(sampleFrequency / 128);

% Downsample the audio signal using the resample function
outputSignal = resample(audioSample, outputSampleRate, sampleFrequency);

% Write the downsampled signal to the output audio file
audiowrite('audio_downsampleby128.wav', outputSignal, outputSampleRate);

% minimum Sample Rate in MATLAB is 1000 Hz.
% Play audio vector over speakers 
if outputSampleRate > 1000
    sound(outputSignal,outputSampleRate,16)
    pause(audioDuration) % pause duration of sample to listen
end

% Perform FFT
x = fft(outputSignal);   % DFT of signal

% Graph Filtered Frequency Domain 
subplot(7,1,7)
plot(real(x));
xlim([20 20000]);     % limit to human hearing range
ylim([0 1500])        % arbitrary amplitude limit
title('Freq-Domain Downsample by 128 Audio Signal')
xlabel('freq (Hz)');
ylabel('Amplitude');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Up Sample by 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Wait for the user to press a button
disp('Press a button to Upsample by 2')
w = waitforbuttonpress; 

% Set the desired output sample rate
outputSampleRate = round(sampleFrequency *2 );

% Upsample the audio signal
outputSignal = resample(audioSample, outputSampleRate, sampleFrequency);

% Write the downsampled signal to the output audio file
audiowrite('audio_upsampleby002.wav', outputSignal, outputSampleRate);

% Play audio vector over speakers 
if outputSampleRate > 1000
    sound(outputSignal,outputSampleRate,16)
    pause(audioDuration) % pause duration of sample to listen
end

% Perform FFT
x = fft(outputSignal);   % DFT of signal

figure;

% Graph Filtered Frequency Domain 
subplot(4,1,1)
plot(real(x));
xlim([20 20000]);     % limit to human hearing range
ylim([0 1500])        % arbitrary amplitude limit
title('Freq-Domain Upsample by 2')
xlabel('freq (Hz)');
ylabel('Amplitude');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Up Sample by 4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Wait for the user to press a button
disp('Press a button to apply Upsample by 4')
w = waitforbuttonpress; 

% Set the desired output sample rate
outputSampleRate = round(sampleFrequency *4 );

% Upsample the audio signal
outputSignal = resample(audioSample, outputSampleRate, sampleFrequency);

% Write the downsampled signal to the output audio file
audiowrite('audio_upsampleby004.wav', outputSignal, outputSampleRate);

% Play audio vector over speakers 
if ((outputSampleRate > 1000)&&(outputSampleRate <256000))
    sound(outputSignal,outputSampleRate,16)
    pause(audioDuration) % pause duration of sample to listen
end

% Perform FFT
x = fft(outputSignal);   % DFT of signal

% Graph Filtered Frequency Domain 
subplot(4,1,2)
plot(real(x));
xlim([20 20000]);     % limit to human hearing range
ylim([0 1500])        % arbitrary amplitude limit
title('Freq-Domain Upsample by 4')
xlabel('freq (Hz)');
ylabel('Amplitude');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Up Sample by 8
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Wait for the user to press a button
disp('Press a button to Upsample by 8')
w = waitforbuttonpress; 

% Set the desired output sample rate
outputSampleRate = round(sampleFrequency *8 );

% Upsample the audio signal
outputSignal = resample(audioSample, outputSampleRate, sampleFrequency);

% Write the downsampled signal to the output audio file
audiowrite('audio_upsampleby008.wav', outputSignal, outputSampleRate);

% Play audio vector over speakers 
if ((outputSampleRate > 1000)&&(outputSampleRate <256000))
    sound(outputSignal,outputSampleRate,16)
    pause(audioDuration) % pause duration of sample to listen
end

% Perform FFT
x = fft(outputSignal);   % DFT of signal

% Graph Filtered Frequency Domain 
subplot(4,1,3)
plot(real(x));
xlim([20 20000]);     % limit to human hearing range
ylim([0 1500])        % arbitrary amplitude limit
title('Freq-Domain Upsample by 8')
xlabel('freq (Hz)');
ylabel('Amplitude');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Up Sample by 16
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Wait for the user to press a button
disp('Press a button to Upsample by 16')
w = waitforbuttonpress; 

% Set the desired output sample rate
outputSampleRate = round(sampleFrequency *16 );

% Upsample the audio signal
outputSignal = resample(audioSample, outputSampleRate, sampleFrequency);

% Write the downsampled signal to the output audio file
audiowrite('audio_upsampleby016.wav', outputSignal, outputSampleRate);

% Play audio vector over speakers 
if ((outputSampleRate > 1000)&&(outputSampleRate <256000))
    sound(outputSignal,outputSampleRate,16)
    pause(audioDuration) % pause duration of sample to listen
end

% Perform FFT
x = fft(outputSignal);   % DFT of signal

% Graph Filtered Frequency Domain 
subplot(4,1,4)
plot(real(x));
xlim([20 20000]);     % limit to human hearing range
ylim([0 1500])        % arbitrary amplitude limit
title('Freq-Domain Upsample by 16')
xlabel('freq (Hz)');
ylabel('Amplitude');

