%% Receiving from Source
rx = comm.SDRuReceiver(...
               Platform ="B200", ...
              SerialNum ="31FD9A5", ...
              CenterFrequency =400*1e6,OutputDataType ="double",Gain=45);
                
rx.EnableBurstMode = true;
rx.NumFramesInBurst = 20;
rx.SamplesPerFrame = 10512; % 2 * length(IF_signal)

rxLog = dsp.SignalSink;
for i=1:5
    [data,overrun] = rx();
    rxLog(data)
end
a=rxLog.Buffer;

figure(1);
plot(real(a));
hold on
plot(imag(a));

IF_frequency=20*1e3;
fs=3*IF_frequency;
Ts=1/fs;
time_vector=0:1:length(a)-1;

baseband_signal=transpose(a).*exp(1j*2*pi*1.3*time_vector);

rxfilter = rcosdesign(0.55,10,8,"sqrt");

filtered_sig=conv(transpose(baseband_signal),rxfilter);

% Freq Sync
freq_sync_signal = coarse_carrier_sync(transpose(filtered_sig),32e6,512,4);

% Time Sync
[signal, eps]=THAL_meyr_oeder_symbol_sync(8,64,freq_sync_signal);
signal = transpose(signal);

% Frame Header 
N_zc = 63; % Length of Zadoff Chu
cf = mod(N_zc,2);
q = 0; % Cyclically Shifting coeff
u = 1; % Root of Zadoff Chu Function
n = 0:N_zc - 1;
frame_header = exp(-1i*pi*u.*n.*(n + cf + 2*q) / N_zc); % Zadoff Chu Sequence as Frame Header
cross_corr = xcorr(frame_header, signal(end:-1:1));
plot(abs(cross_corr))

fh_indices = find(abs(cross_corr) > 10 * mean(abs(cross_corr)));

% Frame and Phase Sync
ser = ones(1,length(fh_indices) - 1);

source_information_data = zeros(length(fh_indices)-1,512);
for i=1:length(fh_indices)-1

    % Extracting Meaningful Data from Package
    information_data = signal(fh_indices(i) + 11: fh_indices(i) + 522);

    % Frame Sync
    parity_phase_shift_sum = 0;
    first_parities = signal(fh_indices(i) + 1: fh_indices(i) + 10);
    
    % Channel Estimation
    parity_amp_avg = sum(abs(first_parities)) / length(first_parities);
    for j=1:2:10
        %parity_phase_shift_sum = parity_phase_shift_sum + atan(imag(first_parities(j))/real(first_parities(j)));
        parity_phase_shift_sum = parity_phase_shift_sum + atan2(imag(first_parities(j)),real(first_parities(j)));
        
    end
    parity_phase_shift_mean = (parity_phase_shift_sum/5);

    channel = parity_amp_avg * exp(-1i * parity_phase_shift_mean);

    % Phase Shift
    information_data = information_data .* exp(-1i * parity_phase_shift_mean);
    %information_data = information_data .* exp(-1i * pi); % for atan function
    information_data = transpose(information_data);
    source_information_data(i,:) = information_data;

    scatterplot(information_data)
end

%% Receiving from Relay
rxLog2 = dsp.SignalSink;
for i=1:5
    [data,overrun] = rx();
    rxLog2(data)
end

a=rxLog2.Buffer;

figure(1);
plot(real(a));
hold on
plot(imag(a));

IF_frequency=20*1e3;
fs=3*IF_frequency;
Ts=1/fs;
time_vector=0:1:length(a)-1;

baseband_signal=transpose(a).*exp(1j*2*pi*1.3*time_vector);

rxfilter = rcosdesign(0.55,10,8,"sqrt");

filtered_sig=conv(transpose(baseband_signal),rxfilter);

release(rx);
% Freq Sync
freq_sync_signal = coarse_carrier_sync(transpose(filtered_sig),32e6,512,4);
scatterplot(freq_sync_signal)
% Time Sync
[signal, eps]=THAL_meyr_oeder_symbol_sync(8,64,freq_sync_signal);
signal = transpose(signal);

% Frame Header 
N_zc = 63; % Length of Zadoff Chu
cf = mod(N_zc,2);
q = 0; % Cyclically Shifting coeff
u = 1; % Root of Zadoff Chu Function
n = 0:N_zc - 1;
frame_header = exp(-1i*pi*u.*n.*(n + cf + 2*q) / N_zc); % Zadoff Chu Sequence as Frame Header
cross_corr = xcorr(frame_header, signal(end:-1:1));
plot(abs(cross_corr))

fh_indices = find(abs(cross_corr) > 10 * mean(abs(cross_corr)));

relay_information_data = zeros(length(fh_indices)-1,512);
for i=1:length(fh_indices)-1

    % Extracting Meaningful Data from Package
    information_data = signal(fh_indices(i) + 11: fh_indices(i) + 522);

    % Frame Sync
    parity_phase_shift_sum = 0;
    first_parities = signal(fh_indices(i) + 1: fh_indices(i) + 10);
    
    % Channel Estimation
    parity_amp_avg = sum(abs(first_parities)) / length(first_parities);
    for j=1:2:10
        %parity_phase_shift_sum = parity_phase_shift_sum + atan(imag(first_parities(j))/real(first_parities(j)));
        parity_phase_shift_sum = parity_phase_shift_sum + atan2(imag(first_parities(j)),real(first_parities(j)));
        
    end
    parity_phase_shift_mean = (parity_phase_shift_sum/5);

    channel = parity_amp_avg * exp(-1i * parity_phase_shift_mean);

    % Phase Shift
    information_data = information_data .* exp(-1i * parity_phase_shift_mean);
    %information_data = information_data .* exp(-1i * pi); % for atan function
    information_data = transpose(information_data);
    relay_information_data(i,:) = information_data;

    scatterplot(information_data)
end
%% EGC and Decision Making
ser = ones(1,length(fh_indices) - 1);

total_information_data = source_information_data+relay_information_data ;

for k = 1:length(fh_indices)-1
    cond1 = real(total_information_data(k,:)) > 0 & imag(total_information_data(k,:)) > 0; % 0;0
    cond2 = real(total_information_data(k,:)) < 0 & imag(total_information_data(k,:)) > 0; % 0;1
    cond3 = real(total_information_data(k,:)) < 0 & imag(total_information_data(k,:)) < 0; % 1;0
    cond4 = real(total_information_data(k,:)) > 0 & imag(total_information_data(k,:)) < 0; % 1;1

    decision = zeros(2,200);
    
    decision(:, cond1) = repmat([0; 0], 1, sum(cond1(:)));
    decision(:, cond2) = repmat([0; 1], 1, sum(cond2(:)));
    decision(:, cond3) = repmat([1; 0], 1, sum(cond3(:)));
    decision(:, cond4) = repmat([1; 1], 1, sum(cond4(:)));
    
    % bit1_err = bitstream(1,:) ~= decision(1,:);
    % bit2_err = bitstream(2,:) ~= decision(2,:);
    % ser(k) = sum(bit1_err | bit2_err) / length(bitstream);

end


    
