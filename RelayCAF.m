rx = comm.SDRuReceiver(...
               Platform ="B200", ...
              SerialNum ="31FD9BD", ...
              CenterFrequency =400*1e6,OutputDataType ="double",Gain=45);
                
rx.EnableBurstMode = true;
rx.NumFramesInBurst = 20;
rx.SamplesPerFrame = 5520; % 2 * length(IF_signal)

rxLog = dsp.SignalSink;
for i=1:4
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

release(rx);

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

% Frame Sync
cross_corr = xcorr(frame_header, signal(end:-1:1));
plot(abs(cross_corr))
fh_indices = find(abs(cross_corr) > (10*mean(abs(cross_corr))));

% CAF Process
tx=comm.SDRuTransmitter("Platform","B200", ...
    "CenterFrequency",400*1e6,"SerialNum","31FD9BD","Gain",45);

for i=2:length(fh_indices)-1  
    
    % Amplification
    first_parities = signal(fh_indices(i) + 1: fh_indices(i) + 10);
    last_parities = signal(fh_indices(i) + 211: fh_indices(i) + 220);
    parities_amp_avg = ((sum(abs(first_parities)) + sum(abs(last_parities))) / (length(first_parities) + length(last_parities)));
    frame = signal(fh_indices(i) - 104: fh_indices(i) + 230);
    amplified_frame = (frame ./ parities_amp_avg);
    
    % Phase Sync
    parity_phase_shift_sum = 0;
    for j=1:2:10
        parity_phase_shift_sum = parity_phase_shift_sum + atan2(imag(first_parities(j)),real(first_parities(j)));
    end
    parity_phase_shift_mean = (parity_phase_shift_sum/5);
    frame = frame .* exp(-1i * parity_phase_shift_mean);

    frame(fh_indices(i) + 1) = 1;

    % Forward
    oversampling_rate = 8; 
    upsampled_data = upsample(transpose(amplified_frame),oversampling_rate);
    
    txfilter = rcosdesign(0.55,10,8,"sqrt");
    
    x=conv(upsampled_data,txfilter); 
    IF_frequency=1.3;
    fs=3*IF_frequency;
    Ts=1/fs;
    time_vector=0:length(x)-1;
  
    IF_signal=x.*exp(-1i*2*pi*IF_frequency*time_vector);
    
    for k=1:1:400
    tx(transpose(IF_signal));
    end

end
release(tx);