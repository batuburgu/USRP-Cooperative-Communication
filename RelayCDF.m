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
    
    % Decoding
    first_parities = signal(fh_indices(i) + 7: fh_indices(i) + 16);   
    
    % Phase Sync
    parity_phase_shift_sum = 0;
    for j=1:2:10
        parity_phase_shift_sum = parity_phase_shift_sum + atan2(imag(first_parities(j)),real(first_parities(j)));
    end
    parity_phase_shift_mean = (parity_phase_shift_sum/5);

    information_data = signal(fh_indices(i) + 17: fh_indices(i) + 216);
    frame_number_data = signal(fh_indices(i) + 1: fh_indices(i) + 6);

    % Phase Shift
    information_data = information_data .* exp(-1i * parity_phase_shift_mean);
    frame_number_data = frame_number_data .* exp(-1i * parity_phase_shift_mean);
    %information_data = information_data .* exp(-1i * pi); % for atan function
    information_data = transpose(information_data);

scatterplot(information_data)

    % Decision Making
    cond1 = real(information_data) > 0 & imag(information_data) > 0; % 0;0
    cond2 = real(information_data) < 0 & imag(information_data) > 0; % 0;1
    cond3 = real(information_data) < 0 & imag(information_data) < 0; % 1;0
    cond4 = real(information_data) > 0 & imag(information_data) < 0; % 1;1
    
    decision = zeros(2,200);
    
    decision(:, cond1) = repmat([0; 0], 1, sum(cond1(:)));
    decision(:, cond2) = repmat([0; 1], 1, sum(cond2(:)));
    decision(:, cond3) = repmat([1; 0], 1, sum(cond3(:)));
    decision(:, cond4) = repmat([1; 1], 1, sum(cond4(:)));
    
    frame_number_stream = zeros(1,length(frame_number_data));
    frame_number_stream(real(frame_number_data) > 0.5) = 1;
    frame_number_stream(1) = 1;

    %retransmission
    M = 4; % Modulation order
    bit_per_symbol = log2(M);
    
    power_of_twos = 0:1:bit_per_symbol-1; % Decimal Value of Bits
    
    decimal_values = 2.^flip(power_of_twos); % Decimal Value of Bit Stack
    index = decimal_values * decision; % Bitstream as Symbol Indexes
    
    % Data Package
    zero_bit_stream = zeros(1,10);
    empty_bit_stream = repmat([1, -1], 1, 16);
    
    frame_header_stream = exp(-1i*pi*u.*n.*(n + cf + 2*q) / N_zc); % Zadoff Chu Sequence as Frame Header
    
    parity_bit_stream = repmat([1, -1], 1, 5);
    
    signals = exp(1j*((2*pi*index/M)+pi/4)); % MPSK Signal Stream     

    data = [zero_bit_stream, empty_bit_stream, frame_header_stream, frame_number_stream, parity_bit_stream, signals, parity_bit_stream, zero_bit_stream];

    % Forward
    oversampling_rate = 8; 
    upsampled_data = upsample(transpose(data),oversampling_rate);
    
    txfilter = rcosdesign(0.55,10,8,"sqrt");
    
    x=conv(upsampled_data,txfilter); 
    IF_frequency=1.3;
    fs=3*IF_frequency;
    Ts=1/fs;
    time_vector=0:length(x)-1;
  
    IF_signal=transpose(x).*exp(-1i*2*pi*IF_frequency*time_vector);
    
    for k=1:1:400
    tx(transpose(IF_signal));
    end

end
release(tx);