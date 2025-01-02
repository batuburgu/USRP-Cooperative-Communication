rx = comm.SDRuReceiver(...
               Platform ="B200", ...
              SerialNum ="31FD9BD", ...
              CenterFrequency =400*1e6,OutputDataType ="double",Gain=45);
                
rx.EnableBurstMode = true;
rx.NumFramesInBurst = 20;
rx.SamplesPerFrame = 3200; % 2 * length(IF_signal)

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

[signal, eps]=THAL_meyr_oeder_symbol_sync(8,64,transpose(filtered_sig));
signal = transpose(signal);
[fine_sync_signal, freq_offset] = fine_freq_phase_syncronization(signal, 32e6, 512);
%%
% Frame Header 
N_zc = 63; % Length of Zadoff Chu
cf = mod(N_zc,2);
q = 0; % Cyclically Shifting coeff
u = 1; % Root of Zadoff Chu Function
n = 0:N_zc - 1;
frame_header = exp(-1i*pi*u.*n.*(n + cf + 2*q) / N_zc); % Zadoff Chu Sequence as Frame Header
cross_corr = xcorr(frame_header, fine_sync_signal(end:-1:1));
plot(abs(cross_corr))

fh_indices = find(abs(cross_corr) > 30);

% Extracting Meaningful Data from Package
meaningful_data = fine_sync_signal(fh_indices(1) + 11: fh_indices(1) + 210);

% Decision Making
cond1 = real(meaningful_data) > 0 & imag(meaningful_data) > 0; % 1;1;
cond2 = real(meaningful_data) < 0 & imag(meaningful_data) > 0; % 1;0
cond3 = real(meaningful_data) < 0 & imag(meaningful_data) < 0; % 0;1
cond4 = real(meaningful_data) > 0 & imag(meaningful_data) < 0; % 1;1

decision = zeros(2,200);

decision(:, cond1) = repmat([0; 0], 1, sum(cond1(:)));
decision(:, cond2) = repmat([1; 0], 1, sum(cond2(:)));
decision(:, cond3) = repmat([0; 1], 1, sum(cond3(:)));
decision(:, cond4) = repmat([1; 1], 1, sum(cond4(:)));

bit1_err = bitstream(1,:) ~= decision(1,:);
bit2_err = bitstream(2,:) ~= decision(2,:);

ser = sum(bit1_err | bit2_err) / length(bitstream);





