M = 4; % Modulation order
bit_per_symbol = log2(M); % Bits coded per symbol

N = 200; % Number of Sent Symbols 

% Frame Header Arguments
N_zc = 63; % Length of Zadoff Chu
cf = mod(N_zc,2);
q = 0; % Cyclically Shifting coeff
u = 1; % Root of Zadoff Chu Function
n = 0:N_zc-1;

bitstream = randi([0 1],bit_per_symbol,N); % Upper Column is the MSB
oversampling_rate = 8; % Pulse Shaping Oversampling Rate

power_of_twos = 0:1:bit_per_symbol-1; % Decimal Value of Bits

decimal_values = 2.^flip(power_of_twos); % Decimal Value of Bit Stack
index = decimal_values * bitstream; % Bitstream as Symbol Indexes

% Data Package
zero_bit_stream = zeros(1,10);
empty_bit_stream = repmat([1, -1], 1, 16);

frame_header_stream = exp(-1i*pi*u.*n.*(n + cf + 2*q) / N_zc); % Zadoff Chu Sequence as Frame Header

parity_bit_stream = repmat([1, -1], 1, 5);

% Modulated Data
signals = exp(1j*((2*pi*index/M)+pi/4)); % MPSK Signal Stream 

data = [zero_bit_stream, empty_bit_stream, frame_header_stream, parity_bit_stream, signals, parity_bit_stream, zero_bit_stream];

autocorr=xcorr(frame_header_stream, data(end:-1:1) );
plot(abs(autocorr))
% Pulse Shaping
upsampled_data = upsample(data,oversampling_rate);

txfilter = rcosdesign(0.55,10,8,"sqrt");

x=conv(upsampled_data,txfilter); % Transmitted Waveform
IF_frequency=1.3;
fs=3*IF_frequency;
Ts=1/fs;
time_vector=0:length(x)-1;

IF_signal=x.*exp(-1i*2*pi*IF_frequency*time_vector);

% spectrum_baseband=fft(x);
% spectrum_IF=fft(IF_signal);
% frequency_vector=linspace(0,length(x),length(spectrum_baseband));
% 
% figure(1)
% plot(abs(spectrum_baseband))
% 
% figure(2)
% plot(abs(spectrum_IF))

tx=comm.SDRuTransmitter("Platform","B200", ...
    "CenterFrequency",400*1e6,"SerialNum","31FD9A5","Gain",20);

% Transmit The Signal 1000 Times
for i=1:1:1000
    tx(transpose(IF_signal));
end

release(tx);