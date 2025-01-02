M = 4; % Modulation order
bit_per_symbol = log2(M); % Bits coded per symbol

N = 200; % Number of Sent Symbols 
bitstream = randi([0 1],bit_per_symbol,N); %First row is the LSB 

oversampling_rate = 8; % Pulse Shaping Oversampling Rate

power_of_twos = 0:1:bit_per_symbol-1; % Decimal Value of Bits

decimal_values = 2.^flip(power_of_twos); % Decimal Value of Bit Stack
index = decimal_values * bitstream; % Bitstream as Symbol Indexes

signals = exp(1j*((2*pi*index/M)+pi/4)); % MPSK Signal Stream 

% Pulse Shaping
signals=upsample(signals,oversampling_rate);

txfilter = rcosdesign(0.55,10,8,"sqrt");

x=conv(signals,txfilter); % Transmitted Waveform

IF_frequency=1.3;
fs=3*IF_frequency;
Ts=1/fs;
% time_vector=0:Ts:((Ts)*(length(x)-1))+1;
time_vector=0:length(x)-1;

% IF_signal=x.*cos(2*pi*IF_frequency*time_vector);
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
    "CenterFrequency",400*1e6,"SerialNum","31FD9A5","Gain",45);
for i=1:1:1000
    tx(transpose(IF_signal));
end

release(tx);
