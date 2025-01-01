%%% data input:1xN row vector, does not work with columm
%%% masterclock: masterclock configuration of SDR
%%% decimation: Decimation configuration of SDR

function [freq_sync_signal] = coarse_carrier_sync(data,masterclock,decimation)
x_power_4 = data.^4;

x_fft=fftshift(fft(x_power_4));
[~,index]=max(abs(x_fft));

fs=masterclock/decimation;  %masterclock/decimation_factor from sdr
Ts=1/fs;

%frequency and time vectors for the corresponding sampling frequency
freq_vector=linspace(-fs/2,fs/2,length(data)); 
time_vector=0:Ts:Ts*(length(data)-1);

f0=freq_vector(index)/4;%should be divided by the order of the modulation scheme, 4 for QPSK

freq_sync_signal=data.*exp(-1j*2*pi*f0*time_vector);
end












