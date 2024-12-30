rx = comm.SDRuReceiver(...
               Platform ="B200", ...
              SerialNum ="31FD9BD", ...
              CenterFrequency =400*1e6,OutputDataType ="double",Gain=45);
              
rx.EnableBurstMode = true;
rx.NumFramesInBurst = 20;
rx.SamplesPerFrame = 5000;
for i=1:1:20
    rxLog = dsp.SignalSink;
    
    [data,overrun] = rx();
    rxLog(data)
end
    a=rxLog.Buffer;


IF_frequency=20*1e3;
fs=3*IF_frequency;
Ts=1/fs;
time_vector=0:1:length(a)-1;


baseband_signal=transpose(a).*exp(1j*2*pi*1.3*time_vector);

rxfilter = rcosdesign(0.55,10,8,"sqrt");

filtered_sig=conv(transpose(baseband_signal),rxfilter);
scatterplot(baseband_signal)
release(rx);

[signal, eps]=THAL_meyr_oeder_symbol_sync(8,64,transpose(filtered_sig));



scatterplot(signal)