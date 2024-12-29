rx = comm.SDRuReceiver(Platform="B200",SerialNum="31FD9A5", ...
    CenterFrequency=400*1e6,Gain=12);

rxLog = dsp.SignalSink;
    for counter = 1:5000
data=step(rx);
rxLog(data);
    end
release(rx)


spectrumScope = spectrumAnalyzer('SampleRate',32000000/512); 
a=rxLog.Buffer;

recieve_filter=comm.RaisedCosineReceiveFilter("InputSamplesPerSymbol",8,"RolloffFactor",0.5,"Shape","Square root","FilterSpanInSymbols",10);

gr_T=recieve_filter.coeffs;

a=conv(a,1,"same");
c=bandpass(transpose(a),[15000 25000],32000000/512);
h=cos(2*pi*2000.*linspace(0,1,length(c)));
e=h.*c;
recieve_data=conv(e,gr_T.Numerator,"same");

scatterplot(recieve_data)

freq_vector= linspace(0,6250,length(transpose(anten_data)));
spectrum=fft(transpose(anten_data));
plot(freq_vector,abs(spectrum));

