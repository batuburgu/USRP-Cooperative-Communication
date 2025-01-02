%inputs:
%data: datas that will be fed into the costas loop, 1xN row vector, didn't
%try as column vector, it works great with row vector 

%masterclock: SDR masterclock configuration
%decimation: SDR decimation configuration


function [fine_sync_signal, freq_ofset]= fine_freq_phase_syncronization(data, masterclock,decimation)

fs=masterclock/decimation;%not really needed here, just needed here to monitor 
% the frequency offset errors, if needed

data=data./max(abs(data)); %normalization. Needed in order for the constellation to work

N=length(data);

%initial values
phase=0;
frequency=0;

%costas feedback loop parameters. they change the responsiveness of the
%loop. Change, play with them if you are brave enough to try and not on a deadline!
alpha=0.132;
beta=0.00932;

%output vectors, initilialized b4 the loop
fine_sync_signal=zeros(1,N);
freq_ofset=zeros(1,N);

%main loop
for i=1:1:N

fine_sync_signal(i)=data(i).*exp(-1j*phase);
error=error_calculation(fine_sync_signal(i));%its in our github folder

frequency=frequency+(beta*error);
phase=phase+frequency+(alpha*error);

freq_ofset(i)=frequency.*fs/(2*pi);

%phase normalization, not at all needed, you can even commend it, it
%doesn't matter
while(phase>=2*pi)
phase=phase-2*pi;
end

while(phase<0)
    phase=phase+2*pi;
end

end

end
