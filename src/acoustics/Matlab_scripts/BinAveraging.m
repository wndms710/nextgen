%**********************************************************************
% Powerspectrum converter (SPL or PSD -> PSD with Constant Bandwidth averaging)
% ARGUMENTS:
% input: PSD or SPL 
% f: frequency
% fstart: first frequency for bin averaging
% fend: last frequency for bin averaging 
% bandwidth: constant bandwidth for bin averaging
% option = 1     -> use if input is PSD (i.e., energy/frequency like Pascal^2/Hz)
% option = other -> use if input is autospectrum (i.e. energy like Pascal^2)
% dB = 1         -> use if input/output in dB
% 
% OUTPUT: 
% fbin: bin frequency, namely  0.5*(fc_u - fc_l)
% fc: center frequency, namely sqrt(fc_u * fc_l)
% fc_l: lower bounds of the frequency bands
% fc_u: upper bounds of the frequency bands
% spectraBin: bin-averaged PSD in energy/frequency
%
% G. A. Brès - Cascade Technologies
%**********************************************************************

function [df,fbin, fc, fc_l, fc_u, PSDbin] = BinAveraging(input, f, fstart, fend, bandwidth, option, dB)
p2ref=4e-10;
%==================================
%convert dB to energy if necessary
if (dB == 1)
    input=p2ref*10.^(input/10);
end
% multiply by df if input is PSD
if (option == 1)
    input = input.*(f(2)-f(1));
end

% number of bins
nbin = floor ((fend-fstart)/(bandwidth));
nstart = ceil(fstart/bandwidth);
if (nstart==0) 
    nstart=1;
end
;for j=1:nbin
    fbin(j)=(nstart+j-1)*bandwidth;
    fc_l(j)= fbin(j) - 0.5*bandwidth;
    fc_u(j)= fbin(j) + 0.5*bandwidth;
    fc(j)=sqrt(fc_u(j)*fc_l(j));
    % find effective frequency range 
    lfreq = find( (f>=fc_l(j)) & (f<fc_u(j)) );
    % compute effective bin width
    df(j)=max(f(lfreq))-min(f(lfreq))+f(2)-f(1);
    % compute bin-averaged PSD
    PSDbin(j)= sum(input(lfreq))/df(j);
end

% convert back to dB if necessary
if (dB==1)
    PSDbin=10*log10(PSDbin/p2ref);
end
