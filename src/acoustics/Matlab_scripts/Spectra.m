%**********************************************************************
% Computation of Power Spectrum with Hanning window
% ARGUMENTS:
% x: time history of quantity of interest (e.g., P) 
% nfft: number of fft points
% fs: sampling frequency
% 
% OUTPUT: 
% fftx: array of complex fft of x (same size than f)
% mx: array of magnitude of the fft of x (same size than f)
% f: array of frequencies
%
%
% G. A. Brès - Cascade Technologies
%**********************************************************************

function [fftx,mx,f] = PowerSpectrum(x,nfft,fs)
%removing average
npts=length(x);
avg=0.0;
for j=1:npts
    avg=avg+x(j);
end
avg=avg/npts;
x=x-avg;

%############## Fast Fourrier Transform ##########################
%hanning window
if size(x,1)==1
    h=transpose(hann(npts,'periodic'));
else
    h=hann(npts,'periodic');
end   


%take the FFFT
fftx = fft(x.*h,nfft);


%remove the duplicate data and divide by input length
NumUniquePts = ceil((nfft+1)/2);
fftx = fftx(1:NumUniquePts)/length(x);

%Power
mx =abs(fftx);
mx = mx.^2;

%account for duplicate frequencies
if rem(nfft,2)
    mx(2:end) = mx(2:end)*2;
else
    mx(2:end-1) = mx(2:end-1)*2;
end

%array of frequencies
f=(0:NumUniquePts-1)*fs/nfft;

%correct amplitude and power for use of hanning window
% see http://www.mathworks.com/matlabcentral/newsreader/view_thread/237833
% amplitude = sum(h)/size(mic,1) = 0.5
% power = sum(h.^2)/size(mic,1) = 3/8
fftx = fftx/(0.5);
mx = mx/(3/8);
     

%##########################################################################
%% FROM http://www.mathworks.com/support/tech-notes/1700/1702.html
%##########################################################################
% Question:
% 
% How can you correctly scale the output of the FFT function to obtain a meaningful power versus frequency plot?
% 
% Answer:
% 
% Assume x is a vector containing your data. A sample vector used with this technical note is a 200 Hz sinusoid signal.
% % Sampling frequency
% Fs = 1024;
% % Time vector of 1 second
% t = 0:1/Fs:1;
% % Create a sine wave of 200 Hz.
% x =sin(2*pi*t*200);
% 
% First, you need to call the FFT function. For the fastest possible ffts, you will want to pad your data with enough zeros to make its length a power of 2. The built-in FFT function does this for you automatically, if you give a second argument specifying the overall length of the fft, as demonstrated below:
% % Use next highest power of 2 greater than or equal to length(x) to calculate fft.
% nfft = 2^(nextpow2(length(x)));
% % Take fft, padding with zeros so that length(fftx) is equal to nfft
% fftx = fft(x,nfft);
% 
% If nfft is even (which it will be, if you use the above two commands above), then the magnitude of the fft will be symmetric, such that the first (1+nfft/2) points are unique, and the rest are symmetrically redundant. The DC component of x is fftx(1) , and fftx(1+nfft/2)> is the Nyquist frequency component of x. If nfft is odd, however, the Nyquist frequency component is not evaluated, and the number of unique points is (nfft+1)/2 . This can be generalized for both cases to ceil((nfft+1)/2) .
% % Calculate the number of unique points
% NumUniquePts = ceil((nfft+1)/2);
% 
% % FFT is symmetric, throw away second half
% fftx = fftx(1:NumUniquePts);
% 
% Next, calculate the magnitude of the fft:
% % Take the magnitude of fft of x mx = abs(fftx);
% 
% Consider the fact that MATLAB does not scale the output of fft by the length of the input:
% % Scale the fft so that it is not a function of the length of x mx = mx/length(x);
% % Now, take the square of the magnitude of fft of x which has been scaled properly. % Take the square of the magnitude of fft of x. mx = mx.^2;
% % Since we dropped half the FFT, we multiply mx by 2 to keep the same energy. % The DC component and Nyquist component, if it exists, are unique and should not % be mulitplied by 2. if rem(nfft, 2) % odd nfft excludes Nyquist point   mx(2:end) = mx(2:end)*2; else   mx(2:end -1) = mx(2:end -1)*2; end
% 
% Now, create a frequency vector:
% % This is an evenly spaced frequency vector with NumUniquePts points. f = (0:NumUniquePts-1)*Fs/nfft;
% 
% Finally, generate the plot with a title and axis labels.
% % Generate the plot, title and labels. plot(f,mx); title('Power Spectrum of a 200Hz Sine Wave'); xlabel('Frequency (Hz)'); ylabel('Power');
% 
% Bringing this all together, you get the following M-file:
% % Sampling frequency Fs = 1024;
% % Time vector of 1 second t = 0:1/Fs:1;
% % Create a sine wave of 200 Hz.
% x = sin(2*pi*t*200);
% % Use next highest power of 2 greater than or equal to length(x) to calculate FFT.
% nfft= 2^(nextpow2(length(x)));
% % Take fft, padding with zeros so that length(fftx) is equal to nfft
% fftx = fft(x,nfft);
% % Calculate the numberof unique points NumUniquePts = ceil((nfft+1)/2);
% % FFT is symmetric, throw away second half
% fftx = fftx(1:NumUniquePts);
% % Take the magnitude of fft of x and scale the fft so that it is not a function of % the length of x mx = abs(fftx)/length(x);
% % Take the square of the magnitude of fft of x.
% mx = mx.^2;
% 
% % Since we dropped half the FFT, we multiply mx by 2 to keep the same energy.
% % The DC component and Nyquist component, if it exists, are unique and should not
% % be mulitplied by 2.
% if rem(nfft, 2) % odd nfft excludes Nyquist point
%   mx(2:end) = mx(2:end)*2;
% else
%   mx(2:end -1) = mx(2:end -1)*2;
% end
% % This is an evenly spaced frequency vector with NumUniquePts points.
% f = (0:NumUniquePts-1)*Fs/nfft;
% % Generate the plot, title and labels.
% plot(f,mx);
% title('Power Spectrum of a 200Hz Sine Wave');
% xlabel('Frequency (Hz)');
% ylabel('Power');
% 
% The resulting plot looks like the following:

