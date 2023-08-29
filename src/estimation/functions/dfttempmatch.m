

function [output CC2] = dfttempmatch(Iorg,Torg,usfac)
% function [output Greg] = dftregistration(buf1ft,buf2ft,usfac);
% Efficient subpixel image registration by crosscorrelation. This code
% gives the same precision as the FFT upsampled cross correlation in a
% small fraction of the computation time and with reduced memory
% requirements. It obtains an initial estimate of the crosscorrelation peak
% by an FFT and then refines the shift estimation by upsampling the DFT
% only in a small neighborhood of that estimate by means of a
% matrix-multiply DFT. With this procedure all the image points are used to
% compute the upsampled crosscorrelation.
% Manuel Guizar - Dec 13, 2007

% Portions of this code were taken from code written by Ann M. Kowalczyk
% and James R. Fienup.
% J.R. Fienup and A.M. Kowalczyk, "Phase retrieval for a complex-valued
% object by using a low-resolution image," J. Opt. Soc. Am. A 7, 450-458
% (1990).

% Citation for this algorithm:
% Manuel Guizar-Sicairos, Samuel T. Thurman, and James R. Fienup,
% "Efficient subpixel image registration algorithms," Opt. Lett. 33,
% 156-158 (2008).

% Inputs
% buf1ft    Fourier transform of reference image,
%           DC in (1,1)   [DO NOT FFTSHIFT]
% buf2ft    Fourier transform of image to register,
%           DC in (1,1) [DO NOT FFTSHIFT]
% usfac     Upsampling factor (integer). Images will be registered to
%           within 1/usfac of a pixel. For example usfac = 20 means the
%           images will be registered within 1/20 of a pixel. (default = 1)

% Outputs
% output =  [error,diffphase,net_row_shift,net_col_shift]
% error     Translation invariant normalized RMS error between f and g
% diffphase     Global phase difference between the two images (should be
%               zero if images are non-negative).
% net_row_shift net_col_shift   Pixel shifts between images
% Greg      (Optional) Fourier transform of registered version of buf2ft,
%           the global phase difference is compensated for.


% Partial-pixel shift

% First upsample by a factor of 2 to obtain initial estimate
% Embed Fourier data in a 2x larger array

if usfac == 1
    %% Zeropadding
    [I, T, fsize] = zeropadding(Iorg,Torg);
    sumT = sum(sum(abs(T)));
    
    %% Tansformation in den Freuqenzbereich
    Ifft = fft2(I);
    Tfft = fft2(T);
    
    %% Kreuzkorrelation über Multiplikation im Freuquenzbereich
    CCfft = fftshift(Ifft).*fftshift(Tfft);
    CC = ifft2(ifftshift(CCfft))./sumT;
    
    %% Kontrolle
%     figure
%     subplot(1,3,1)
%     imshow(I,[]);
%     subplot(1,3,2)
%     imshow(T,[]);
%     subplot(1,3,3)
%     imshow(CC,[]);
    
    %% Maximum finden
    [max1,loc1] = max(real(CC));
    [max2,loc2] = max(max1);
    rloc=loc1(loc2);cloc=loc2;
    
    %% Maximale Korrelation
    CCmax=CC(rloc,cloc);
    
    %% Phasenverschiebung:
    diffphase=atan2(imag(CCmax),real(CCmax));
    
    %% Zeilen- und Spaltenverschiebung
    [m,n] = size(CC);
    md2 = fix(m/2);
    nd2 = fix(n/2);
    if rloc > md2
        row_shift = rloc - 1 - m;
    else
        row_shift = rloc - 1;
    end
    if cloc > nd2
        col_shift = cloc - 1 - n;
    else
        col_shift = cloc - 1;
    end
    row_shift = row_shift + 1;
    col_shift = col_shift + 1;
    
    %% Output
    output=[real(CCmax),diffphase,row_shift,col_shift];
    CC2 = CC;
elseif usfac > 1
    %% Zeropadding
    [I, T, fsize, Ta] = zeropadding(Iorg,Torg);
    sumT = sum(sum(abs(T)));
        
    %% Tansformation in den Freuqenzbereich
    Ifft = fft2(I);
    Tfft = fft2(T);
    
    %% Kenngrößen fürs Upsampling
    [m,n]=size(Ifft);
    mlarge=m*2;
    nlarge=n*2;
    CCfft=zeros(mlarge,nlarge);
    
    %% Kreuzkorrelation über Multiplikation im Frequenzbereich
    CCfft(m+1-fix(m/2):m+1+fix((m-1)/2),n+1-fix(n/2):n+1+fix((n-1)/2)) = ...
        fftshift(Ifft).*fftshift(Tfft);
    
    % Compute crosscorrelation and locate the peak
    CC = ifft2(ifftshift(CCfft))./(sumT/2/2);
    
    %% Kontrolle
%     figure
%     subplot(1,3,1)
%     imshow(I,[]);
%     subplot(1,3,2)
%     imshow(T,[]);
%     subplot(1,3,3)
%     imshow(CC,[]);
    
    %% Maximum finden
    [max1,loc1] = max(real(CC));
    [max2,loc2] = max(max1);
    rloc=loc1(loc2);cloc=loc2;
    
    %% Maximale Korrelation
    CCmax=CC(rloc,cloc);
      
    %% Zeilen- und Spaltenverschiebung
    % Obtain shift in original pixel grid from the position of the
    % crosscorrelation peak
    [m,n] = size(CC); md2 = fix(m/2); nd2 = fix(n/2);
    if rloc > md2
        row_shift = rloc - m - 1;
    else
        row_shift = rloc - 1;
    end
    if cloc > nd2
        col_shift = cloc - n - 1;
    else
        col_shift = cloc - 1;
    end
    row_shift=row_shift/2 + 1;
    col_shift=col_shift/2 + 1;
    
    CC2 = CC;
    %% Weiteres Upsampling
    % If upsampling > 2, then refine estimate with matrix multiply DFT
    Tffta = fft2(Ta); 
    
    if usfac > 2
        %% Funktioniert
        %%% DFT computation %%%
        % Initial shift estimate in upsampled grid
        row_shift = round(row_shift*usfac)/usfac;
        col_shift = round(col_shift*usfac)/usfac;
        dftshift = fix(ceil(usfac*3)/2); %% Center of output array at dftshift+1
        
        % Matrix multiply DFT around the current shift estimate
        CC = conj(dftups(Tffta.*conj(Ifft),ceil(usfac*3),ceil(usfac*3),usfac,...
            dftshift-row_shift*usfac,dftshift-col_shift*usfac))/(md2*nd2*sumT);
        
        % Locate maximum and map back to original pixel grid
        [max1,loc1] = max(real(CC));
        [max2,loc2] = max(max1);
        rloc = loc1(loc2); cloc = loc2;
        CCmax = CC(rloc,cloc);
        
        rloc = rloc - dftshift - 1;
        cloc = cloc - dftshift - 1;
        row_shift = row_shift + rloc/usfac;
        col_shift = col_shift + cloc/usfac;
    end
    
    %% Phasenverschiebung:
    diffphase=atan2(imag(CCmax),real(CCmax));
    
    %% output array
    output=[real(CCmax),diffphase,row_shift,col_shift];
end

% % Compute registered version of buf2ft
% if (nargout > 1)&&(usfac > 0),
%     [nr,nc]=size(buf2ft);
%     Nr = ifftshift([-fix(nr/2):ceil(nr/2)-1]);
%     Nc = ifftshift([-fix(nc/2):ceil(nc/2)-1]);
%     [Nc,Nr] = meshgrid(Nc,Nr);
%     Greg = buf2ft.*exp(i*2*pi*(-row_shift*Nr/nr-col_shift*Nc/nc));
%     Greg = Greg*exp(i*diffphase);
% elseif (nargout > 1)&&(usfac == 0)
%     Greg = buf2ft*exp(i*diffphase);
% end
return

function out=dftups(in,nor,noc,usfac,roff,coff)
% function out=dftups(in,nor,noc,usfac,roff,coff);
% Upsampled DFT by matrix multiplies, can compute an upsampled DFT in just
% a small region.
% usfac         Upsampling factor (default usfac = 1)
% [nor,noc]     Number of pixels in the output upsampled DFT, in
%               units of upsampled pixels (default = size(in))
% roff, coff    Row and column offsets, allow to shift the output array to
%               a region of interest on the DFT (default = 0)
% Recieves DC in upper left corner, image center must be in (1,1)
% Manuel Guizar - Dec 13, 2007
% Modified from dftus, by J.R. Fienup 7/31/06

% This code is intended to provide the same result as if the following
% operations were performed
%   - Embed the array "in" in an array that is usfac times larger in each
%     dimension. ifftshift to bring the center of the image to (1,1).
%   - Take the FFT of the larger array
%   - Extract an [nor, noc] region of the result. Starting with the
%     [roff+1 coff+1] element.

% It achieves this result by computing the DFT in the output array without
% the need to zeropad. Much faster and memory efficient than the
% zero-padded FFT approach if [nor noc] are much smaller than [nr*usfac nc*usfac]

[nr,nc]=size(in);
% Set defaults
if exist('roff')~=1, roff=0; end
if exist('coff')~=1, coff=0; end
if exist('usfac')~=1, usfac=1; end
if exist('noc')~=1, noc=nc; end
if exist('nor')~=1, nor=nr; end
% Compute kernels and obtain DFT by matrix products
kernc=exp((-i*2*pi/(nc*usfac))*( ifftshift([0:nc-1]).' - floor(nc/2) )*( [0:noc-1] - coff ));
kernr=exp((-i*2*pi/(nr*usfac))*( [0:nor-1].' - roff )*( ifftshift([0:nr-1]) - floor(nr/2)  ));
out=kernr*in*kernc;
return

function [x, m, fsize, ma] = zeropadding(x, m)
xsize = size(x);
msize = size(m);

fsize = xsize;
% Neue Dimensionen für Zeropadding
mmid = floor(msize/2);
xsize = xsize + mmid;   % border to avoid edge effects

% Zeropadding
x = exindex(x, 1:xsize(1), 1:xsize(2), {0});
m = exindex(m, 1:xsize(1), 1:xsize(2), {0});
ma = m;

% Zyklische Anorndung des Templates ??
% recentre m so that y(1,1) corresponds to mask centred on x(1,1)
mc = 1 + 2*mmid;
me = mc + xsize - 1;
m = exindex(m, mc(1):me(1), mc(2):me(2), 'circular');
return
