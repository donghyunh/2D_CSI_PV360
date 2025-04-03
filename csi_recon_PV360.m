% load Bruker PV360 CSI data
% 
% Note:
%
% Bruker fid_proc.64 is not a raw k-space data set. 
% The spectral domain is in kt (i.e. FID). 
% But the spatial domains are in x and y (Fourier transform on the console). 
%
% I want to zero fill all domains. 
% So, after I readin fid_proc.64, I first fft the x and y domains to generate a true raw k-space matrix in (kt, kx, ky, nt), named "raw".
% Then I do line broadening of the kt domain and zero fill all domains.
% Finally, I fft all domains to generate a (f, x, y, nt) matrix, named "rc123".
%
%
% Salva's 2024 rat1_1 E18 CSI
% dimy "phase" (=12) and dimx "read" (=12) on the scannr (both phase encodings)
% FOV = (60, 60)

close all
clear all
clc

%cdata = read_2dseq;
cdata = read_2dseq('CSI/pdata/1');
[dimf dimy dimx dummy1 dummy2 dummy3 nt] = size(cdata)

% load Bruker "fid" data
fileID = fopen('CSI/pdata/1/fid_proc.64');
fid = fread(fileID, inf, 'double');
fclose(fileID);

% Make complex fid array cfid and reshape it to rcfid
% rcfid is (kt, x, y) NOT (kt, kx, ky)
cfid = fid(1:2:end)+i*fid(2:2:end);
rcfidxy = reshape(cfid, [dimf, dimx, dimy, nt]);  
rcfid = rcfidxy;

% Make raw data (kt, kx, ky, nt)
k2rcfid = fftshift(fft(fftshift(rcfid, 2), [], 2), 2);
raw = fftshift(fft(fftshift(k2rcfid, 3), [], 3), 3);

% Plot a raw FID to check blank points
figure
imagesc(squeeze(abs(raw(100, :, :, 1)))) % To identify the voxel(s) at the center of the k-space
plot(abs(raw(1:200, 7, 7, 1))); % Make a plot to identify the blank points
plot(sum(sum(abs(raw(:, 6:8, 7, 1)),2),3));

% shift blank
blank = 78; % receiver blanking points in the beginning of FID
rawb = zeros(size(raw)); 
rawb(1:(dimf-blank), :) = raw((blank+1):dimf, :);

% Define other recon parameters (per method)
cf = 100.66; % working frequency in MHz ($PVM_FrqWork)
cppm = 172; % center frequency in ppm ($PVM_FrqWorkPpm)
ppm = 40.6; % spectral BW in ppm ($PVM_SpecSW)
udimf = 2 * dimf;
udimx = 3 * dimx;
udimy = 3 * dimy;
lb = 50; % line broadening in Hz

% line broadening
tp = (1/(ppm*cf))*(1:dimf);
lbf = exp(-lb*tp);
for jf = 1:nt
    for jy = 1:dimy
        tmp = rawb(:,:,jy,jf);
        tmplb = (lbf.*tmp')';
        rawblb(:,:,jy,jf) = tmplb;
    end
end

% Visualize the line broadening
figure
plot(sum(sum(abs(rawb(:, 4:5, 4:5, 1)),2),3), 'b');  
hold
plot(sum(sum(abs(rawblb(:, 4:5, 4:5, 1)),2),3), 'r');
hold off

% Reconstruct all dimensions with zero fills
rc1=fftshift((ifft(rawblb, udimf,1)), 1);
rc12=fftshift(ifft(rc1,udimx, 2), 2);
rc123=fftshift(ifft(rc12,udimy, 3), 3);

% plot spectrum in NMR convention
figure
xppm = (ppm/udimf) * (1:udimf);
xppmf = fliplr(xppm) - (ppm/2-cppm);
plot(xppmf, mean(mean(mean(abs(rc123(:,:,:,:)), 2),3), 4)); set(gca,'xdir','reverse') 
xlow = cppm-ppm/2
xhi = cppm+ppm/2
xlim([floor(cppm-ppm/2) ceil(cppm+ppm/2)]);
xlabel('Spectral Bandwidth (ppm)');
ylabel('MR Signal (a.u.)');

%%___________________________%% image overlay %% __________________________

% read proton image 
% E20: highres anatomical images

%pimg = read_2dseq;
pimg = read_2dseq('T2/pdata/1');
[pdimx pdimy dummy1 dummy2 pslc] = size(pimg)
igray = (squeeze(pimg(:,:,dummy1,dummy2,12)))'; % Pick the middle slice to represent the CSI slice

