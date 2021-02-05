%% Huawei Challenge - Squinted SAR - Massimiliano Sica 10558133

clc, clear, close all
load('SquintedSAR_Data.mat');

%% First Step: Range Compression by means of Fourier Transform

c = 3e8;                                                                  %Speed of ligth
[m,n,ch] = size(D_SB);                                                    %get data matrix dimensions 
theta=deg2rad(35);                                                        %Radar squint angle
lookingAngle=deg2rad(0);                                                  %Angle at which we want to look, useful later 
lambda=c/f0;                                                              %Carrier wavelength 

%Find spacing between virtual antennas (considering only the first two)

virtual_1=(Tx_pos(1,1) + Rx_pos(1,1))/2;
virtual_2=(Tx_pos(1,1) + Rx_pos(2,1))/2;
dv=abs(virtual_2-virtual_1);                                               %needed for compensation

%FFT
conditioner = 1;                                                           % Change this value to see how the Fourier transform chnages with N
OSF = 2;                                                                   % Oversampling factor
N = conditioner * OSF * 2^(nextpow2(m));                                   % Number of FFT sample, OSF times the neares frequency of 2
fourier_ch=zeros(N,n,ch);

for channel=[1 2 3 4 9 10 11 12]
    fprintf("Doing %d / %d \n", channel,ch);
    fourier_ch (:,:,channel) = fft(D_SB(:,:,channel), N);                  % FFT along fast times with N points
    fourier_ch (:,:,channel)= conj(fourier_ch(:,:,channel));
    fourier_ch (:,:,channel) = fourier_ch (:,:,channel).*exp(+1j*4*pi*sin(theta-lookingAngle)*dv*(channel-1)/lambda); %Compensation
end
fourier=sum(fourier_ch,3);

% Compute the axis
df = fs/N;                                                                 %Fourier interspacing
freqAxis = (0:N-1)*df;                                                     %Fourier frequency  axis
fastTimeAxis = freqAxis/K;                                                 %Fast time axis
rangeAxis = fastTimeAxis*c/2;                                              %Conversion to range axis (slant range)

fourier = fourier(1:N/2,:);                                                %Cut half of the spectrum that is symmetric and recompute axis
freqAxis = freqAxis(1:N/2);
fastTimeAxis = fastTimeAxis(1:N/2);
rangeAxis = rangeAxis(1:N/2);                                              %Sampling space not so large otherwise aliasing

imagesc(abs(fourier))

%% Second Step: generation of the backprojection grid in Cartesian coordinates


B=5e9;                                                                      %from 76 to 81 Ghz
rho_r = c/(2*B);                                                            %range resolution

rho_x=rho_r;                                                                %impose azimutal resolution equal to the range one
dx = 0.75*rho_x;                                                            % grid spacing along the motion of the car
dy = 0.75*rho_r;                                                            % and in the orthogonal direction

x_axis = 0:dx:max(Sx);
y_axis = 0:dy:12;

modifier = 1;
Ls = (lambda/(2*rho_x))*max(y_axis) * modifier;                             %Sintetic aperture
angleResolution = lambda/(2*Ls);                                            %Angular resolution

% Backprojection grid

[X,Y] = ndgrid(x_axis, y_axis);


%% Backprojection loop 

Z = zeros(size(X));                                                         %z value of the pixels
I = zeros(size(X));

step = 1;                                                                   %if 100 is equivalent to PRF around 10hz
for ii =1:step:n   %5000 to check steady subaprture.          
    fprintf("Doing %d / %d \n", ii,n);
    delta_x = X-Sx(ii);
    delta_y = Y-Sy(ii);
    delta_z = Z-Sz(ii);
    D = sqrt(delta_x.^2 + delta_y.^2 + delta_z.^2);                         %Distance of the pixels
    anglesOfView = asin(delta_x./D);                                        %Angle of the pixels wrt car
    ind = abs(anglesOfView-lookingAngle) < (angleResolution);               %Indexes of the grid where the condition is satisfied
    %{
    temp=zeros(size(X));    %To visualize looking angle
    temp(ind)=1;
    figure()
    imagesc(temp);
    pause
    %}
    interpolatedData = interp1(rangeAxis, fourier(:,ii), D(ind));             %Data interpolation
    rephasedData = interpolatedData.*exp(+1j.*(4*pi/lambda)*D(ind));          %Rephase
    videophaseData = rephasedData.*exp(-1j.*((4*pi*K*(D(ind).^2))/(c^2)));    %Video phase compensation
    I(ind) = I(ind) + videophaseData;                                         %Summing (average)
end

absImage = abs(I);
Im = quantile(absImage, 0.99, 'all');                                         %To cut very high value, to scale the image
absImage(absImage > Im) = Im;

figure; imagesc(y_axis, x_axis, absImage); colorbar; axis xy equal tight
colormap('gray'), xlabel("Orthogonal to motion"), ylabel('Motion direction');
