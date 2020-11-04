%% MIMO RADAR design
clear all
close all
clc

%% creating the positions of the center of the array
center_x = 0;

center_y = 0;

%% computing the number of antennas (according to res) and dx
freq = 9.6e9;

c = 3e8;

lambda= c/freq;

rho_x = 2;

R=sqrt((35-center_x)^2 + (150-center_y)^2);

L_tot=(lambda*R)/(2*rho_x);                   

theta_max = atan(35/100);

dx=(lambda)/(4*sin(theta_max));           
N=ceil((L_tot/dx)+1);

N=55;

L_tot=(N-1)*dx;

%% Creating virtual array

positions = linspace(-floor(N/2)*dx, floor(N/2)*dx, N);

%% Creating Real array: 5 antennas tx, 11 receiving

Ntx=5;

Nrx=11;

pos_rx=zeros(1,Nrx);

%Note: this part of code works for odd numbers of transmitters, due to indexes reasons 

central_tx=0;
index_central = ceil(Ntx/2);

% From third transmitter, that since is the center of simmetry of the Tx
% array we place in zero, we compute the position of the Rx antennas
% accordingly to the virtual antennas from 23 to 33;

for i=1:length(pos_rx)
    pos_rx(i)=2*positions(((index_central-1)*Nrx)+i)-central_tx; %inverted relation of midpoint
end

dx_rx=abs(pos_rx(2)-pos_rx(1));

%get second transmitter, in the same but inverse way

next_tx=2*positions((index_central*Nrx)+1)-pos_rx(1);

%derive dtx and create the transmitter array

dx_tx=abs(next_tx-central_tx);

pos_tx=linspace(0,(Ntx-1)*dx_tx,Ntx);

pos_tx=pos_tx-(dx_tx*(index_central-1));


%% Check if Real array end up to virtual one

check=zeros(1,N);

for i = 1:length(pos_tx)
    for j = 1:length(pos_rx)
        check(((i-1)*length(pos_rx))+j)=(pos_tx(i)+pos_rx(j))/2;
    end
end

dx2=abs(check(1)-check(2));

figure(1)
scatter(pos_tx,zeros(1,length(pos_tx))), hold on
scatter(pos_rx,zeros(1,length(pos_rx)),'r'), hold on
scatter(check,ones(1,length(check)))
%assert (round(dx,5)==round(dx2,5), 'Dx of virtual array not equal to the one initially calculated!');
%assert (min(round(positions,5)==round(check,5)), 'Virtual array from real one not equal to the one initially calculated!');

%% Single target without noise: DOA estimaton

type = 1; %to select one option in the switch

switch type
    case 1
        %place target randomly in the FOV
        pos_x = -35 + (35+35).*rand(1,1);
        pos_y =  100 + (150-100).*rand(1,1);
    case 2
        %good case
        pos_x = -0.7165;
        pos_y = 122.2793;
    case 3
        %perfect broadside
        pos_x = 0;
        pos_y = 145;
    case 4
        %endfire
        pos_x = 35;
        pos_y = 100.1;
end

%tx and rx real antennas were placed in the previus section
%Calculare the real angle of arrival with geometry

real_theta = atan(pos_x/pos_y);
real_theta_deg=rad2deg(real_theta);

%Propagating the signal

R0=sqrt((pos_x-pos_tx).^2 + (pos_y-0)^2);
R1=sqrt((pos_x-pos_rx).^2 + (pos_y-0)^2);

rho_a = randn + 1i*randn;

signal=zeros(1,Ntx*Nrx);

for i=1:Ntx
    for j=1:Nrx
        signal(((i-1)*Nrx)+j)=(rho_a/(R0(i)+R1(j)))*exp(-1i*2*pi*(R0(i)+R1(j))/lambda);
    end
end

fs=1/dx;

N = 4096;

df = fs/N;

freqAxis = (-N/2:N/2-1)*df;

freqDomain = fftshift(fft(signal, N));

figure(2)

plot(freqAxis, abs(freqDomain).^2);

[M,I]= max(abs(freqDomain).^2);

max_freq = freqAxis(I);

est_theta=asin((lambda*max_freq)/2);

est_theta_deg=rad2deg(est_theta);

%% DOA estimation Multiple Targets.

while 1
    pos_x1 = -35 + (35+35).*rand(1,1);
    pos_y1 =  100 + (150-100).*rand(1,1);
    pos_x2 = -35 + (35+35).*rand(1,1);
    pos_y2 =  100 + (150-100).*rand(1,1);
    if sqrt(((pos_x1-pos_x2)^2)+((pos_y1-pos_y2)^2)) >rho_x %< to see what happens when they are under the resolution threshold
        break
    end
end

real_theta1 = atan(pos_x1/pos_y1);
real_theta1_deg=rad2deg(real_theta1);
real_theta2 = atan(pos_x2/pos_y2);
real_theta2_deg=rad2deg(real_theta2);

R01=sqrt((pos_x1-pos_tx).^2 + (pos_y1-0)^2);
R11=sqrt((pos_x1-pos_rx).^2 + (pos_y1-0)^2);
R02=sqrt((pos_x2-pos_tx).^2 + (pos_y2-0)^2);
R12=sqrt((pos_x2-pos_rx).^2 + (pos_y2-0)^2);

rho_a1 = randn + 1i*randn;
rho_a2 = randn + 1i*randn;

signal1=zeros(1,Ntx*Nrx);
signal2=zeros(1,Ntx*Nrx);

for i=1:Ntx
    for j=1:Nrx
        signal1(((i-1)*Nrx)+j)=(rho_a/(R01(i)+R11(j)))*exp(-1i*2*pi*(R01(i)+R11(j))/lambda);
        signal2(((i-1)*Nrx)+j)=(rho_a/(R02(i)+R12(j)))*exp(-1i*2*pi*(R02(i)+R12(j))/lambda);
    end
end

%We add the two signals since the FFT is linear

signal_tot=signal1+signal2;

fs=1/dx;

N = 4096;

df = fs/N;

freqAxis = (-N/2:N/2-1)*df;

freqDomain = fftshift(fft(signal_tot, N));

figure(3)

plot(freqAxis, abs(freqDomain).^2);

%we find the peacks

[M,I] = findpeaks(abs(freqDomain).^2); %to find highest peacks

[Ma,In]=maxk(M,2);  %to extract k=2 highest peacks

indexes=zeros(1,length(In));

for k=1:length(In)
    indexes(k)=I(In(k));
end

max_freq=freqAxis(indexes);

est_theta_all=asin((lambda*max_freq)/2);

est_theta_all_deg=rad2deg(est_theta_all);



%% DOA of Wideband Signal

load handel.mat

filename = 'handel.wav';
audiowrite(filename,y,Fs);
clear y Fs 

[y,Fs] = audioread('handel.wav'); %we get y matrix sampled at Fs frequency, since we have one channel, we have m*1 matrix for y

figure(4)
periodogram(y);
 
info= audioinfo('handel.wav');

[y,Fs] = audioread('handel.wav'); %we get y matrix sampled at Fs frequency, since we have one channel, we have m*1 matrix for y
sound(y,Fs);

prop_vel= 343;
attenuation =0+1.*rand(1,1);

pos_tar_x = -35 + (35+35).*rand(1,1);
pos_tar_y =  100 + (150-100).*rand(1,1);
target_pos = [pos_tar_x;pos_tar_y;0];
antenna1 = [-8.18;0;0];
antenna2 = [8.18;0;0];
real_angle=atand(pos_tar_x/pos_tar_y);

R1_wb = sqrt((target_pos(1)-antenna1(1)).^2 + (target_pos(2)-antenna1(2))^2);
delay1 = R1_wb/prop_vel;

R2_wb = sqrt((target_pos(1)-antenna2(1)).^2 + (target_pos(2)-antenna2(2))^2);
delay2 = R2_wb/prop_vel;


Y1 = delayseq(y, delay1, Fs)*attenuation;
Y2 = delayseq(y, delay2, Fs)*attenuation;

figure(5)
periodogram([y,Y1])
figure(6)
periodogram([y,Y2])

t_plot1_ds=0:seconds(1/Fs):seconds(info.Duration);
t_plot1_ds=t_plot1_ds(1:end-1);
figure(7)
plot(t_plot1_ds,transpose(real(y)))
hold on;
plot(t_plot1_ds,transpose(real(Y1)))
legend('original','propagated')

figure(8)
plot(t_plot1_ds,transpose(real(y)))
hold on;
plot(t_plot1_ds,transpose(real(Y2)))
legend('original','propagated')

figure(9)
plot(t_plot1_ds,transpose(real(Y1)))
hold on;
plot(t_plot1_ds,transpose(real(Y2)))
legend('Rx at mic 1','Rx at mic 2')

C = xcorr(Y1,Y2);
[~,index]=max(C);
delay_time=(index-ceil(length(C)/2))/Fs;

wb_est_theta_ds = asind((delay_time)*prop_vel/abs(antenna2(1)-antenna1(1)));


