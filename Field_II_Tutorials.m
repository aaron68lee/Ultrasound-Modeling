clc; clear all; close all;
path(path, 'C:\Users\yoshi\OneDrive\Documents\MATLAB\Field_II_ver_3_30_windows.tar\Field_II_ver_3_30_windows')
field_init(0); % field_init(0) suppresses output


%% General Program Setup
%{
    set geometric parameters
        element width does not include kerf
    set exitation wave pulses
    
    % for trasmission
    set_sampling
    xdc_linear_array
    xdc_impulse
    xdc_excitation
    
    % for receiving 
    xdc_linear_array
    xdc_impulse

    % use beamforming to reconstruct image from multiple scan lines
    xdc_center_focus
    xdc_apodization
    xdc_

    y = abs(hilbert(x))
    y = 20 log10(x)
    y = interpLateral(x, 20)
    sync time between scan lines
    separate data from code
    
    xdc_free for both transmitter and receiver
%}

example = 8;
switch(example)

    case 1
%% 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%% Point-spread Function Heat Graph %%%%%%%%%%%%%%%%%%%%%%%%%%%

%  Define a concave piston transducer,
%  set its impulse response and excitation
%  and calculate its point spread function
%
%  Note that the field_init routine must be
%  calculated before the routine is called
%
%  Version 1.0, June 29, 2001 by JAJ

%  Set initial parameters

R=8/1000;             %  Radius of transducer [m]
Rfocus=80/1000;       %  Geometric focus point [m]
ele_size = 1/100000;      %  Size of mathematical elements [m]
f0=9e6;               %  Transducer center frequency [Hz]
fs=100e6;             %  Sampling frequency [Hz]

%  Define the transducer

Th = xdc_concave (R, Rfocus, ele_size);

%  Set the impulse response and excitation of the emit aperture

impulse_response=sin(2*pi*f0*(0:1/fs:2/f0));
impulse_response=impulse_response.*hanning(max(size(impulse_response)))';
xdc_impulse (Th, impulse_response);

excitation=sin(2*pi*f0*(0:1/fs:2/f0));
xdc_excitation (Th, excitation);

%  Calculate the pulse echo field and display it

xpoints=(-10:0.2:10);
[RF_data, start_time] = calc_hhp (Th, Th, [xpoints; zeros(1,101); 30*ones(1,101)]'/1000);

%  Make a display of the envelope

figure(1)
env=abs(hilbert(RF_data(1:5:600,:)));
env=20*log10(env/max(max(env)));
[N,M]=size(env);
env=(env+60).*(env>-60) - 60;
mesh(xpoints, ((0:N-1)/fs + start_time)*1e6, env)
ylabel('Time [\mus]')
xlabel('Lateral distance [mm]')
title('Pulse-echo field from 8 mm concave transducer at 30 mm')
axis([-10 10 38.41 39.6 -60 0])
view([-14 80])

%  Create a smaller image for the thumbnail

figure(2)
imagesc(env)
colormap(hot)
axis off

%  Free the aperture after use

xdc_free (Th)

    case 2
%% 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Point-spread Fucntion Phantoms %%%%%%%%%%%%%%%%%%%%%%%%%%


folder  = 'C:\Users\yoshi\OneDrive\Documents\MATLAB\psf_example.tar\RFrqtc';
list    = dir(fullfile(folder, '*.m'));
nFile   = length(list);
success = false(1, nFile);

for k = 1:nFile
  file = list(k).name;
  try
    run(fullfile(folder, file));
    success(k) = true;
  catch
    fprintf('failed: %s\n', file);
  end
end

    case 3
%% 3 %%%%%%%%%%%%%%%%%%%%%% Calculation of Intensities and Peak Pressure %%%%%%%%%%%%%%%%%%%%%%%%%%


%  Simulate an array using Field II 
%  Calculate the intensity profile along the acoustical axis for
%  the transducer 
% 
%  Version 1.0, 27/6-1997, JAJ 
%  Version 1.1, 27/3-2000, JAJ: Transducer impulse response added
%  Version 1.1, 13/8-2007, JAJ: Change of display during calculation
 
% 
%  Array: 65 elements, width: lambda/2, Height: 5 mm, 
%  Type:  Elevation focused, linear array 

fetal=1;                 %  Whether to use cardiac or fetal intensities

%  Set values for the intensity 

if (fetal==1)

  %  For fetal 
  
  f0=5e6;                  %  Transducer center frequency [Hz] 
  M=2;                     %  Number of cycles in pulse 
  Ispta=0.170*100^2;       %  Fetal intensity: Ispta [w/m^2] 
  %Ispta=0.046*100^2;      %  Fetal intensity In Situ: Ispta [w/m^2] 
  Itype='Fetal';           %  Intensity type used 
else
 
 %  For cardiac 
 
  f0=3e6;                  %  Transducer center frequency [Hz] 
  M=8;                     %  Number of cycles in pulse 
  Ispta=0.730*100^2;       %  Cardiac intensity: Ispta [w/m^2] 
  Itype='Cardiac';         %  Intensity type used 
  end

%  Generate the transducer apertures for send and receive 
 
fs=250e6;                %  Sampling frequency [Hz] 
c=1540;                  %  Speed of sound [m/s] 
lambda=c/f0;             %  Wavelength 
width=lambda/2;          %  Width of element 
element_height=5/1000;   %  Height of element [m] 
kerf=lambda/10;          %  Kerf [m] 
focus=[0 0 60]/1000;     %  Fixed focal point [m] 
elefocus=1;              %  Whether to use elevation focus 
Rfocus=40/1000;          %  Elevation focus [m] 
N_elements=65;           %  Number of physical elements 

  
Z=1.480e6;          %  Characteristic acoustic impedance [kg/(m^2 s)] 
Tprf=1/5e3;         %  Pulse repetition frequency [s] 
Tp=M/f0;            %  Pulse duration [s] 
 
P0=sqrt(Ispta*2*Z*Tprf/Tp);   %  Calculate the peak pressure 
 
%  Set the attenuation to 5*0.5 dB/cm, 0.5 dB/[MHz cm] around 
%  f0 and use this: 
 
set_field ('att',2.5*100); 
set_field ('Freq_att',0.5*100/1e6); 
set_field ('att_f0',f0); 
set_field ('use_att',0);          %  Set this flag to one when including attenuation
 
%  Set the sampling frequency 
 
set_sampling(fs); 
 
%  Make the aperture for the rectangles 
 
if (elefocus == 0) 
  ape = xdc_linear_array (N_elements, width, element_height, kerf, 10, 10, focus); 
else 
  ape = xdc_focused_array (N_elements, width, element_height, kerf, Rfocus, 10, 10, focus); 
  end 
 
%  Set the excitation of the aperture 
 
excitation=sin(2*pi*f0*(0:1/fs:M/f0)); 
excitation=excitation.*hanning(length(excitation))'; 
xdc_excitation (ape, excitation); 
 
%  Set the impulse response of the aperture 
 
impulse=sin(2*pi*f0*(0:1/fs:1/f0)); 
impulse=impulse.*hanning(length(impulse))'; 
xdc_impulse (ape, excitation); 
 
%  Find the scaling factor from the peak value 
 
point=[0 0 0]/1000; 
zvalues=(2:2:100)/1000; 
index=1; 
I=0; 
disp('Finding calibration...') 
for z=zvalues 
  point(3)=z; 
  [y,t] = calc_hp(ape,point); 
  I(index)=sum(y.*y)/(2*Z)/fs/Tprf; 
  index=index+1; 
  end 
I_factor=Ispta/max(I); 
 
%  Set the correct scale factor 
 
scale_factor=sqrt(I_factor); 
excitation=scale_factor*excitation; 
xdc_excitation (ape, excitation); 
 
%  Make the calculation in elevation 
 
disp('Finding pressure and intensity.. ') 
point=[0 0 0]/1000; 
zvalues=(1:1:100)/1000; 
index=1; 
I=0; 
Ppeak=0; 
for z=zvalues 
  if rem(z*1000,10)==0
    disp(['Calculating at distance ',num2str(z*1000),' mm'])
    end
  point(3)=z; 
  [y,t] = calc_hp(ape,point); 
  I(index)=sum(y.*y)/(2*Z)/fs/Tprf; 
  Ppeak(index)=max(y); 
  index=index+1; 
  end 
Pmean=sqrt(I*2*Z*Tprf/Tp); 
 
%  Plot the calculated response 
 
figure(1)
subplot(211) 
plot(zvalues*1000,I*1000/(100^2)) 
xlabel('Axial distance [mm]') 
ylabel('Intensity: Ispta  [mW/cm^2]') 
if (elefocus == 0) 
  title(['Focus at ',num2str(focus(3)*1000),' mm, No elevation focus (',Itype,')']) 
else 
  title(['Focus at ',num2str(focus(3)*1000),' mm, elevation focus at ',num2str(Rfocus*1000),' mm (',Itype,')']) 
end 
subplot(212) 
plot(zvalues*1000,Ppeak/1e6) 
xlabel('Axial distance [mm]') 
ylabel('Peak pressure [MPa]') 
 
%  Do the calculation for a single element 
 
figure(2)
xdc_apodization (ape, 0, [zeros(1,floor(N_elements/2)) 1 zeros(1,floor(N_elements/2))]);  
Ppeak_single=0; 
Isingle=0; 
index=1; 
zvalues=0; 
z=0.001/1000; 
factor=(10/1000/z)^(1/200); 
for index=1:200 
  point(3)=z; 
  zvalues(index)=z; 
  [y,t] = calc_hp(ape,point); 
  z=z*factor; 
  Isingle(index)=sum(y.*y)/(2*Z)/fs/Tprf; 
  Ppeak_single(index)=max(y); 
  index=index+1; 
  end 
Pmean=sqrt(I*2*Z*Tprf/Tp); 
clf 
plot(zvalues*1000,Ppeak_single/1e3) 
axis([-0.1 max(zvalues)*1000 0 1.2*max(Ppeak_single)/1e3]) 
xlabel('Axial distance [mm]') 
ylabel('Peak pressure [kPa]') 
if (elefocus == 0) 
  title(['Single element. Focus at ',num2str(focus(3)*1000),' mm, No elevation focus (',Itype,')']) 
else 
  title(['Single element. Focus at ',num2str(focus(3)*1000),' mm, elevation focus at ',num2str(Rfocus*1000),' mm (',Itype,')']) 
end 
 
%  Release the aperture 
 
xdc_free(ape);

    case 4
%% 4 %%%%%%%%%%%%%%%%%%%%%% Phantom Cyst OLD NON-FUNCTIONAL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

folder  = 'C:\Users\yoshi\OneDrive\Documents\MATLAB\cyst_phantom_example.tar\QGP0ka';
list    = dir(fullfile(folder, '*.m'));
nFile   = length(list);
success = false(1, nFile);

for k = 1:nFile
  file = list(k).name;
  try
    run(fullfile(folder, file));
    success(k) = true;
  catch
    fprintf('failed: %s\n', file);
  end
end

    case 5
%% 5 %%%%%%%%%%%%%%% Artificial Point Scatterers and Phantom Cyst - SIMPLE %%%%%%%%%%%%%%%%%%

%  Compress the data to show 60 dB of
%  dynamic range for the cyst phantom image
%
%  version 1.3 by Joergen Arendt Jensen, April 1, 1998.
%  version 1.4 by Joergen Arendt Jensen, August 14, 2007.
%          Calibrated 50 dB display made

% make circle cyst phantoms for later simulations
[phantom_positions, phantom_amplitudes] = cyst_phantom(100000);
save pht_data.mat phantom_positions phantom_amplitudes

f0=3.5e6;                 %  Transducer center frequency [Hz]
fs=100e6;                 %  Sampling frequency [Hz]
c=1540;                   %  Speed of sound [m/s]
no_lines=50;              %  Number of lines in image
image_width=40/1000;      %  Size of image sector
d_x=image_width/no_lines; %  Increment for image

%  Read the data and adjust it in time 

path = 'C:\Users\yoshi\OneDrive\Documents\MATLAB\rf_data';

%%

%  Generate the transducer apertures for send and receive

f0=3.5e6;                %  Transducer center frequency [Hz]
fs=100e6;                %  Sampling frequency [Hz]
c=1540;                  %  Speed of sound [m/s]
lambda=c/f0;             %  Wavelength [m]
width=lambda;            %  Width of element
element_height=5/1000;   %  Height of element [m]
kerf=0.05/1000;          %  Kerf [m]
focus=[0 0 70]/1000;     %  Fixed focal point [m]
N_elements=192;          %  Number of physical elements
N_active=64;             %  Number of active elements 

%  Set the sampling frequency

set_sampling(fs);

%  Generate aperture for emission

xmit_aperture = xdc_linear_array (N_elements, width, element_height, kerf, 1, 10,focus);

%  Set the impulse response and excitation of the xmit aperture

impulse_response=sin(2*pi*f0*(0:1/fs:2/f0));
impulse_response=impulse_response.*hanning(max(size(impulse_response)))';
xdc_impulse (xmit_aperture, impulse_response);

excitation=sin(2*pi*f0*(0:1/fs:2/f0));
xdc_excitation (xmit_aperture, excitation);

%  Generate aperture for reception

receive_aperture = xdc_linear_array (N_elements, width, element_height, kerf, 1, 10,focus);

%  Set the impulse response for the receive aperture

xdc_impulse (receive_aperture, impulse_response);

%   Load the computer phantom

if ~exist('pht_data.mat')
  disp('Scatterer positions should be made by the script mk_pht')
  disp('before this script can be run')
  return
else
  load pht_data
  end

%  Set the different focal zones for reception

focal_zones=[30:20:200]'/1000;
Nf=max(size(focal_zones));
focus_times=(focal_zones-10/1000)/1540;
z_focus=60/1000;          %  Transmit focus

%  Set the apodization

apo=hanning(N_active)';

%   Do linear array imaging

no_lines=50;                    %  Number of lines in image
image_width=40/1000;            %  Size of image sector
d_x=image_width/no_lines;       %  Increment for image

% Do imaging line by line

for i= [1:no_lines]

  if (i > 0)
      %  Test if the file for the line exist.
      %  Skip the simulation, if the line exits and
      %  go the next line. Else make the simulation

      file_name=['rf_data/rf_ln',num2str(i),'.mat'];

      if ~exist(file_name)

        %  Save a file to reserve the calculation

        cmd=['save rf_data/rf_ln',num2str(i),'.mat i'];
        eval(cmd);

        disp(['Now making line ',num2str(i)])

        %  The the imaging direction

        x= -image_width/2 +(i-1)*d_x;

        %   Set the focus for this direction with the proper reference point

        xdc_center_focus (xmit_aperture, [x 0 0]);
        xdc_focus (xmit_aperture, 0, [x 0 z_focus]);
        xdc_center_focus (receive_aperture, [x 0 0]);
        xdc_focus (receive_aperture, focus_times, [x*ones(Nf,1), zeros(Nf,1), focal_zones]);

        %  Calculate the apodization 

        N_pre  = round(x/(width+kerf) + N_elements/2 - N_active/2);
        N_post = N_elements - N_pre - N_active;
        apo_vector=[zeros(1,N_pre) apo zeros(1,N_post)];
        xdc_apodization (xmit_aperture, 0, apo_vector);
        xdc_apodization (receive_aperture, 0, apo_vector);

        %   Calculate the received response

        [rf_data, tstart]=calc_scat(xmit_aperture, receive_aperture, phantom_positions, phantom_amplitudes);

        %  Store the result

        cmd=['save rf_data/rf_ln',num2str(i),'.mat rf_data tstart'];
        disp(cmd)
        eval(cmd);
      else
        disp(['Line ',num2str(i),' is being made by another machine.'])
      end
  end
end

%   Free space for apertures

xdc_free (xmit_aperture)
xdc_free (receive_aperture)

disp('You should now run make_image to display the image')

%%
min_sample=0;
for i=1:no_lines

  %  Load the result
  
  %cmd=['load rf_ln',num2str(i),'.mat'];
  cmd = ['load old_rf_data/rf_ln', num2str(i), '.mat'];
  disp(cmd)
  eval(cmd)
  
  %  Find the envelope
  
  rf_env=abs(hilbert([zeros(round(tstart*fs-min_sample),1); rf_data]));
  env(1:max(size(rf_env)),i)=rf_env;
end

%  Do logarithmic compression

D=10;         %  Sampling frequency decimation factor
dB_range=50;  % Dynamic range for display in dB

disp('Finding the envelope')
log_env=env(1:D:max(size(env)),:)/max(max(env));
log_env=20*log10(log_env);
log_env=127/dB_range*(log_env+dB_range);

%  Make an interpolated image

disp('Doing interpolation')
ID=20;
[n,m]=size(log_env);
new_env=zeros(n,m*ID);

for i=1:n
  new_env(i,:)=interp(log_env(i,:),ID);
end
[n,m]=size(new_env);
  
fn=fs/D;
clf
image(((1:(ID*no_lines-1))*d_x/ID-no_lines*d_x/2)*1000,((1:n)/fn+min_sample/fs)*1540/2*1000,new_env)
disp(d_x/ID-no_lines*d_x/2);
xlabel('Lateral distance [mm]')
ylabel('Axial distance [mm]')
colormap(gray(128))
axis('image')
axis([-20 20 35 90])

    case 6
%% 6 %%%%%%%%%%%%%%% Phased Array Imaging Graph %%%%%%%%%%%%%%%%%%%%%%%%%
        
% Example of use of the new Field II program running under Matlab
%
% This example shows how a phased array B-mode system scans an image
%
% This script assumes that the field_init procedure has been called
%
% Example by Joergen Arendt Jensen, Nov. 28, 1995.
% Generate the transducer apertures for send and receive

f0=3e6; % Transducer center frequency [Hz]
fs=100e6; % Sampling frequency [Hz]
c=1540; % Speed of sound [m/s]
lambda=c/f0; % Wavelength
element_height=5/1000; % Height of element [m]
kerf=0.1/1000; % Kerf [m]
focus=[0 0 70]/1000; % Fixed focal point [m]

% Generate aperture for emission
emit_aperture = xdc_linear_array (128, lambda/2, element_height, kerf, 1, 1,focus);

% Set the impulse response and excitation of the emit aperture
impulse_response=sin(2*pi*f0*(0:1/fs:2/f0));
impulse_response = impulse_response.*hanning(max(size(impulse_response)))';
xdc_impulse (emit_aperture, impulse_response);
excitation=sin(2*pi*f0*(0:1/fs:2/f0));
xdc_excitation (emit_aperture, excitation);

% Generate aperture for reception
receive_aperture = xdc_linear_array (128, lambda/2, element_height, kerf, 1, 1,focus);
% Set the impulse response for the receive aperture
xdc_impulse (receive_aperture, impulse_response);

% Do phased array imaging
point_position=[0 0 70]/1000; % Position of the point to be imaged
no_lines=50; % Number of A-lines in image
sector=20 * pi/180; % Size of image sector

d_theta=sector/no_lines; % Increment in angle for 90 deg. image
% Pre-allocate some storage
image_data=zeros(800,no_lines);
theta= -sector/2;
for i=1:no_lines
    % Set the focus for this direction
    xdc_focus (emit_aperture, 0, [70*sin(theta) 0 70*cos(theta)]/1000);
    xdc_focus (receive_aperture, 0, [70*sin(theta) 0 70*cos(theta)]/1000);
    % Calculate the received response
    [v, t1]=calc_scat(emit_aperture, receive_aperture, point_position, 1);
    % Store the result
    image_data(1:max(size(v)),i) = v';
    times(i) = t1;
    % Steer in another angle
    theta = theta + d_theta;
end
% Here the display of the data is inserted
plot(image_data)

    case 7
%% 7 %%%%%%%%%%%%%%% Single Element Transducer %%%%%%%%%%%%%%%%%%%%%

    elements = 100;
    width = 18.5/1000;
    height = 13/1000;
    kerf = 0;
    num_x = 1;
    num_y = 1;
    focus = [0 0 60];
    f0 = 2.5e6; % transducer center freq Hz
    fs = 400e6; % sampling freq Hz
    c = 1540; % speed of sound m/s
    M = 2; % cycles per pulse
    dt = 1.5; % pulse length
    
    Th = xdc_linear_array(elements, width, height, kerf, num_x, num_y, focus); % transmitter
    Rh = xdc_linear_array(elements, width, height, kerf, num_x, num_y, focus); % reciever same as trasmitter
    
    set_sampling(fs);
    
    excitation = sin(2*pi*f0*(0:1/fs:M/f0)); 
    impulse = sin(2*pi*f0*(0:1/fs:1/f0)); 
    
    xdc_excitation(Th, excitation);
    xdc_excitation(Rh, excitation);
    % impulse=impulse.*hanning(length(impulse))'; 
    
    xdc_impulse(Th, impulse);
    xdc_impulse(Rh, impulse);
    
    t_ir = -2/f0:1/fs:2/f0;
    Bw = 0.6; % bandwidth
    %tc = gauspuls('cutoff',f0, Bw,[],-40); 
    %t = -tc : 1e-7 : tc; 
    [yi,yq,ye] = gauspuls(t_ir, f0, Bw); 

    impulse_response = gauspuls(t_ir, f0, Bw);
    
    plot(t_ir, impulse_response);
    xlabel('time(s)');
    ylabel('Pulse');
    title('Impulse Response vs. Time');
    figure
    
    point = [30, 10, 50]/1000;
    [h, t] = calc_h(Th, point);
    plot((0:length(h)-1)/fs+t,h)
    xlabel('time(s)');
    ylabel('Spatial Pulse Response');
    title('Spatial Impulse Response vs. Time');
    figure
    
    % Expansion to Array transducer
    elements = 64;
    width = 18.5/1000;
    height = 13/1000;
    pitch = 29/1000;
    kerf = 0.020/1000;
    Rfocus = 60/1000;
    subX = 5;
    subY = 15;
    num_x = 1;
    num_y = 1;
    focus = [0 0 60];
    f0 = 2.5e6; % transducer center freq Hz
    fs = 400e6; % sampling freq Hz
    c = 1540; % speed of sound m/s
    M = 2; % cycles per pulse
    dt = 1.5; % pulse length
    
    [h, t] = calc_h(Th, point);
    
    show_xdc(Th);
    
    case 8
%% 8 %%%%%%%%%%%%%%% Phantom Cyst Example - MODULATION %%%%%%%%%%%%%%%%%%

%  Compress the data to show 60 dB of
%  dynamic range for the cyst phantom image
%
%  version 1.3 by Joergen Arendt Jensen, April 1, 1998.
%  version 1.4 by Joergen Arendt Jensen, August 14, 2007.
%          Calibrated 50 dB display made

% make circle cyst phantoms for later simulations
[phantom_positions, phantom_amplitudes] = cyst_phantom(10000);
save pht_data2.mat phantom_positions phantom_amplitudes

f0=3.5e6;                 %  Transducer center frequency [Hz]
fs=100e6;                 %  Sampling frequency [Hz]
c=1540;                   %  Speed of sound [m/s]
no_lines=50;              %  Number of lines in image
image_width=40/1000;      %  Size of image sector
d_x=image_width/no_lines; %  Increment for image

%  Read the data and adjust it in time 

path = 'C:\Users\yoshi\OneDrive\Documents\MATLAB\rf_data2';

%%

%  Generate the transducer apertures for send and receive

f0=3.5e6;                %  Transducer center frequency [Hz]
fs=100e6;                %  Sampling frequency [Hz]
c=1540;                  %  Speed of sound [m/s]
lambda=c/f0;             %  Wavelength [m]
width=lambda;            %  Width of element
element_height=5/1000;   %  Height of element [m]
kerf=0.05/1000;          %  Kerf [m]
focus=[0 0 70]/1000;     %  Fixed focal point [m]
N_elements=192;          %  Number of physical elements
N_active=16;             %  Number of active elements 

%  Set the sampling frequency

set_sampling(fs);

%  Generate aperture for emission

xmit_aperture = xdc_linear_array (N_elements, width, element_height, kerf, 1, 10,focus);

%  Set the impulse response and excitation of the xmit aperture

impulse_response=sin(2*pi*f0*(0:1/fs:2/f0));
impulse_response=impulse_response.*hanning(max(size(impulse_response)))';
xdc_impulse (xmit_aperture, impulse_response);

excitation=sin(2*pi*f0*(0:1/fs:2/f0));
xdc_excitation (xmit_aperture, excitation);

%  Generate aperture for reception

receive_aperture = xdc_linear_array (N_elements, width, element_height, kerf, 1, 10,focus);

%  Set the impulse response for the receive aperture

xdc_impulse (receive_aperture, impulse_response);

%   Load the computer phantom

if ~exist('pht_data2.mat')
  disp('Scatterer positions should be made by the script mk_pht')
  disp('before this script can be run')
  return
else
  load pht_data2
  end

%  Set the different focal zones for reception

focal_zones = [30:20:200]'/1000;
Nf=max(size(focal_zones));
focus_times=(focal_zones-10/1000)/1540;
z_focus=60/1000;          %  Transmit focus

%  Set the apodization

apo=hanning(N_active)';

%   Do linear array imaging

%no_lines=50;                    %  Number of lines in image
image_width=40/1000;            %  Size of image sector
d_x=image_width/no_lines;       %  Increment for image

% Do imaging line by line

for i= [1:no_lines]

  if (i > 0)
      %  Test if the file for the line exist.
      %  Skip the simulation, if the line exits and
      %  go the next line. Else make the simulation

      file_name=['rf_data2/rf_ln',num2str(i),'.mat'];

      if ~exist(file_name)

        %  Save a file to reserve the calculation

        cmd=['save rf_data2/rf_ln',num2str(i),'.mat i'];
        eval(cmd);

        disp(['Now making line ',num2str(i)])

        %  The the imaging direction

        x= -image_width/2 +(i-1)*d_x;

        %   Set the focus for this direction with the proper reference point

        xdc_center_focus (xmit_aperture, [x 0 0]);
        xdc_focus (xmit_aperture, 0, [x 0 z_focus]);
        xdc_center_focus (receive_aperture, [x 0 0]);
        xdc_focus (receive_aperture, focus_times, [x*ones(Nf,1), zeros(Nf,1), focal_zones]);

        %  Calculate the apodization 

        N_pre  = round(x/(width+kerf) + N_elements/2 - N_active/2);
        N_post = N_elements - N_pre - N_active;
        apo_vector=[zeros(1,N_pre) apo zeros(1,N_post)];
        xdc_apodization (xmit_aperture, 0, apo_vector);
        xdc_apodization (receive_aperture, 0, apo_vector);

        %   Calculate the received response

        [rf_data, tstart]=calc_scat(xmit_aperture, receive_aperture, phantom_positions, phantom_amplitudes);

        %  Store the result

        cmd=['save rf_data2/rf_ln',num2str(i),'.mat rf_data tstart'];
        disp(cmd)
        eval(cmd);
      else
        disp(['Line ',num2str(i),' is being made by another machine.'])
      end
  end
end

%   Free space for apertures

xdc_free (xmit_aperture)
xdc_free (receive_aperture)

disp('You should now run make_image to display the image')

%%
min_sample=0;
for i=1:no_lines

  %  Load the result
  
  %cmd=['load rf_ln',num2str(i),'.mat'];
  cmd = ['load rf_data2/rf_ln', num2str(i), '.mat'];
  disp(cmd)
  eval(cmd)
  
  %  Find the envelope
  
  % comment this out later
  tstart = 4.0790e-5;
  rf_env=abs(hilbert([zeros(round(tstart*fs-min_sample),1); rf_data]));
  env(1:max(size(rf_env)),i)=rf_env;
end

%  Do logarithmic compression

D=10;         %  Sampling frequency decimation factor
dB_range=50;  % Dynamic range for display in dB

disp('Finding the envelope')
log_env=env(1:D:max(size(env)),:)/max(max(env));
log_env=20*log10(log_env);
log_env=127/dB_range*(log_env+dB_range);

%  Make an interpolated image

disp('Doing interpolation')
ID=20;
[n,m]=size(log_env);
new_env=zeros(n,m*ID);

for i=1:n
  new_env(i,:)=interp(log_env(i,:),ID);
end
[n,m]=size(new_env);

fn=fs/D;
clf
if (d_x/ID-no_lines*d_x/2 <= 0)
    disp('');
    disp('Change Input Parameters due to error');
    fprintf("\n%d\n", d_x/ID-no_lines*d_x/2);
end

image(((1:(ID*no_lines-1))*d_x/ID-no_lines*d_x/2)*1000,((1:n)/fn+min_sample/fs)*1540/2*1000,new_env)
xlabel('Lateral distance [mm]')
ylabel('Axial distance [mm]')
colormap(gray(128))
axis('image')
axis([-20 20 35 90])

    case 9
%% 9 %%%%%%%%%% Flow Data Generation %%%%%%%%%%%%%%%%%%%%%%

        % Example of use of the new Field II program running under Matlab
    %
    % This example shows how flow can simulated
    %
    % This script assumes that the field_init procedure has been called
    %
    % Example by Joergen Arendt Jensen, March 22, 2011.
    % Generate the transducer apertures for send and receive
    f0=3e6; % Transducer center frequency [Hz]
    fs=100e6; % Sampling frequency [Hz]
    c=1540; % Speed of sound [m/s]
    lambda=c/f0; % Wavelength
    element_height=5/1000; % Height of element [m]
    kerf=0.1/1000; % Kerf [m]
    focus=[0 0 70]/1000; % Fixed focal point [m]
    % Generate aperture
    aperture = xdc_linear_array (128, lambda/2, element_height, kerf, 1, 1,focus);
    % Set the impulse response and excitation of the emit aperture
    impulse_response=sin(2*pi*f0*(0:1/fs:2/f0));
    impulse_response=impulse_response.*hanning(max(size(impulse_response)))';
    xdc_impulse (aperture, impulse_response);
    excitation=sin(2*pi*f0*(0:1/fs:8/f0));
    xdc_excitation (aperture, excitation);
    % Set the seed of the random number generator
    randn('seed',sum(100*clock))
    % Initialize the ranges for the scatterers
    % Notice that the coordinates are in meters
    x_range=0.015; % x range for the scatterers [m]
    y_range=0.015; % y range for the scatterers [m]
    z_range=0.015; % z range for the scatterers [m]
    z_offset=0.70; % Offset of the mid-point of the scatterers [m]
    R=0.005; % Radius of blood vessel [m]
    % Set the number of scatterers. It should be roughly
    % 10 scatterers per resolution cell
    c=1540; % Ultrasound propagation velocity [m/s]
    f0=3e6; % Center frequency of transducer [Hz]
    lambda=c/f0;
    N=round(10*x_range/(5*lambda)*y_range/(5*lambda)*z_range/(lambda*2));
    disp([num2str(N),'Scatterers'])
    % Generate the coordinates and amplitude
    % Coordinates are rectangular within the range.
    % The amplitude has a Gaussian distribution.
    x=x_range*(rand(1,N)-0.5);
    y=y_range*(rand(1,N)-0.5);
    z=z_range*(rand(1,N)-0.5);
    % Find which scatterers that lie within the blood vessel
    
    r=(y.^2+z.^2).^0.5;
    within_vessel= (r < R)';
    % Assign an amplitude and a velocity for each scatterer
    v0=0.5; % Largest velocity of scatterers [m/s]
    velocity=v0*(1-(r/R).^2).*within_vessel';
    blood_to_stationary= 0.1; % Ratio between amplitude of blood to stationary tissue
    amp=randn(N,1).*((1-within_vessel) + within_vessel*blood_to_stationary);
    % Calculate a suitable Tprf
    theta=45/180*pi;
    f_max=2*v0*cos(theta)/c*f0;
    fprf = 3*f_max;
    Tprf=1/fprf; % Time between pulse emissions [s]
    Nshoots=128; % Number of shoots
    % Find the response by calling field
    for i=1:Nshoots
        i;
        % Generate the rotated and offset block of sample
        xnew=x*cos(theta)+z*sin(theta);
        znew=z*cos(theta)-x*sin(theta) + z_offset;
        scatterers=[xnew; y; znew;]' ;
        % Calculate the received response
        [v, t1]=calc_scat(aperture, aperture, scatterers, amp);
        % Store the result
        image_data(1:max(size(v)),i)=v;
        times(i) = t1;
        % Propagate the scatterers and alias them
        % to lie within the correct range
        x=x + velocity*Tprf;
        outside_range= (x > x_range/2);
        x=x - x_range*outside_range;
    end
    % Here the display of the data is inserted
    plot(image_data)

    case 10
%% 10 %%%%%%%%%% Transducer with Apodization %%%%%%%%%%%%%%%%%%%%%%%%%%

    % Set initial parameters
    
f0=3e6; % Transducer center frequency [Hz]
fs=100e6; % Sampling frequency [Hz]
c=1540; % Speed of sound [m/s]
lambda=c/f0; % Wavelength [m]
height=5/1000; % Height of element [m]
width=lambda; % Width of element [m]
kerf=width/4; % Distance between transducer elements [m]
N_elements=10; % Number of elements
no_sub_x=4; % Number of sub-divisions in x-direction of elements.
no_sub_y=10; % Number of sub-divisions in y-direction of elements.

focus=[0 0 40]/1000; % Initial electronic focus

% Define the transducer
Th = xdc_linear_array (N_elements, width, height, kerf, ...
no_sub_x, no_sub_y, focus);

% Set the apodization for the individual mathematical elements
element_no=(1:N_elements)';
hann=hanning(no_sub_y)';
apo=ones(N_elements,1)*reshape(ones(no_sub_x,1)*hann, 1, no_sub_x*no_sub_y);

ele_apodization (Th, element_no, apo);
show_xdc(Th);

    case 100
%% 100 %%%%%% Sketch %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Set initial parameters
    R=8/1000; % Radius of transducer
    ele_size=1/1000; % Size of mathematical elements
    % Define the transducer
    Th = xdc_piston (R,ele_size);
    xdc_show(Th);

end

% free Th and Rh if necessary
try 
    xdc_free(Th);
    xdc_free(Rh);
catch
    disp("No transmitters and receivers need to be freed");
end

%% Function Definitions $$$$

% Show the transducer surface in a surface plot
% Calling: show_xdc(Th)
% Argument Th - Transducer handle
% Return: Plot of the transducer surface on the current figure
% Bote this version onlys shows the defined rectangles
% Version 1.1, June 29, 1998, JAJ

function res = show_xdc (Th)
    % Do it for the rectangular elements
    colormap(cool(128));
    data = xdc_get(Th,'rect');
    [N,M]=size(data);
    % Do the actual display
    for i=1:M
        x=[data(11,i), data(20,i); data(14,i), data(17,i)]*1000;
        y=[data(12,i), data(21,i); data(15,i), data(18,i)]*1000;
        z=[data(13,i), data(22,i); data(16,i), data(19,i)]*1000;
        c=data(5,i)*ones(2,2);
        hold on
        surf(x,y,z,c)
    end
    % Put som axis legends on
    Hc = colorbar;
    view(3)
    xlabel('x [mm]')
    ylabel('y [mm]')
    zlabel('z [mm]')
    grid
    axis('image')
    hold off

end

% Artificial Point Scatterers and Phantom Cyst - SIMPLE
function [positions, amp] = cyst_phantom (N)

    x_size = 40/1000; % Width of phantom [m]
    y_size = 10/1000; % Transverse width of phantom [m]
    z_size = 50/1000; % Height of phantom [m]
    z_start = 30/1000; % Start of phantom surface [m];
    
    % Create the general scatterers
    x = (rand (N,1)-0.5)*x_size;
    y = (rand (N,1)-0.5)*y_size;
    z = rand (N,1)*z_size + z_start;
    
    % Generate the amplitudes with a Gaussian distribution
    amp=randn(N,1);
    
    % Make the cyst and set the amplitudes to zero inside
    r=5/1000; % Radius of cyst [m]
    xc=0/1000; % Place of cyst [m]
    zc=25/1000+z_start;
    inside = ( ((x-xc).^2 + (z-zc).^2) < r^2 );
    amp = amp .* (1-inside);
    
    % Place the point scatterers in the phantom
    dz=z_size/10;
    for i= N-9:N
        x(i) = -15/1000;
        y(i) = 0;
        z(i) = z_start + (i-N+9)*dz;
        amp(i) = 100;
    end
    % Return the variables
    positions=[x y z];
end


%% old phantom generating function
%{


function [positions, amp] = cyst_phantom (N, width, t_width, height, start, size)

    x_size = width/1000; %50/1000;   %  Width of phantom [mm]
    y_size = t_width/1000; %10/1000;   %  Transverse width of phantom [mm]
    z_size = height/1000; %60/1000;   %  Height of phantom [mm]
    z_start = start/1000; %30/1000;  %  Start of phantom surface [mm];

    %  Create the general scatterers

    x = (rand (N,1)-0.5)*x_size;
    y = (rand (N,1)-0.5)*y_size;
    z = rand (N,1)*z_size + z_start;

    %  Generate the amplitudes with a Gaussian distribution

    amp=randn(N,1);

    %  Make the cyst and set the amplitudes to zero inside
    %% Dark circles
    
    %  6 mm cyst
    r=6/2/1000;      % 6 Radius of cyst [mm]
    xc= 10/1000;     %  Place of cyst [mm]
    zc=10/1000+z_start;  

    inside = ( ((x-xc).^2 + (z-zc).^2) < r^2);
    amp = amp .* (1-inside); 

    %  5 mm cyst
    r=5/2/1000;      %  Radius of cyst [mm]
    zc=20/1000+z_start;  

    inside = ( ((x-xc).^2 + (z-zc).^2) < r^2);
    amp = amp .* (1-inside); 

    %  4 mm cyst
    r=4/2/1000;      %  Radius of cyst [mm]
    zc=30/1000+z_start;  

    inside = ( ((x-xc).^2 + (z-zc).^2) < r^2);
    amp = amp .* (1-inside); 

    %  3 mm cyst
    r=3/2/1000;      %  Radius of cyst [mm]
    zc=40/1000+z_start;  

    inside = ( ((x-xc).^2 + (z-zc).^2) < r^2);
    amp = amp .* (1-inside); 

    %  2 mm cyst
    r=2/2/1000;      %  Radius of cyst [mm]
    zc=50/1000+z_start;  

    inside = ( ((x-xc).^2 + (z-zc).^2) < r^2);
    amp = amp .* (1-inside); 

    
    
    %  Make the high scattering region and set the amplitudes to 10 times inside
    %% White circles
    
    %  6 mm region
    r= 10/2/1000; %5/2/1000;       %  Radius of cyst [mm]
    xc= -5/1000; %-5    %  Place of cyst [mm]
    zc=50/1000+z_start;  

    inside = ( ((x-xc).^2 + (z-zc).^2) < r^2) ;
    amp = amp .* (1-inside) + 10*amp .* inside; 

    %  5 mm region
    r=4/2/1000;       %  Radius of cyst [mm]
    zc=40/1000+z_start;  

    inside = ( ((x-xc).^2 + (z-zc).^2) < r^2) ;
    amp = amp .* (1-inside) + 10*amp .* inside; 

    
    %  4 mm region
    r=3/2/1000;       %  Radius of cyst [mm]
    zc=30/1000+z_start;  

    inside = ( ((x-xc).^2 + (z-zc).^2) < r^2) ;
    amp = amp .* (1-inside) + 10*amp .* inside; 

    %  3 mm region
    r=2/2/1000;       %  Radius of cyst [mm]
    zc=20/1000+z_start;  

    inside = ( ((x-xc).^2 + (z-zc).^2) < r^2) ;
    amp = amp .* (1-inside) + 10*amp .* inside; 

    %  2 mm region
    r=1/2/1000;       %  Radius of cyst [mm]
    zc=10/1000+z_start;  

    inside = ( ((x-xc).^2 + (z-zc).^2) < r^2) ;
    amp = amp .* (1-inside) + 10*amp .* inside; 
    
    %%  Place the point scatterers in the phantom (point targets)

    for i=N-5:N
      x(i) = -15/1000;
      y(i) = 0;
      z(i) = z_start + (10+5*10)/1000 + (i-N)*10/1000;
      amp(i) = 20;
    end

    %  Return the variables

    positions=[x y z];
end

%}
