clc; clear all; close all;
path(path, '/Users/aaron68lee/Documents/MATLAB/Field_II_ver_3_30_mac')
field_init(0); % field_init(0) suppresses output

rfRemake = 0; % toggle set to 1 if want to overwrite existing RF Data

%% Sample Display Transducer
%{
f0=3e6; % Transducer center frequency [Hz]
fs=100e6; % Sampling frequency [Hz]
c=1540; % Speed of sound [m/s]
lambda=c/f0; % Wavelength [m]
height=5/1000; % Height of element [m]
width=lambda; % Width of element [m]
kerf=width/4; % Distance between transducer elements [m]
N_elements=10; % Number of elements
no_sub_x=6; % Number of sub-divisions in x-direction of elements.
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
%show_xdc(Th)
%}

%% Iterate through all combinations of point-line-elements

points = [100000]; %[1e5/125]; %[10, 100, 1000];
lines = [10, 20, 30, 40, 50]; % default 5th element
elements = [16, 32, 64, 128, 192];
active = [4, 16, 24, 48, 64]; %64];

% 2nd element default value
widthCyst = [50, 50, 75]; 
t_width = [10, 10, 30];
height = [70, 60, 90];
start = [10, 30, 45];
startX = -15; % location of first Cysts X pos in grid

% for generate grid cyst parameters grid dimensions (cols, rows)
numX = 4;
numY = 6;

%{
widthCyst = [50, 50, 75];
t_width = [10, 10, 30];
height = [60, 60, 90];
start = [30, 30, 45];
%}

sizeCyst = [2, 6, 100];
f = [3.5e6]; %, 3.5e6]; %5e6, 20e6];

fprintf("\n");

count = 0;

% current iterations: 1 x 2 x 1 x 1 x 3 x 1 x 1 x 1
% == 6

for freq = 1:length(f) % center freq
    for i = length(points) % points
        for j = 5%length(lines) % lines
            for k = 4 %length(elements) % transducer elements
                for l = 1 %1:length(widthCyst) % cyst width
                    for mm = 1%length(t_width) % cyst transverse width
                        for nn = 1%length(height) % cyst height
                            for st = 1%length(start) % cyst start point Z
                                for ss = 2 %1:length(sizeCyst) % cyst radius size
                                    
                                    count = count + 1;
                                    disp(["Iteration: ", num2str(count)]);
                                    
                                    %% Make Phantom Cysts for Simulations
                                    %  Make the scatteres for a simulation and store
                                    %  it in a file for later simulation use
                                    numPoints = points(i);
                                    % N, width, t_width, height, start, size
                                    [phantom_positions, phantom_amplitudes] = cyst_phantom(numPoints, widthCyst(l), t_width(mm), height(nn), start(st), sizeCyst(ss));
                                    save pht_data.mat phantom_positions phantom_amplitudes

                                    %%  Generate the transducer apertures for send and receive

                                    f0= f(freq);   %3.5e6;                %  Transducer center frequency [Hz]
                                    fs=100e6;                %  Sampling frequency [Hz]
                                    c=1540;                  %  Speed of sound [m/s]
                                    lambda=c/f0;             %  Wavelength [m]
                                    width=lambda;            %  Width of element
                                    element_height=5/1000;   %  Height of element [m]
                                    kerf=0.05/1000;          %  Kerf [m]
                                    ff = 20; %widthCyst(l) + start(st) - t_width(mm);
                                    focus = [0 0 ff]/1000; %; [0 0 70]/1000     %  Fixed focal point [m]
                                    N_elements = elements(k);          %  Number of physical elements
                                    N_active = active(k);             %  Number of active elements 

                                    %  Set the sampling frequency

                                    set_sampling(fs);

                                    %  Generate aperture for emission

                                    xmit_aperture = xdc_linear_array (N_elements, width, element_height, kerf, 1, 10,focus);
                                    %show_xdc(xmit_aperture);
                                    %figure

                                    %  Set the impulse response and excitation of the xmit aperture

                                    impulse_response=sin(2*pi*f0*(0:1/fs:2/f0));
                                    impulse_response=impulse_response.*hanning(max(size(impulse_response)))';
                                    xdc_impulse (xmit_aperture, impulse_response);

                                    excitation=sin(2*pi*f0*(0:1/fs:2/f0));
                                    xdc_excitation (xmit_aperture, excitation);

                                    %  Generate aperture for reception

                                    receive_aperture = xdc_linear_array (N_elements, width, element_height, kerf, 1, 10,focus);
                                    %show_xdc(receive_aperture);
                                    %figure

                                    %  Set the impulse response for the receive aperture

                                    xdc_impulse (receive_aperture, impulse_response);

                                    %%   Load the computer phantom

                                    if ~exist('pht_data.mat')
                                      disp('Scatterer positions should be made by the script mk_pht')
                                      disp('before this script can be run')
                                      return
                                    else
                                      load pht_data
                                      end

                                    %  Set the different focal zones for reception

                                    focal_zones = [30:10:200]'/1000; % change step size to 10, default: 20
                                    Nf=max(size(focal_zones));
                                    focus_times=(focal_zones - 5/1000)/1540; % default -10
                                    z_focus = ff; %height(nn); %30; %60/1000;   %  Transmit focus

                                    %  Set the apodization

                                    apo=hanning(N_active)';

                                    %   Do linear array imaging

                                    no_lines = lines(j);    %  Number of lines in image
                                    image_width = widthCyst(l)/1000;  % max width is 60 for default transducer; Size of image sector
                                    d_x = image_width/no_lines;       %  Increment for image

                                    % Do imaging line by line

                                    %% RF Data Generation
                                    
                                    if(rfRemake) % only remake lines if toggle true

                                        no_lines = lines(j);
                                        for I=[1:no_lines]

                                          %  Test if the file for the line exist.
                                          %  Skip the simulation, if the line exits and
                                          %  go the next line. Else make the simulation

                                          file_name=['rf_data/rf_ln',num2str(I),'.mat'];

                                          if ~exist(file_name) || exist(file_name)

                                            %  Save a file to reserve the calculation

                                            cmd=['save rf_data/rf_ln',num2str(I),'.mat i'];
                                            eval(cmd);

                                            disp(['Now making line ',num2str(I)])

                                            %  The the imaging direction

                                            x = -image_width/2 +(I-1)*d_x;

                                            %   Set the focus for this direction with the proper reference point

                                            xdc_center_focus (xmit_aperture, [x 0 0]);
                                            xdc_focus (xmit_aperture, 0, [x 0 z_focus]);
                                            xdc_center_focus (receive_aperture, [x 0 0]);
                                            xdc_focus (receive_aperture, focus_times, [x*ones(Nf,1), zeros(Nf,1), focal_zones]);

                                            %%  Calculate the apodization 

                                            % ADD APODIZATION BACK FOR %%%%%%%%
                                            % FREQ_SAMPLE < 3.5MHz &
                                            % Transducer elements < 192

                                            N_pre  = round(x/(width+kerf) + N_elements/2 - N_active/2);
                                            N_post = N_elements - N_pre - N_active;
                                            apo_vector=[zeros(1,N_pre) apo zeros(1,N_post)];
                                            %xdc_apodization (xmit_aperture, 0, apo_vector);
                                            %xdc_apodization (receive_aperture, 0, apo_vector);

                                            %   Calculate the received response

                                            %phantom_positions = phantom_positions(1:10, :);
                                            %phantom_amplitudes = phantom_amplitudes(1:10);

                                            [rf_data, tstart] = calc_scat(xmit_aperture, receive_aperture, phantom_positions, phantom_amplitudes);

                                            %rf_data = rf_data(1:1000);
                                            %  Store the result

                                            cmd=['save rf_data/rf_ln',num2str(I),'.mat rf_data tstart'];
                                            disp(cmd)
                                            eval(cmd);
                                          else
                                            disp(['Line ',num2str(I),' is being made by another machine.'])
                                          end
                                        end

                                        %   Free space for apertures

                                        xdc_free (xmit_aperture)
                                        xdc_free (receive_aperture)

                                        disp('You should now run make_image to display the image')
                                    end
                                    
                                    %% Make Image 
                                    %  Compress the data to show 60 dB of
                                    %  dynamic range for the cyst phantom image
                                    %
                                    %  version 1.3 by Joergen Arendt Jensen, April 1, 1998.
                                    %  version 1.4 by Joergen Arendt Jensen, August 14, 2007.
                                    %          Calibrated 50 dB display made

                                    %f0=3.5e6;                 %  Transducer center frequency [Hz]
                                    %fs=100e6;                 %  Sampling frequency [Hz]
                                    c=1540;                   %  Speed of sound [m/s]
                                    %no_lines=50;              %  Number of lines in image
                                    %image_width = 40/1000;      %  Size of image sector
                                    %d_x = image_width/no_lines; %  Increment for image

                                    %  Read the data and adjust it in time 

                                    min_sample=0;
                                    for I=1:no_lines

                                      %  Load the result

                                      cmd=['load rf_data/rf_ln',num2str(I),'.mat'];
                                      disp(cmd)
                                      eval(cmd)

                                      %  Find the envelope

                                      rf_env=abs(hilbert([zeros(round(tstart*fs-min_sample),1); rf_data]));
                                      env(1:max(size(rf_env)),I) = rf_env;
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
                                    for I=1:n
                                      new_env(I,:)=interp(log_env(I,:),ID);
                                    end
                                    [n,m]=size(new_env);

                                    fn=fs/D;
                                    clf
                                    
                                    %% Image Plot Starts
                                    
                                    %subplot(1, length(lines), length(lines)*(j-1) + k);
                                    %(1:(ID*no_lines-1))*d_x/ID-no_lines*d_x/2)*1000;
                                    hold on
                                    
                                    
                                    locX = [];
                                    locY = [];
                                    Xf = lines(j) * 20;
                                    Yf = length(new_env(:, 1));
                                    inc = Xf / widthCyst(l);
                                    
                                    %% Temp get indices of splices for
                                    % testing
                                    for xi = 1:numX
                                        locX(end + 1) = (startX + widthCyst(l)/2)*inc + 10*inc*(xi-1);
                                        %new_env(:, locX(xi)) = 1000; % vertical lines
                                    end
                                    
                                    space = 130;
                                    ratio = 13;
                                    for yi = 1:numY
                                        locY(end + 1) = space*(yi+1) + 10;
                                        %{
                                        new_env(locY(yi), :) = 1000; % horizontal lines bounds at [140, 140 + space*height(nn)/10]
                                        new_env(locY(yi) - ratio*4, :) = 1000;
                                        new_env(locY(yi) + ratio*4, :) = 1000;
                                        %}
                                        % y bounds at [140, 1050] ratio
                                        % 1:13
                                        % y dist at locY
                                    end
                                    
                                    %% Generate Image
                                    
                                    image(((1:(ID*no_lines-1))*d_x/ID-no_lines*d_x/2)*1000,((1:n)/fn+min_sample/fs)*1540/2*1000,new_env);                       
                                    xlabel('Lateral distance [mm]')
                                    ylabel('Axial distance [mm]')
                                    title(['FreqC: ', num2str(f0), ' Points: ', num2str(points(i)), ', Lines: ', num2str(no_lines), ', Elements: ', num2str(N_elements), ' Cyst Width: ', num2str(widthCyst(l))]);
                                    %title(string);
                                    subtitle(['Figure: ', num2str(count)]); 
                                    colormap(gray(128))
                                    axis('image')
                                    endCoor = height(nn) + start(st);
                                    axis([-widthCyst(l)/2, widthCyst(l)/2, start(st) - 5, endCoor + 5]);
                                    
                                    %% Annotate The Ultrasound Image
                                    
                                    plot(0, ff, 'r+', 'MarkerSize', 15, 'LineWidth', 2);
                                    
                                    for xi = 1:numX
                                        for yi = 1:numY
                                            bound1 = locY(yi) - ratio*xi;
                                            bound2 = locY(yi) + ratio*xi;
                                            plot(locX(xi)/20 - 25, bound1/13, 'g+', 'MarkerSize', 5, 'LineWidth', 2);
                                            plot(locX(xi)/20 - 25, bound2/13, 'b+', 'MarkerSize', 5, 'LineWidth', 2);
                                           
                                            %data2 = rawEnv(bound1-100:bound2+100, locX(xi));
                                        end
                                    end
                                    %plot(index1, halfMax, 'r+', 'MarkerSize', 15, 'LineWidth', 2);
                                    %axis([-20, 20, start(st), endCoor])
                                    %axis([-20 20 35 90])
                                    hold off
                                    
                                    figure
  
                                    hold on
                                    circles(4, 6, height(nn), start(st), new_env, locX, locY); % simple grid representation
                                    hold off
                                    
                                    %% Save Matrix as Image
                                    %{
                                    I = new_env; %single(new_env);
                                    Nan = isnan(I);
                                    Sum1 = sum(Nan, 'all');
                                    INF = isinf(I);
                                    Sum2 = sum(INF, 'all');
                                    I(isnan(I)) = 0;
                                    I(isinf(I)) = 100;
                                    I = I(140:1050, :);
                                    %I = 255*rescale(I);
                                    
                                    I = 255*(I - min(I(:))) ./ (max(I(:)) - min(I(:))); %scale values between 0 and 255
                                    I = cast(I,'uint8');
                                    imshow(I);
                                    imwrite(I,'myImage.png')
                                    
                                    figure
                                    
                                    for times = 1:3
                                        I = imsharpen(I);
                                        I = imadjust(I);
                                        I = histeq(I);
                                    end
                                    
                                    imshow(I);
                                    imwrite(I,'myImage2.png')
                                    %I = rgb2gray(I);
                                    
                                    [centers,radii] = imfindcircles(I,[15 55],'ObjectPolarity','dark', ...
                                    'Sensitivity',0.86,'Method','twostage');
                                
                                    viscircles(centers, radii);
                                    
                                    hold on
                                    %}
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

%% Generate Labels
% overlay FWHM data over circles in given current image

function genLabels()
   
    numX = 4;
    numY = 6;
    for i = 1:numX
        for j = 1:numY
            
        end
    end
    
end

%% Show Transducer Element


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

%% Full Width Half Max
function Fwhm = fwhm(data)
        % Find the half max value.
        
    error = 10; % correction for noise regions
    halfMax = (min(data) + max(data)) / 2;
    
    % Find where the data first drops below half the max.
    %index1 = find(data >= halfMax, 1, 'first');
    index1 = find(data <= halfMax, 1, 'first'); % consider using data(error:end)
    if(length(index1) == 0)
        index1 = 0;
    end
    % Find where the data last rises above half the max.
    %index2 = find(data >= halfMax, 1, 'last'); %index2 = find(data <= halfMax, 1, 'last'); %
    section = data(index1:end);
    index2 = find(section >= halfMax, 1, 'last');
    if(length(index2) == 0)
        index2 = length(data);
    end
    
    Fwhm = index2-index1 + 1; % FWHM in indexes.
    
    %{
    hold on
    plot(index1, halfMax, 'r+', 'MarkerSize', 15, 'LineWidth', 2);
    plot(index2, halfMax, 'r+', 'MarkerSize', 15, 'LineWidth', 2);
    plot(index1, data(index1), 'b+', 'MarkerSize', 10, 'LineWidth', 2); % actual locations
    plot(index2, data(index2), 'b+', 'MarkerSize', 10, 'LineWidth', 2);
    text(Fwhm/2 + index1, halfMax, ['Fwhm: ' num2str(Fwhm)]);
    hold off
    %}
    
    % OR, if you have an x vector
    %fwhmx = x(index2) - x(index1);
end

%% Phantom generating function batch resolving

function [positions, amp] = cyst_phantom (N, width, t_width, height, start, size)

    x_size = width/1000; %50/1000;   %  Width of phantom [mm]
    y_size = t_width/1000; %10/1000;   %  Transverse width of phantom [mm]
    z_size = height/1000; %60/1000;   %  Height of phantom [mm]
    z_start = start/1000; %30/1000;  %  Start of phantom surface [mm];

    %  Create the general scatterers

    x = (rand (N,1)-0.5)*x_size;
    y = (rand (N,1)-0.5)*y_size;
    z = rand (N,1)*z_size + z_start;

    Zf = height + start;
    numPointScat = 5;
    increment = height / numPointScat;
    
    %  Generate the amplitudes with a Gaussian distribution

    amp=randn(N,1);

    %% Dark Circles
    %{
    %  6 mm cyst
    r=6/2/1000;      % 6 Radius of cyst [mm]
    xc= 10/1000;     %  Place of cyst [mm]
    zc=1*increment/1000+z_start;  

    inside = ( ((x-xc).^2 + (z-zc).^2) < r^2);
    amp = amp .* (1-inside); 

    %  5 mm cyst
    r=5/2/1000;      %  Radius of cyst [mm]
    zc=2*increment/1000+z_start;  

    inside = ( ((x-xc).^2 + (z-zc).^2) < r^2);
    amp = amp .* (1-inside); 

    %  4 mm cyst
    r=4/2/1000;      %  Radius of cyst [mm]
    zc=3*increment/1000+z_start;  

    inside = ( ((x-xc).^2 + (z-zc).^2) < r^2);
    amp = amp .* (1-inside); 

    %  3 mm cyst
    r=3/2/1000;      %  Radius of cyst [mm]
    zc=4*increment/1000+z_start;  

    inside = ( ((x-xc).^2 + (z-zc).^2) < r^2);
    amp = amp .* (1-inside); 

    %  2 mm cyst
    r=2/2/1000;      %  Radius of cyst [mm]
    zc=5*increment/1000+z_start;  

    inside = ( ((x-xc).^2 + (z-zc).^2) < r^2);
    amp = amp .* (1-inside); 
    %}
    
    %% Old code white ===========================
    %{
    %  6 mm region
    r= 5/2/1000; %5/2/1000;       %  Radius of cyst [mm]
    xc= -5/1000; %-5    %  Place of cyst [mm]
    zc=5*increment/1000+z_start;  

    inside = ( ((x-xc).^2 + (z-zc).^2) < r^2) ;
    amp = amp .* (1-inside) + 10*amp .* inside; 

    %  5 mm region
    r=4/2/1000;       %  Radius of cyst [mm]
    zc=4*increment/1000+z_start;  

    inside = ( ((x-xc).^2 + (z-zc).^2) < r^2) ;
    amp = amp .* (1-inside) + 10*amp .* inside; 

    
    %  4 mm region
    r=3/2/1000;       %  Radius of cyst [mm]
    zc=3*increment/1000+z_start;  

    inside = ( ((x-xc).^2 + (z-zc).^2) < r^2) ;
    amp = amp .* (1-inside) + 10*amp .* inside; 

    %  3 mm region
    r=2/2/1000;       %  Radius of cyst [mm]
    zc=2*increment/1000+z_start;  

    inside = ( ((x-xc).^2 + (z-zc).^2) < r^2) ;
    amp = amp .* (1-inside) + 10*amp .* inside; 

    %  2 mm region
    r=1/2/1000;       %  Radius of cyst [mm]
    zc=1*increment/1000+z_start;  

    inside = ( ((x-xc).^2 + (z-zc).^2) < r^2) ;
    amp = amp .* (1-inside) + 10*amp .* inside; 
    %}
    
    %% Grid circles
    
    %  Make the high scattering region and set the amplitudes to 10 times inside
    
    numX = 4;
    numY = 6;
    radii = [2, 4, 6, 8];
    Xstart = -15;
    increment = (height - start) / numY;
    
    for i = 1:numX
        for j = 1:numY
            % start from top to bottom
            r = radii(i)/2/1000; %5/2/1000;       %  Radius of cyst [mm]
            xc = (Xstart + 10*(i-1))/1000; %-5    %  Place of cyst [mm] incremented by 10 per iteration
            zc = (j)*increment/1000 + z_start;  % units of "increment" determined by bounds

            inside = ( ((x-xc).^2 + (z-zc).^2) < r^2) ;
            %amp = amp .* (1-inside) + 10*amp .* inside; % white circles
            amp = amp .* (1-inside); % dark circles
            
        end
    end
    
    
    %%  Place the point scatterers in the phantom (point targets)

    % old code modified auto-scaling
    
    for i = N - numPointScat:N
      x(i) = -20/1000;
      y(i) = 0;
      z(i) = z_start + height/1000 + (i-N)*increment / 1000; %(i-N)*10/1000; (10+5*10)/1000
      amp(i) = 20;
    end
    
    
    % make single point scatterer in corner to normalize color map
    %{
    for i=N-1:N
      x(i) = -20/1000;
      y(i) = 0;
      z(i) = 15; %z_start + (10+5*10)/1000 + (i-N)*10/1000;
      amp(i) = 20;
    end
    %}

    %  Return the variables 

    positions=[x y z];
end

%% Creates a grid of circles with FWHM Overlay
function circles(numX, numY, height, start, rawEnv, locX, locY)

    xc = [];
    zc = [];
    r = [];
    centers = [];
    
    numX = 4;
    numY = 6;
    ratio = 13;
    
    radii = [2, 4, 6, 8];
    Xstart = -15;
    increment = (height - start) / numY;
    
    for i = 1:numX
        for j = 1:numY
            % start from top to bottom
            r(end + 1) = radii(i)/2/1000; %5/2/1000;       %  Radius of cyst [mm]
            xc(end + 1) = (Xstart + 10*(i-1))/1000; %-5    %  Place of cyst [mm] incremented by 10 per iteration
            zc(end + 1) = (j-1)*increment/1000 + start/1000;  % units of "increment" determined by bounds

            %plot(xc(i), 50/1000, 'b+', 'MarkerSize', 30, 'LineWidth', 2);
            %centers((i-1)*numX + j*numY, 1) = xc(end);
            %centers((i-1)*numX + j*numY, 2) = zc(end);
        end
    end
   
    
    centers(:, 1) = xc;
    centers(:, 2) = zc;
    viscircles(centers, r);
    offset = 2/1000;

    %% generate FWHM overlay on sample grid
    samples = 3; % number of vertical lines to call FWHM on
    r = []; % samples long array storing results of FWHM
    variation = [-5, 0, 5];
    radii = [1, 2, 3, 4];
    
    for xi = 1:numX
        for yi = 1:numY
            bound1 = locY(yi) - ratio*radii(xi);
            bound2 = locY(yi) + ratio*radii(xi);
            r = [];
            for k = 1:samples
                data = rawEnv(bound1:bound2, locX(xi) + variation(k));
                data2 = rawEnv(bound1-100:bound2+100, locX(xi) + variation(k));
                result = fwhm(data);
                r(end+1) = result/13;
                %value = ['' num2str(round(fwhm(data)/13, 3))];
            end
            value = ['' num2str(round(mean(r), 3))];
            t = text(xc((xi-1)*6 + yi) + offset, zc((xi-1)*6 + yi), value);
            t.Color = [0 0 0]; % RGB in range (0, 1)
        end
    end
    
    % Experimental
    %{
    offset = 2/1000;
    for i = 1:24
        data = rawEnv(799, 868:972)
        value = ['' num2str(round(fwhm(data), 3))];
        t = text(xc(i) + offset, zc(i), value);
        t.Color = [0 0 0]; % RGB in range (0, 1)
    end
    %}
    
    
    %plot(xc, zc, 'b+', 'MarkerSize', 30, 'LineWidth', 2);
    grid on;
    axis equal

end

%% Create a Circle at a Location

function circleCoor(x, y, r)

    xc = [];
    yc = [];
    radius = [];
    centers = [];
    
    radius(end + 1) = r;
    xc(end + 1) = x;
    yc(end + 1) = y;
    
    centers(:, 1) = xc;
    centers(:, 2) = yc;
    viscircles(centers, radius);
    axis equal

end

%% Model Cyst Phantom

%  Create a computer model of a cyst phantom. The phantom contains
%  fiven point targets and 6, 5, 4, 3, 2 mm diameter waterfilled cysts, 
%  and 6, 5, 4, 3, 2 mm diameter high scattering regions. All scatterers 
%  are situated in a box of (x,y,z)=(50,10,60) mm and the box starts 
%  30 mm from the transducer surface.
%
%  Calling: [positions, amp] = cyst_phantom (N);
%
%  Parameters:  N - Number of scatterers in the phantom
%
%  Output:      positions  - Positions of the scatterers.
%               amp        - amplitude of the scatterers.
%
%  Version 2.2, April 2, 1998 by Joergen Arendt Jensen

% enter parameters in [mm]
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

    Zf = height + start;
    numPointScat = 5;
    increment = (Zf - 35) / numPointScat;
    
    %  Generate the amplitudes with a Gaussian distribution

    amp=randn(N,1);

    
    %  Make the cyst and set the amplitudes to zero inside
    %% Dark circles
    
    %  6 mm cyst
    r=6/2/1000;      % 6 Radius of cyst [mm]
    xc= 10/1000;     %  Place of cyst [mm]
    zc=1*increment/1000+z_start;  

    inside = ( ((x-xc).^2 + (z-zc).^2) < r^2);
    amp = amp .* (1-inside); 

    
    
    %  5 mm cyst
    r=5/2/1000;      %  Radius of cyst [mm]
    zc=2*increment/1000+z_start;  

    inside = ( ((x-xc).^2 + (z-zc).^2) < r^2);
    amp = amp .* (1-inside); 

    %  4 mm cyst
    r=4/2/1000;      %  Radius of cyst [mm]
    zc=3*increment/1000+z_start;  

    inside = ( ((x-xc).^2 + (z-zc).^2) < r^2);
    amp = amp .* (1-inside); 

    %  3 mm cyst
    r=3/2/1000;      %  Radius of cyst [mm]
    zc=4*increment/1000+z_start;  

    inside = ( ((x-xc).^2 + (z-zc).^2) < r^2);
    amp = amp .* (1-inside); 

    %  2 mm cyst
    r=2/2/1000;      %  Radius of cyst [mm]
    zc=5*increment/1000+z_start;  

    inside = ( ((x-xc).^2 + (z-zc).^2) < r^2);
    amp = amp .* (1-inside); 
    
    %  Make the high scattering region and set the amplitudes to 10 times inside
    %% White circles
    
    %  6 mm region
    r= 5/2/1000; %5/2/1000;       %  Radius of cyst [mm]
    xc= -5/1000; %-5    %  Place of cyst [mm]
    zc=5*increment/1000+z_start;  

    inside = ( ((x-xc).^2 + (z-zc).^2) < r^2) ;
    amp = amp .* (1-inside) + 10*amp .* inside; 

    %  5 mm region
    r=4/2/1000;       %  Radius of cyst [mm]
    zc=4*increment/1000+z_start;  

    inside = ( ((x-xc).^2 + (z-zc).^2) < r^2) ;
    amp = amp .* (1-inside) + 10*amp .* inside; 

    
    %  4 mm region
    r=3/2/1000;       %  Radius of cyst [mm]
    zc=3*increment/1000+z_start;  

    inside = ( ((x-xc).^2 + (z-zc).^2) < r^2) ;
    amp = amp .* (1-inside) + 10*amp .* inside; 

    %  3 mm region
    r=2/2/1000;       %  Radius of cyst [mm]
    zc=2*increment/1000+z_start;  

    inside = ( ((x-xc).^2 + (z-zc).^2) < r^2) ;
    amp = amp .* (1-inside) + 10*amp .* inside; 

    %  2 mm region
    r=1/2/1000;       %  Radius of cyst [mm]
    zc=1*increment/1000+z_start;  

    inside = ( ((x-xc).^2 + (z-zc).^2) < r^2) ;
    amp = amp .* (1-inside) + 10*amp .* inside; 
    
    %%  Place the point scatterers in the phantom (point targets)
 
    for i=N - numPointScat:N
      x(i) = -15/1000;
      y(i) = 0;
      z(i) = z_start + height/1000 + (i-N)*increment / 1000; %(i-N)*10/1000; (10+5*10)/1000
      amp(i) = 20;
    end
    
    
    %  Return the variables

    positions=[x y z];
    end

%}