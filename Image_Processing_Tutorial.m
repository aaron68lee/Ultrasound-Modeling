clc; clear all; close all;

%% Image Processing in MATLAB %%%%

path = "C:\Users\yoshi\OneDrive\Documents\MATLAB\ImageExamples\Ganyu.jpg";
image = imread(path);

% display original image
imshow(image)
title('Original Image');
figure

% display grayscale image
grayImage = rgb2gray(image);
imshow(grayImage);
title('Grayscale Image');
figure


% display binary image
meanVal = mean(grayImage,"all");
binaryImage = grayImage >= meanVal;
title('Binary Image');
imshow(binaryImage)

