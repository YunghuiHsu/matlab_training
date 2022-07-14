%%
clear, close all
I = imread('pout.tif'); % read
imshow(I); % show

%%
imageinfo('pout.tif'); % missing funtionality in Octave

%%
imtool('pout.tif'); % missing funtionality in Octave

%%
I = imread('pout.tif'); I2 = histeq(I);
subplot(1,4,1); imhist(I);
subplot(1,4,2); imshow(I);
subplot(1,4,3); imshow(I2);
subplot(1,4,4); imhist(I2);

%%
I = imread('rice.png'); 
J = imread('cameraman.tif'); 
K = imadd(I,J); 
subplot(1,3,1); imshow(I); 
subplot(1,3,2); imshow(K);
subplot(1,3,3); imshow(J); 

%%
I=imread('rice.png'); 
subplot(1,2,1); imshow(I);
J=immultiply(I, 1.5);
subplot(1,2,2); imshow(J);

%%
I = imread('rice.png'); 
subplot(1,2,1); imshow(I); 
J = imrotate(I, 35, 'bilinear'); 
subplot(1,2,2); imshow(J);
size(I)
size(J) 

%%
I = imread('rice.png'); 
imhist(I);

%%
I = imread('rice.png'); 
level=graythresh(I); 
bw=im2bw(I, level);
round(level*255)
subplot(1,2,1); imshow(I); 
subplot (1,2,2); imshow(bw);

%%
I = imread('rice.png'); 
BG = imopen(I, strel('disk', 15));
% BG = imopen(I, strel('ball', 15, 15)); % for Octave-online
imshow(BG);

%%
I = imread('rice.png'); 
subplot(1,3,1); imshow(I); 
BG = imopen(I, strel('disk', 15));
% BG = imopen(I, strel('ball', 15, 15)); % for Octave-online
subplot(1,3,2); imshow(BG);
I2 = imsubtract(I, BG);
subplot(1,3,3); imshow(I2);

%%
I = imread('rice.png'); 
level = graythresh(I); bw = im2bw(I, level); 
subplot (1,2,1); imshow(bw); 

BG = imopen(I, strel('disk', 15));
% BG = imopen(I, strel('ball', 15, 15)); % for Octave-online
I2 = imsubtract(I, BG); level = graythresh(I2); 
bw2 = im2bw(I2, level); 
subplot(1,2,2); imshow(bw2);

%%
I = imread('rice.png');
BG = imopen(I, strel('disk', 15));
% BG = imopen(I, strel('ball', 15, 15)); % for Octave-online
I2 = imsubtract(I, BG); level = graythresh(I2); 
BW = im2bw(I2, level); 
[labeled, numObjects] = bwlabel(BW, 8);

grains = zeros(numObjects, 1);
for ii = 1: numObjects
    grains(ii) = sum(labeled(:) == ii);
end

figure;
set(gcf, 'Color', [1 1 1])
subplot(1, 2, 1);
imshow(BW);
subplot(1, 2, 2);
histogram(grains);
% hist(grains); % for Octave-online

%%
I = imread('rice.png'); 
BG = imopen(I, strel('disk', 15));
% BG = imopen(I, strel('ball', 15, 15)); % for Octave-online
I2 = imsubtract(I, BG); level = graythresh(I2); 
BW = im2bw(I2, level); 
[labeled, numObjects] = bwlabel(BW, 8);
graindata = regionprops(labeled, 'basic');
graindata(51)

%%

%%
originalBW = imread('circles.png'); 
figure;
set(gcf, 'Color', [1 1 1])
subplot( 1, 3, 1); imshow(originalBW);
erodedBW = imerode(originalBW, strel('disk', 7));
%erodedBW = imerode(originalBW, strel('square', 7)); % for Octave-online

subplot( 1, 3, 2); imshow(erodedBW);
erodedBW = imerode(originalBW, strel('disk', 11));
%erodedBW = imerode(originalBW, strel('square', 11)); % for Octave-online
subplot( 1, 3, 3); imshow(erodedBW);

%%
bw = imread('text.png');
subplot( 1, 3, 1); imshow(bw);
dilatedBW = imdilate(bw, strel('line', 11, 90));
subplot( 1, 3, 2); imshow(dilatedBW);
dilatedBW = imdilate(bw, strel('line', 11, 0));
subplot( 1, 3, 3); imshow(dilatedBW);

%%
I = imread('rice.png');
BG = imopen(I, strel('disk', 15));
% BG = imopen(I, strel('ball', 15, 15)); % for Octave-online
I2 = imsubtract(I, BG); level = graythresh(I2); 
BW = im2bw(I2, level); 
figure;
set(gcf, 'Color', [1 1 1])
subplot( 1, 2, 1); imshow(BW);

J = imopen(BW, strel('diamond',2));
subplot( 1, 2, 2); imshow(J);

%%
I = imread('eight.tif'); 
subplot( 1, 3, 1); imshow(I);
J = imnoise(I,'salt & pepper',0.02);
subplot( 1, 3, 2); imshow(J);
K = medfilt2(J); 
subplot( 1, 3, 3); imshow(K);

%%
J = imnoise(imread('eight.tif'),'salt & pepper',0.02); 
subplot( 2, 3, 1); imshow(J);
for i=1:5
    J = uint8(filter2(fspecial('average', 3), J));
    subplot( 2, 3, i+1); imshow(J);
end

%%
I = imread('cameraman.tif');
subplot(1, 4, 1); imshow(I);
J = uint8(filter2(fspecial('sobel'), I));
subplot(1, 4, 2); imshow(J);
K = uint8(filter2(fspecial('sobel')', I));
subplot(1, 4, 3); imshow(K);
L = J + K; subplot(1, 4, 4); imshow(L);
