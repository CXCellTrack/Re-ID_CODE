
basedir = 'E:\人脸识别数据集\ChokePoint\P1E_S1\P1E_S1_C1\';

im_bg = im2double( imread([basedir, '00000000.jpg']) ); % imshow(im_bg)
im_1 = im2double( imread([basedir, '00000210.jpg']) ); 
im_sub = abs(im_1 - im_bg);
im_sub = rgb2gray(im_sub);


threshhold = 0.07;
tmp = im2bw(im_sub, threshhold);
tmp = imopen(tmp, strel('disk',3));

figure
subplot(1,3,1); imshow(im_bg); 
subplot(1,3,2); imshow(im_1); 
subplot(1,3,3); imshow(tmp); 


figure
% pic1 = imread('C:\Users\Administrator\Desktop\230.jpg');
% pic1(pic1<248) = 0;
% pic1 = imfill(pic1, 'holes');
% pic1 = imclose(pic1, strel('disk',3));
% pic1 = bwareaopen(pic1, 50); imshow(pic1);

subplot(3,1,1); imshow(imread([basedir, '00000210.jpg']));
subplot(3,1,2); imshow(imread([basedir, '00000220.jpg'])); 
subplot(3,1,3); imshow(imread([basedir, '00000230.jpg'])); 







