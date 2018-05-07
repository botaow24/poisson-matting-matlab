im = imread('im1.jpg');
mask = imread('mask.bmp');
r = matting(im,mask);

 figure(1),imshow(r);