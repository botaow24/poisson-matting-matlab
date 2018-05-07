function [ r ] = matting( im1,mask )
%MATTING Summary of this function goes here
%   Detailed explanation goes here
im1 = im2double(im1);
im1 = rgb2gray(im1);


f_mask = mask > 200;
b_mask = mask < 50;
bt_mask = (~f_mask) & (~b_mask);
[~,f_map] = bwdist(f_mask);
[~,b_map] = bwdist(b_mask);
F = zeros(size(mask));
B = zeros(size(mask));
F(:) = im1(f_map);
B(:) = im1(b_map);
r = solver(F-B,im1,bt_mask,f_mask);
for asd = 1:3
    f_mask ((r > 0.95) & (bt_mask)) = true;
    b_mask ((r < 0.05) & (bt_mask)) = true;
    bt_mask = (~f_mask) & (~b_mask);
    if(any(bt_mask) == false)
        ['break']
        break;
    end
    
    [~,f_map] = bwdist(f_mask);
    [~,b_map] = bwdist(b_mask);
    F(:) = im1(f_map);
    B(:) = im1(b_map);
    r = solver(F-B,im1,bt_mask,f_mask);
end

end

function result = solver(FD,im1,mask,f_mask)
FD = imgaussfilt(FD,3.0);
convMask = [ 0,-1, 0;
            -1, 4,-1
            0,-1, 0];

m = im1 ;

m = conv2(m,convMask,'same');% laplacian
s = sum(mask(:));                            %number of pixels

b = zeros(s,1);
imTable = zeros(size(mask));      % index look up table
[i,j] = find(mask);
index = find(mask);
for k =1 :s                                  % fill the look up table
    imTable(i(k),j(k)) = k;
end
spi = zeros(5*s,1);                          % data for the sparse matrix
spj = zeros(5*s,1);
spvalue = zeros(5*s,1);
count = 1;

for k = 1 : s
    spi(count) = k;
    spj(count) = k;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% diagonal
    spvalue(count) = 4; 
    count = count +1;
    
    i0 = i(k);
    j0 = j(k);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % test neighbor 1
    temp = imTable(i0-1,j0);
    if temp>0                          % not border
        spi(count) = k;
        spj(count) = temp;
        spvalue(count) = -1;
        count = count +1;
    else
        b(k) = b(k) + f_mask(i0-1,j0)* FD(i0-1,j0);     % border
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % test neighbor 2
    temp = imTable(i0+1,j0);
    if temp>0
        spi(count) = k;
        spj(count) = temp;
        spvalue(count) = -1;
        count = count +1;
    else
        b(k) = b(k) + f_mask(i0+1,j0)* FD(i0+1,j0);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % test neighbor 3
    temp = imTable(i0,j0-1);
    if temp>0
        spi(count) = k;
        spj(count) = temp;
        spvalue(count) = -1;
        count = count +1;
    else
        b(k) = b(k) + f_mask(i0,j0-1)* FD(i0,j0-1);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % test neighbor 4
    temp = imTable(i0,j0+1);
    if temp>0
        spi(count) = k;
        spj(count) = temp;
        spvalue(count) = -1;
        count = count +1;
    else
        b(k) = b(k) + f_mask(i0,j0+1) * FD(i0,j0+1);
    end
    
    c1 = m(i0,j0);
    b(k) = b(k) +c1;           % add the laplacian image
    
end

count = count - 1;
A =  sparse(spi(1:count),spj(1:count),spvalue(1:count),s,s);
f = A\b;
% fill the pixel
result = double(f_mask);
for k = 1 :s
    result(index(k)) = f(k) / FD(index(k));
end
end
