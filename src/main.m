clear all; clc;
%% Testing SinglePore
dir_src = '..\data\SinglePore';
extension = '.bmp'; %In this case optional because the folder has only
% images with same extension
fprintf("Running code 'poredistribution_yang_optimized' for dataset '%s'\n",dir_src);
C_sp = load_volume(dir_src,extension,'fit');
[C0_opt0, C1_opt0, Re_opt0] = poredistribution_yang_optimized(C_sp);

figure;
slice = round(size(C_sp,3)/2);
subplot(2,2,1); imshow(C_sp(:,:,slice)); title('C Image');
subplot(2,2,2); plot(Re_opt0); title('Pore Distribution');
subplot(2,2,3); imshow(C0_opt0(:,:,slice), []); title('C0');
subplot(2,2,4); imshow(C1_opt0(:,:,slice), []); title('C1');

%% Testing CT_01
dir_src = '..\data\CT_01';
%extension = '.tif';
fprintf("Running code 'poredistribution_yang_optimized' for dataset '%s'\n",dir_src);

C_ct1 = load_volume(dir_src, extension,'fit');
% Invert the pore logic (see image, the pore are black (false/0 value) but must be white (true/1 value)
C_ct1 = ~C_ct1;
[C0_opt1, C1_opt1, Re_opt1] = poredistribution_yang_optimized(C_ct1);

figure;
slice = round(size(C_ct1,3)/2);
subplot(2,2,1); imshow(C_ct1(:,:,slice)); title('C Image');
subplot(2,2,2); plot(Re_opt1); title('Pore Distribution');
subplot(2,2,3); imshow(C0_opt1(:,:,slice), []); title('C0');
subplot(2,2,4); imshow(C1_opt1(:,:,slice), []); title('C1');

%% Testing CT_02
dir_src = '..\data\CT_02';
%extension = '.tif';
fprintf("Running code 'poredistribution_yang_optimized' for dataset '%s'\n",dir_src);
C_ct2 = load_volume(dir_src, extension,'fit'); 
% Invert the pore logic (see image, the pore are black (false/0 value) but must be white (true/1 value)
C_ct2 = ~C_ct2;  
[C0_opt2, C1_opt2, Re_opt2] = poredistribution_yang_optimized(C_ct2);

figure;
slice = round(size(C_ct2,3)/2);
subplot(2,2,1); imshow(C_ct2(:,:,slice)); title('C Image');
subplot(2,2,2); plot(Re_opt2); title('Pore Distribution');
subplot(2,2,3); imshow(C0_opt2(:,:,slice), []); title('C0');
subplot(2,2,4); imshow(C1_opt2(:,:,slice), []); title('C1');