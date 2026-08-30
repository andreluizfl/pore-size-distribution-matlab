clear all; clc;

dir_src = '..\data\CT_01';
extension = '.bmp';

C = load_volume(dir_src,extension,[],-1,false);

sizes = 10:10:100;
repeats = 5;

ts = zeros([length(sizes),2]);
check_res = false([length(sizes),6]);

i = 1;
for N = sizes
    C_sub = remap_volume(C, N, 'top-left-front');
    t_ori=0;t_opt=0;
    for r = 1:repeats
        tic;
        [C0_ori, C1_ori, Re_ori] = poredistribution_yang_original(C_sub);
        t_ori = t_ori+toc;

        tic;
        [C0_opt, C1_opt, Re_opt] = poredistribution_yang_optimized(C_sub);
        t_opt = t_opt+toc;

        check_res(i,:) = [isequal(C0_ori,C0_opt), isequal(C1_ori,C1_opt), isequal(Re_ori,Re_opt), isequal(C0_ori,C0_opt), isequal(C1_ori,C1_opt), isequal(Re_ori,Re_opt)];

    end
    ts(i,:) = [t_ori, t_opt]/repeats;
    i=i+1;
end

A = cat(2,sizes', ts);
T = array2table(A);
T.Properties.VariableNames(1:3) = {'n','original','opt_seq'};
writetable(T,'..\results\ts.csv')

fprintf("Checking if all results are identical to the original: %s\n",string(all(check_res(:))));

figure; grid on; hold on;
plot(ts(:,1),'-r'); 
plot(ts(:,2),'-g');
title('Time Complexity between 2 algorithms');
ylabel('Time (s)');
xlabel('N');
legend('Original','Optimized');
saveas(gcf,'..\results\time_complexity.png');

slice = sizes(end)/2;
figure;
subplot(2,2,1); imshow(C(:,:,slice)); title('C Image');
subplot(2,2,2); plot(Re_ori); title('Pore Distribution');
subplot(2,2,3); imshow(C0_ori(:,:,slice), []); title('C0');
subplot(2,2,4); imshow(C1_ori(:,:,slice), []); title('C1');
saveas(gcf,'..\results\results_original_alg.png');