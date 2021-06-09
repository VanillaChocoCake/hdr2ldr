clc;
clear;
hdr = hdrread("hdr/snowman.hdr");
imwrite(hdr, "hdr.png");
delta = 1e-6;
gamma = 1 / 1.5;
Lw = 0.27 .* hdr(:, :, 1) + 0.67 .* hdr(:, :, 2) + 0.06 .* hdr(:, : ,3);
Red = hdr(:,:,1) ./ Lw;
Green = hdr(:,:,2) ./ Lw;
Blue = hdr(:,:,3) ./ Lw;
[height, width] = size(Lw);
Lw = Lw + delta;
average_Lw = exp(sum(sum(log(Lw))) ./ (height .* width));
a = 0.1;
L = a .* Lw ./ average_Lw;
ratio = 1.2;
alpha_1 = 0.35;
alpha_2 = 1.6 * alpha_1;
epsilon = 0.05;
Phi = 8.0;
Scale = 8;
V1 = zeros(height, width, Scale+1);
V2 = V1;
for  i = 1:Scale+1
    s = ratio ^ (i-1);
    sigma_1 = alpha_1 * s / sqrt(2);
    sigma_2 = alpha_2 * s / sqrt(2);
    V1(:, :, i) = imgaussfilt(L, sigma_1);
    V2(:, :, i) = imgaussfilt(L, sigma_2);
end
V = zeros(height, width, Scale+1);
for i = 1:Scale+1
    s = ratio^(i-1);
    V(:, :, i) = (V1(:, :, i) - V2(:, :, i)) ./ ((2 ^ Phi) * a ./ (s ^ 2) + V1(:, :, i));
end
Vsm = V(:, :, Scale);
for k = 1:Scale+1
    for j = 1:height
        for i = 1:width
            if(V(j, i, k) > epsilon)
                Vsm(j, i) = V1(j, i, k);
                break;
            end
        end
    end
end
Ld = L ./ (1 + Vsm);
Ld(Ld > 1) = 1;
Ld(Ld < 0) = 0;
imwrite(Ld, "Ld.png");
res = zeros(height, width, 3);
res(:, :, 1) = Ld .* Red;
res(:, :, 2) = Ld .* Green;
res(:, :, 3) = Ld .* Blue;
res = res .^ gamma;
figure;
imshow(res, []);
imwrite(res, "res_belgium.png");


    