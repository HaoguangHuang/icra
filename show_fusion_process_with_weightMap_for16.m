function show_fusion_process_for16
% figure(1) is to show camera coordinate systems

ax_seq = [1 2 3];
fileNum = 10;
lidarRotation = []; %记录前i帧相对global coo的变换
lidarTranslation = [];
% alpha = 0;
% beta = 0;
% gamma = 0;
% R = rodrigues([alpha beta gamma]);
%读取10帧.mat数据
for i=1:fileNum
    a(i) = load(num2str(i, 'cloud_cluster_%d_2.txt.mat'));
    [tmp R t] = transformFormat(a(i).data_out);%tmp是3*n的向量，包含原点
    %对frame(i).data除了原点外的数据，全部旋转（等效于坐标轴旋转）
    if 0
        tmp_rotated = rotate(tmp, R , t);
        frame(i).data = change_axis(tmp_rotated); %frame是元素为.mat的向量
    else 
        frame(i).data = change_axis(tmp);
    end
    lidarRotation = [lidarRotation; R];
    lidarTranslation = [lidarTranslation; t]
end

K = [1000, 0,  0; ...
     0,  100, 0; ...
     0,  0,  1]; 
H = 64; W = 300;
for i = 1:fileNum
    Orig = frame(i).data(ax_seq,1)*1000;   %世界坐标
    Ps3D = frame(i).data(ax_seq,2:end)*1000;
    Ps2D = K*Ps3D; z = Ps3D(3,:);
    Ps2D(1,:) = Ps2D(1,:)./z;
    Ps2D(2,:) = Ps2D(2,:)./z;
%     figure(2), plot(Ps2D(1,:), Ps2D(2,:),'k.'),grid on
    
    us = Ps2D(2,:);  vs = Ps2D(1,:); 
    u_min = min(us); u_max = max(us);
    v_min = min(vs); v_max = max(vs);
    Ps2D_new = Ps2D;
    Ps2D_new(1,:) = Ps2D_new(1,:) - v_min + 10;
    Ps2D_new(2,:) = Ps2D_new(2,:) - u_min + 10;
    D_map = produce_depth(Ps2D_new, H, W);
    imwrite(uint16(D_map), num2str(i,'depth_beforeInterpolate_withoutRotate_%d_2.png'));
    figure(11),imshow(mat2gray(D_map)), colormap(jet)
    
    % compute the fusion weights
    BW = uint8(255*logical(D_map));
    ROI = find_ROI(BW);
    psf = fspecial('averag',7);
    Blur_BW = imfilter(BW,psf,'replicate');
    Blur_BW1 = double( double(Blur_BW) / norm(double(Blur_BW(:)),2) ) * 10;
    figure(10), imshow(BW)
    figure(20), imshow(uint8(mat2gray(Blur_BW1)*255),'border','tight')
    mask = zeros(size(BW)); mask(ROI(1):ROI(2), ROI(3):ROI(4)) = 1;
    Blur_BW1(logical(BW)) = 1;  Blur_BW1(~mask) = 0;
    figure(30), imshow(mat2gray(Blur_BW1))
    imwrite(uint8(mat2gray(Blur_BW1)*255), num2str(i,'weight_beforeInterpolate_withoutRotate_%d_2.png'));
    
    % depth interpolation
    D_in = D_map(ROI(1):ROI(2), ROI(3):ROI(4));
    D_out = depth_interp(D_in);
    
    P0 = mean(Ps3D, 2); % center of points
    norm_vec = P0 - Orig; norm_vec = norm_vec / norm(norm_vec)*300;
    x_axis = Ps3D(:,end-22) - Ps3D(:,end-29); x_axis = x_axis / norm(x_axis)*300;
    y_axis = cross(x_axis, norm_vec); y_axis = -y_axis / norm(y_axis)*200;
    x_ps = [Orig Orig + x_axis];    y_ps = [Orig Orig + y_axis];    z_ps = [Orig Orig + norm_vec];    
    if i==1
        figure(1), clf, hold on
        plot3(x_ps(1,:), x_ps(2,:),x_ps(3,:),'g-','LineWidth',2)  % x-axis
        plot3(y_ps(1,:), y_ps(2,:),y_ps(3,:),'b-','LineWidth',2)  % y-axis
        plot3(z_ps(1,:), z_ps(2,:),z_ps(3,:),'k-','LineWidth',2)  % z-axis
        plot3(Orig(1), Orig(2), Orig(3), 'ro','markersize',8,'Linewidth',2)
        plot3(Ps3D(1,:), Ps3D(2,:), Ps3D(3,:), 'r.'), 
        legend('X axis', 'Y axis', 'Z axis', 'Origin', 'Object points')
   end
    % plot3(frame(i).data(1,:), frame(i).data(2,:), frame(i).data(3,:), 'r.')
end
figure(1), % axis square 
grid on, hold off
view([-6 -36])

function D = produce_depth(Ps2D, H, W)
D = zeros(H, W);
N = size(Ps2D,2); % point number
zs = Ps2D(3,:);
vs = round(Ps2D(1,:));  us = round(Ps2D(2,:));
for i=1:N
    D(us(i),vs(i)) = zs(i);
end

function b = change_axis(a)
b = a;              % x = x
b(2,:) = -a(3,:);   % y = -z
b(3,:) = a(2,:);    % z = y;

function rect = find_ROI(BW)
[H, W] = size(BW);
i_min = 100; i_max = 0;
j_min = 100; j_max = 0;
for i=1:H
    for j=1:W
        if BW(i,j)
            if i<i_min, i_min = i; end;
            if i>i_max, i_max = i; end;
            if j<j_min, j_min = j; end;
            if j>j_max, j_max = j; end;
        end
    end
end
rect = [i_min, i_max, j_min, j_max];

function D_out = depth_interp(D_in)
[H, W] = size(D_in);
D_out = zeros(H, W);
% for i=1:H
%     D_i = D_in(i, :);
%     d_test = D_i(D_i>0);
%     fft_d = dct(D_i, W); fft_d(abs(fft_d)<200) = 0;
%     d_rec = idct(fft_d, W);
%     D_out(i,:) = d_rec;
% end

%把cloud_cluster_x_y.txt.mat的列存储格式，变成行存储格式
function [tmp R t] = transformFormat(a)    %取出a的前三列
%Origin
tmp(1,1) = a(1,4);
tmp(2,1) = a(2,4);
tmp(3,1) = a(3,4);
R = a(1:3, 1:3);
t = a(1:3, 4);
[r c] = size(a);%r-3就是所有点的个数
for i = 4:r
    if a(i, 7) == 1
        tmp_c = size(tmp, 2); %取已放进tmp_c的个数
        tmp(1, tmp_c+1) = a(i, 1);
        tmp(2, tmp_c+1) = a(i ,2);
        tmp(3, tmp_c+1) = a(i, 3);
    end
end



%等价于改变每帧观察的视角
function rotated_tmp = rotate(tmp, R, t)
vnum = size(tmp, 2); %取列，含原点
for i = 1:vnum
    v = [tmp(1,i); tmp(2, i); tmp(3, i)];
    v_rotated = R*v ;
   rotated_tmp(1:3, i) = v_rotated;
end


