function show_fusion_process_with_weightMap_for64(rfx, rfy, W, H, fileNum, series)
% figure(1) is to show camera coordinate systems

ax_seq = [1 2 3];
%fileNum = 10;
lidar_j_DeviationAngle = 0.2 * pi /180;%����
kitti_lidar_j_DeviationAngle = 0.08 * pi /180;

lidar_i_DeviationAngle = 2*pi/180;  %��ֱ����
kitti_lidar_i_DeviationAngle = 0.4 * pi/180;

lidarRotation_inv = []; %��¼ǰi֡���global coo�ı任
lidarTranslation = [];
calibRotation_inv = [];
%  lidarSrc_address = 'E:\Code\ICRA_dir\inputData_fromLidar\dataFromKitti\3\';
 %lidarSrc_address = 'E:\Code\ICRA_dir\inputData_fromLidar\lidar_from_ccy\';
lidarSrc_address = 'E:\Code\ICRA_dir\inputData_fromLidar\kitti_motor_5to15\';
 
%��ȡ10֡.mat����
for i=1:fileNum
     a(i) = load([lidarSrc_address, 'cloud_cluster_',num2str(series),'_', num2str(i), '.txt.mat']);
%     a(i) = load([lidarSrc_address, 'cloud_cluster_',num2str(i),'_', num2str(series), '.txt.mat']);
    
    %R���״��������굽��ǰ֡�������ת
    %t���״��������굽��ǰ֡�����ƽ��
    %tmp��4*n������������ԭ�㡣ǰ������xyz�����������״�Ƕ�j
    [tmp R t] = transformFormat(a(i).data_out); 
   
    %����״ﵱǰλ�����������е�λ��
    lidarPos_withChangeAxis(1:3, i) = getLidarPos(tmp, i);
    
    %tmp��4*n������������ԭ�㡣ǰ������xyz�����������״�Ƕ�j
    %tmp_rotated��frame(i).data����ԭ��������ݣ�ȫ����ת����Ч����������ת��
    tmp_rotated = rotate(tmp, R , t);
    frame_withoutCalib(i).data = change_axis(tmp_rotated); %frame��Ԫ��Ϊ.mat������
   
    %�����ǰframe��Ӧ��������global coo�����ת�ǹ�ϵ
    j_min = min(frame_withoutCalib(i).data(4, 2:end));
    j_max = max(frame_withoutCalib(i).data(4, 2:end));
%     zRotationAngle = (j_min+(j_min - j_max)/2)*lidar_j_DeviationAngle/2; 

%     yRotationAngle =-( j_min)*lidar_j_DeviationAngle;%Խ��Խ����
%     xRotationAngle = 10 * lidar_i_DeviationAngle;
    yRotationAngle =-j_min*kitti_lidar_j_DeviationAngle;%Խ��Խ����
    xRotationAngle = -10*kitti_lidar_i_DeviationAngle; %Խ��Խ����
   
    Rotation = rodrigues([xRotationAngle yRotationAngle 0]);%��y����ת
    frame(i).data = calib(frame_withoutCalib(i).data, Rotation);
    
    %ע���ʱ��R,t��change_axis֮ǰ��
    lidarRotation_inv = [lidarRotation_inv; inv(R)];  %�����ṩ��R,t������������һ���������ʱ��R,tû��change_axis
    lidarTranslation = [lidarTranslation; t'];  %t�Ǻ��Ŵ��
    calibRotation_inv = [calibRotation_inv; inv(Rotation)];%����R,t�����ע���ʱ��change_axis֮��
    
    
    %������ϵ�Rotation��translation
    %lidarRotation_withChangeAxis = getTotalRotation(R, Rotation);
    %lidarTranslation_withChangeAxis(1:3, 1:3) = lidarTranslation;
end


%���lidarRotation_inv, lidarTranslation, calibRotation_inv��һ��txt�ļ�
outputTransformation(lidarRotation_inv, lidarTranslation, calibRotation_inv, fileNum, rfx, rfy, W, H,  series);


% K = [500, 0,  0; ...
%      0,  200, 0; ...
%      0,  0,  1]; 
K = [rfx, 0,  0; ...
     0,  rfy, 0; ...
     0,  0,  1]; 
cx = W/2; cy = H/2;
 %H = 96; W = 320;
for i = 1 :fileNum
    Orig = frame(i).data(ax_seq,1)*1000;   %��������
    Ps3D = frame(i).data(ax_seq,2:end)*1000;
    Ps2D = K*Ps3D; 
    z = Ps3D(3,:);
    Ps2D(1,:) = Ps2D(1,:)./z + cx;
    Ps2D(2,:) = Ps2D(2,:)./z + cy;
%      Ps2D(1,:) = Ps2D(1,:)./z;
%     Ps2D(2,:) = Ps2D(2,:)./z;
    figure, plot(Ps2D(1,:), Ps2D(2,:),'k.'),grid on
 
   
     us = Ps2D(2,:);  vs = Ps2D(1,:); 
     u_min = min(us); u_max = max(us);
     v_min = min(vs); v_max = max(vs);
    Ps2D_new = Ps2D;  %���Ps2D_new�������ת���ǲ������
    Ps2D_new(1,:) = Ps2D_new(1,:) - v_min + 10;
    Ps2D_new(2,:) = Ps2D_new(2,:) - u_min + 10;
if 0
    D_map = produce_depth(Ps2D_new, H, W);
else
    D_map = produce_depth2(Ps2D, H, W);  %Ps2D����ԭ��
end
    

   depthMap_output_fileAddress = 'E:\Code\ICRA_dir\outputData\depthMap_from_lidarSrc_PNG\';
    %depthMap_output_fileName =num2str(i,'depth_afterInterpolate_withCalib_rfx_500_rfy_200_W_96_H_320_%d_2.png');
   depthMap_output_fileName = ['depth_afterInterpolate_withCalib_rfx_', num2str(rfx), ...
       '_rfy_', num2str(rfy), '_W_', num2str(W), '_H_', num2str(H), '_fileNum_', num2str(fileNum), '_', num2str(i), '_series_', ...
       num2str(series), '.png'];
   imwrite(uint16(D_map), [depthMap_output_fileAddress, depthMap_output_fileName]);
    figure ,imshow(mat2gray(D_map)), %colormap(jet);
    colormap(summer);
%end
   
  
    % compute the fusion weights
    BW = uint8(255*logical(D_map));
    ROI = find_ROI(BW);
    psf = fspecial('averag',7);
    Blur_BW = imfilter(BW,psf,'replicate');
    Blur_BW1 = double( double(Blur_BW) / norm(double(Blur_BW(:)),2) ) * 10;
    %figure(110), imshow(BW), title('������ù�');
    %figure(120), imshow(uint8(mat2gray(Blur_BW1)*255),'border','tight')
    mask = zeros(size(BW)); mask(ROI(1):ROI(2), ROI(3):ROI(4)) = 1;
    Blur_BW1(logical(BW)) = 1;  Blur_BW1(~mask) = 0;
    %figure(130), imshow(mat2gray(Blur_BW1))
    %imwrite(uint8(mat2gray(Blur_BW1)*255), num2str(i,'weight_afterInterpolate_%d_2.png'));
    
    % depth interpolation
    D_in = D_map(ROI(1):ROI(2), ROI(3):ROI(4));
    D_out = depth_interp(D_in);
    
    P0 = mean(Ps3D, 2); % center of points
    norm_vec = P0 - Orig; norm_vec = norm_vec / norm(norm_vec)*300;
    y_axis = Ps3D(:,end-29) - Ps3D(:,end-22); y_axis = y_axis / norm(y_axis)*250;
    
    
    x_axis = cross(y_axis, norm_vec); x_axis = x_axis / norm(x_axis)*600;
    x_ps = [Orig Orig + x_axis];    y_ps = [Orig Orig + y_axis];    z_ps = [Orig Orig + norm_vec];    
    if i==-1
        figure(111), clf, hold on
        plot3(x_ps(1,:), x_ps(2,:),x_ps(3,:),'g--','LineWidth',2)  % x-axis
        plot3(y_ps(1,:), y_ps(2,:),y_ps(3,:),'b--','LineWidth',2)  % y-axis
        plot3(z_ps(1,:), z_ps(2,:),z_ps(3,:),'k--','LineWidth',2)  % z-axis
        plot3(Orig(1), Orig(2), Orig(3), 'ro','markersize',8,'Linewidth',2)
        plot3(Ps3D(1,:), Ps3D(2,:), Ps3D(3,:), 'r.'), 
        legend('X axis', 'Y axis', 'Z axis', 'Origin', 'Object points')
   end
     %plot3(frame(i).data(1,:), frame(i).data(2,:), frame(i).data(3,:), 'r.')
end
%% figure(121), % axis square 
%grid on, hold off
%view([-18 -40])



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
b(4,:) = a(4,:)


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

%��cloud_cluster_x_y.txt.mat���д洢��ʽ������д洢��ʽ
%    tmp = [
%       orig(1), data_1(1), data_2(1), ... ;
%       orig(2), data_1(2), data_2(2), ... ;
%       orig(3), data_1(3), data_2(3), ... ;
%       j_deviation,  j_deviation_1, j_deviation_2, ... ;
%               ]
function [tmp R t] = transformFormat(a)    %ȡ��a��ǰ����
%ԭ�㣬���ԭ�㲻��(0 0 0)
tmp(1,1) = a(1,4);%��ǰ֡���״���������ת���״ﵱǰ֡�������תR
tmp(2,1) = a(2,4);
tmp(3,1) = a(3,4);
tmp(4,1) = 0;%�״�ԭ���ƫ��j=0
[r c] = size(a);
tmp(1,2:r-2) = a(4:end,1);
tmp(2,2:r-2) = a(4:end,2);
tmp(3,2:r-2) = a(4:end,3);
tmp(4,2:r-2) = a(4:end,5);
R = a(1:3, 1:3);
t = a(1:3, 4);


%��ÿ֡����������ת������ǰ֡����
function rotated_tmp = rotate(tmp, R, t)
vnum = size(tmp, 2); %ȡ�У���ԭ��
for i = 1:vnum
    v_3d = [tmp(1,i); tmp(2, i); tmp(3, i)];
    if i == 1
        rotated_tmp(1:3, i) = [0; 0; 0]; %origin
    else
        v_rotated = inv(R)*(v_3d - t) ;
        rotated_tmp(1:3, i) = v_rotated;
    end
    rotated_tmp(4, i) = tmp(4, i); %���� j
end

%���״������׼�㼯
function frame = calib(frame_withoutCalib, Rotation)
v_num = size(frame_withoutCalib, 2);
Rotation_inv = inv(Rotation);
frame(1, 1) = 0;
frame(2, 1) = 0;
frame(3, 1) = 0;
frame(4, 1) = 0;
for i = 2:v_num
    v1 = [frame_withoutCalib(1,i); frame_withoutCalib(2, i); frame_withoutCalib(3, i)];
%     v2 = Rotation_inv * v1;
    v2 = Rotation * v1;
    frame(1:3, i) = v2;
    frame(4, i) = frame_withoutCalib(4, i);
end

