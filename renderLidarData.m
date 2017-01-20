function renderLidarData

    fx_ = 1000;
    fy_ = 100;
    cx_ = 1000;
    cy_ = 0;
    anglePerJ = 0.2 %单位度
    pixel_len = 1; %像素的物理长度,单位mm
    

%读取mat格式数据
load('cloud_cluster_2_2.txt.mat');

I = data_out;
%寻找差值最大的j
    [r c] = size(I);
    max_j = max(I(2 : r, 5), [], 1);
    min_j = min(I(2 : r, 5), [], 1);

    col = max_j - min_j;
    row = 64;

    %把坐标原点转换到当前雷达
   lidar_data_g = I(2:r, :);
   src = I(1, :);
   lidar_data_curr(:,1) = lidar_data_g(:,1) - src(1);
   lidar_data_curr(:,2) = lidar_data_g(:,1) - src(2);
   lidar_data_curr(:,3) = lidar_data_g(:,1) - src(3);
   
    center_j = floor(col/2);
    
    %开始render
    for point_num = 1 : r-1
        v = lidar_data_g(point_num, 4); %v = i;
        Angle = anglePerJ*(lidar_data_g(point_num, 5) - center_j);%单位度
        Angle_r = pi / 180 * Angle;%弧度单位
        d_mm = cos(Angle_r) * lidar_data_g(point_num, 6);%真正的深度值
        u_mm = tan(Angle_r) * fx_;
        u = u_mm / pixel_len + cx_;
%         image(floor(u), floor(v)) = d_mm;
        u_store(point_num) = u;
        v_store(point_num) = v;
        d_store(point_num) = d_mm;
    end

    %输出image
    figure('NumberTitle','off','Name','render后的数据的三维表示');
%     imshow(image);
plot3(u_store, v_store, d_store, 'r.');
 title('render后的数据的三维表示'); 
xlabel('x轴');
ylabel('y轴');

figure('Name','render前的数据（已插值）')
title('render前的数据（已插值）');
plot3(I(2:r, 1), I(2:r, 2), I(2:r, 3), 'b*');


figure('Name', 'render后的数据的(u,v)');
title('render后的数据的(u,v)');
plot(u_store, v_store, 'r.');


