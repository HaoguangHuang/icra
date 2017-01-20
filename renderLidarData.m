function renderLidarData

    fx_ = 1000;
    fy_ = 100;
    cx_ = 1000;
    cy_ = 0;
    anglePerJ = 0.2 %��λ��
    pixel_len = 1; %���ص�������,��λmm
    

%��ȡmat��ʽ����
load('cloud_cluster_2_2.txt.mat');

I = data_out;
%Ѱ�Ҳ�ֵ����j
    [r c] = size(I);
    max_j = max(I(2 : r, 5), [], 1);
    min_j = min(I(2 : r, 5), [], 1);

    col = max_j - min_j;
    row = 64;

    %������ԭ��ת������ǰ�״�
   lidar_data_g = I(2:r, :);
   src = I(1, :);
   lidar_data_curr(:,1) = lidar_data_g(:,1) - src(1);
   lidar_data_curr(:,2) = lidar_data_g(:,1) - src(2);
   lidar_data_curr(:,3) = lidar_data_g(:,1) - src(3);
   
    center_j = floor(col/2);
    
    %��ʼrender
    for point_num = 1 : r-1
        v = lidar_data_g(point_num, 4); %v = i;
        Angle = anglePerJ*(lidar_data_g(point_num, 5) - center_j);%��λ��
        Angle_r = pi / 180 * Angle;%���ȵ�λ
        d_mm = cos(Angle_r) * lidar_data_g(point_num, 6);%���������ֵ
        u_mm = tan(Angle_r) * fx_;
        u = u_mm / pixel_len + cx_;
%         image(floor(u), floor(v)) = d_mm;
        u_store(point_num) = u;
        v_store(point_num) = v;
        d_store(point_num) = d_mm;
    end

    %���image
    figure('NumberTitle','off','Name','render������ݵ���ά��ʾ');
%     imshow(image);
plot3(u_store, v_store, d_store, 'r.');
 title('render������ݵ���ά��ʾ'); 
xlabel('x��');
ylabel('y��');

figure('Name','renderǰ�����ݣ��Ѳ�ֵ��')
title('renderǰ�����ݣ��Ѳ�ֵ��');
plot3(I(2:r, 1), I(2:r, 2), I(2:r, 3), 'b*');


figure('Name', 'render������ݵ�(u,v)');
title('render������ݵ�(u,v)');
plot(u_store, v_store, 'r.');


