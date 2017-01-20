%对箱子数据提前进行ICP
function ICP_xiangzi
    fileNum = 15; %要处理的箱子数据的文件个数
    
    for i = 1:fileNum-1
        address1 = ['E:\Code\ICRA_dir\inputData_fromLidar\lidar_from_ccy\cloud_cluster_',num2str(i),'_2.txt.mat'];
        address2 = ['E:\Code\ICRA_dir\inputData_fromLidar\lidar_from_ccy\cloud_cluster_',num2str(i+1),'_2.txt.mat'];
        a(i).data = load(address1);
        a(i+1).data = load(address2);
        dst = a(i).data.data_out(4:end, 1:3);
        src = a(i+1).data.data_out(4:end, 1:3);
        dst = dst'; src = src';
        [R T] = icp(dst, src, 50); %迭代10次
        display(R);
        display(T);
        src_aligned = R*src + repmat(T,1,length(src));
    %plot看R,T的效果
    
    figure;
    hold on;
    plot3(src(1,:), src(2,:), src(3,:), '.r');
    plot3(dst(1,:), dst(2,:), dst(3,:), '.g');
    plot3(src_aligned(1,:), src_aligned(2,:), src_aligned(3,:), '.b');
    title('kitti的序列3的第7帧与第8帧ICP的结果')
    hold off;   
%% 写入txt文件
    fileAddress = 'E:\Code\ICRA_dir\outputData\ICP_of_xiangzi_TXT\';
        filename = 'xiangzi.txt';
            file_in = fopen([fileAddress, filename], 'at');
    %写入R
    for r = 1:3
        for c = 1:3
            if c == 3
                fprintf(file_in, '%6.4f\n', R(r, c)); %第一个要从新一行开始写
            else
                 fprintf(file_in, '%6.4f\0', R(r, c));
            end
        end
    end
    %写入t
    fprintf(file_in, '%6.4f %6.4f %6.4f\n', T );
    fclose(file_in);

    end
end


