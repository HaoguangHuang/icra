%可视化从雷达获得的原始数据


function readOriginData
    %读取10帧原始数据，坐标为雷达世界坐标
    fileNum = 10;
    for i = 1 : fileNum
        a(i) = load(num2str(i, 'cloud_cluster_%d_2.txt.mat'));
             %改原始数据的格式
        [r c] = size(a(i).data_out);
        b(i).data_out(1:r-3, 1:3) = a(i).data_out(4:r, 1:3);
    end
    
         
    %可视化10帧雷达数据
    %figure;
    hold on;
    for i = 1 : fileNum
        [r c] = size(b(i).data_out);  
        plot3(b(i).data_out(: ,1), b(i).data_out(: ,2), b(i).data_out(: , 3), 'r.', 'markersize', 5); 
    end
    xlabel('x', 'FontSize', 30); ylabel('y', 'FontSize', 30); zlabel('z', 'FontSize', 30);
    hold off;
    
    
end



