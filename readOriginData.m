%���ӻ����״��õ�ԭʼ����


function readOriginData
    %��ȡ10֡ԭʼ���ݣ�����Ϊ�״���������
    fileNum = 10;
    for i = 1 : fileNum
        a(i) = load(num2str(i, 'cloud_cluster_%d_2.txt.mat'));
             %��ԭʼ���ݵĸ�ʽ
        [r c] = size(a(i).data_out);
        b(i).data_out(1:r-3, 1:3) = a(i).data_out(4:r, 1:3);
    end
    
         
    %���ӻ�10֡�״�����
    %figure;
    hold on;
    for i = 1 : fileNum
        [r c] = size(b(i).data_out);  
        plot3(b(i).data_out(: ,1), b(i).data_out(: ,2), b(i).data_out(: , 3), 'r.', 'markersize', 5); 
    end
    xlabel('x', 'FontSize', 30); ylabel('y', 'FontSize', 30); zlabel('z', 'FontSize', 30);
    hold off;
    
    
end



