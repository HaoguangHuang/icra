%首先把总的矩阵整理好，按照每帧lidarRotation_inv, lidarTranslation, calibRotation_inv，共7*fileNum行
%{   
格式：  一共7行
[  R
    t
  calib ]
%}
       
function outputTransformation(lidarRotation_inv, lidarTranslation, calibRotation_inv, fileNum, rfx, rfy, W, H, series)  
    totalTransformation = [];
    for i =1:fileNum
        totalTransformation = [totalTransformation; lidarRotation_inv(3*i-2:3*i, 1:3); lidarTranslation(i, 1:3); calibRotation_inv(3*i-2:3*i, 1:3)];
    end
    
    %把矩阵输出到txt
%     rfx = '_1200_';
%     rfy = '400_H_';
%     y_angle = 'jmin_m10';%m是minors
%     x_angle = 'm2';
%     filename = ['Transformation_', 'frame_', num2str(fileNum), rfx, rfy, num2str(H)];
    fileAddress = 'E:\Code\ICRA_dir\outputData\Transformation_with_lidarRotation_lidarTranslatoin_calibRotation_TXT\';
    %filename = 'Transformation_withCalib_rfx_1200_rfy_400_W_320_H_96_fileNum_10_series_2.txt';
    filename = ['Transformation_withCalib_rfx_', num2str(rfx), '_rfy_', num2str(rfy), '_W_', num2str(W), '_H_', num2str(H), ...
        '_fileNum_', num2str(fileNum), '_series_', num2str(series), '.txt'];
    
    %dlmwrite([fileAddress, filename], totalTransformation);
    file_in = fopen([fileAddress, filename], 'wt');
    [r c] = size(totalTransformation);
    for i = 1:r
        for j = 1:c
            if j == c
                fprintf(file_in, '%g\n', totalTransformation(i, j));
            else
                fprintf(file_in, '%g\0', totalTransformation(i, j));
            end
        end
    end
    fclose(file_in);
end



