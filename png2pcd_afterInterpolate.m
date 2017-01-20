function png2pcd_afterInterpolate(rfx, rfy, W, H, fileNum, series)
%fileNum = 10;
%这里的png2pcd是修改过的
%format 1代表保存为binary的pcd文件，0代表保存为asii的pcd文件。现在的ICRA kinfu代码只能读取binary格式的pcd文件
filename_prefix = ['E:\Code\ICRA_dir\code\vs\png2pcd\build\Debug\png2pcd.exe ', ...
'E:\Code\ICRA_dir\outputData\depthMap_from_lidarSrc_PNG\depth_afterInterpolate_withCalib_rfx_' , num2str(rfx), ...
'_rfy_', num2str(rfy), '_W_', num2str(W), '_H_', num2str(H), '_fileNum_', num2str(fileNum), '_'];
filename_midfix = ['_series_', num2str(series), '.png E:\Code\ICRA_dir\outputData\png_to_pcd_PCD\depth_afterInterpolate_withCalib_rfx_', ...
    num2str(rfx), '_rfy_', num2str(rfy), '_W_', num2str(W), '_H_', num2str(H), '_fileNum_', num2str(fileNum), '_'];
filename_surfix = ['_series_', num2str(series), '.pcd --intensity_type FLOAT -format 1'];
    for i = 1 : fileNum
        dos([filename_prefix, num2str(i), filename_midfix, num2str(i), filename_surfix]);
    end
end