function png2pcd_afterInterpolate(rfx, rfy, W, H, fileNum, series)
%fileNum = 10;
%�����png2pcd���޸Ĺ���
%format 1������Ϊbinary��pcd�ļ���0������Ϊasii��pcd�ļ������ڵ�ICRA kinfu����ֻ�ܶ�ȡbinary��ʽ��pcd�ļ�
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