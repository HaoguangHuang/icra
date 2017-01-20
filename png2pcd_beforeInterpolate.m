fileNum = 10;
filename_prefix = 'E:\Code\C++\pcl\png2pcd\build\Debug\png2pcd.exe C:\Users\hhg\Desktop\ICRA\depth_beforeInterpolate_withoutRotate_';
filename_midfix = '_2.png C:\Users\hhg\Desktop\ICRA\png2pcd\depth_beforeInterpolate_withoutRotate_';
filename_surfix = '_2.pcd --intensity_type FLOAT -format 1';
for i = 1 : fileNum
    dos([filename_prefix, num2str(i), filename_midfix, num2str(i), filename_surfix]);
end