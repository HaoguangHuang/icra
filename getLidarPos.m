%此处的输入frame并没有经过x'=x, y'=-z, z=y的坐标变换。也就是原点在lidar global coo下的坐标
function lidarPos = getLidarPos(frame, i)
    lidarPos(1, 1) = frame(1, 1);
    lidarPos(2, 1) = -frame(3, 1);
    lidarPos(3, 1) = frame(2, 1);
end 