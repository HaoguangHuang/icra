%�˴�������frame��û�о���x'=x, y'=-z, z=y������任��Ҳ����ԭ����lidar global coo�µ�����
function lidarPos = getLidarPos(frame, i)
    lidarPos(1, 1) = frame(1, 1);
    lidarPos(2, 1) = -frame(3, 1);
    lidarPos(3, 1) = frame(2, 1);
end 