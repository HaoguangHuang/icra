in = 0;
for i= 1:4
    for j = (i+1):4
        if (i ==1 &&  j==2)
            min = norm(a(i,:) - a(j,:));
        end
        tmp = norm(a(i,:) - a(j,:))
        if (tmp < min)
            min = tmp;
        end
    end
end
min

%% 看kitty所有数据在世界坐标中的位置
figure;
MarkerSpecs = {'+', 'o'};
ColorSpecs = {'r', 'g', 'b', 'k'};
Specs = {'+r', '+g', '+b', '+k', '.r', '.g', '.b', '.k'};
for i = 1:8
    hold on;
    point_size = size(a(i).data_out,2) - 3;
    plot3(a(i).data_out(4:end, 1), a(i).data_out(4:end, 2),a(i).data_out(4:end, 3), Specs{i});
end
plot3(0,0,0,'oy')
hold off;

%%  计算转换到当前帧世界坐标后，每一帧数据到对应原点的最大距离
for i = 1:8
    pointSize = size(frame(i).data,2)-1; %取列长度，剔除原点
    distance = [];
    for j = 2:(pointSize+1)
    distance = [distance norm(frame(i).data(1:3, j))];
    end
    maxD = max(distance);
    minD = min(distance);
    disp(['第',num2str(i),'帧雷达数据到当前帧原点的最大/最小距离为: ']);
    disp(maxD); disp(minD);
end

%% 看calib之后，数据是否被对齐到原点主轴位置了
i = 8;
figure;
plot3(0,0,0, 'r.', 'markersize', 30);
hold on;
pointSize = size(frame(i).data,2) - 1;
for j = 2:(pointSize)   %第一个点是原点，不plot
    plot3(frame(i).data(1,j), frame(i).data(2,j), frame(i).data(3,j), '+b');
    plot3(frame_withoutCalib(i).data(1,j), frame_withoutCalib(i).data(2,j), frame_withoutCalib(i).data(3,j), '+g');
end
xlabel('this is X');
ylabel('this is Y');
zlabel('this is Z');
view([109, -88]);



