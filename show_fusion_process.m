function show_fusion_process

% centers = load('data\centers.txt');
ax_seq = [1 2 3];
for i=1:5
    a = load(num2str(i, 'C:\\Users\\hhg\\Desktop\\data_without_trans\\cloud_cluster_%d_0.txt'));
    frame(i).data = change_axis(a);  %#ok<AGROW>
%     frame(i).center = centers(i,ax_seq)';                                    %#ok<AGROW>
end  
% c = linspace(0.4, 0.9, 10);
K = [1000, 0,  150; ...
     0,  100, 15; ...
     0,  0,  1];  % related to the depth
for i = [1]
    Orig = frame(i).data(ax_seq,1)*1000;
    Ps3D = frame(i).data(ax_seq,2:end)*1000;
    %把坐标转换到雷达原点
    Ps3D(1,:) = Ps3D(1, :) - Orig(1);
    Ps3D(2,:) = Ps3D(2, :) - Orig(2);
    Ps3D(3,:) = Ps3D(3, :) - Orig(3);
    
    Ps2D = K*Ps3D; z = Ps3D(3,:);
    Ps2D(1,:) = Ps2D(1,:)./z;
    Ps2D(2,:) = Ps2D(2,:)./z;
    figure(2), plot(Ps2D(1,:), Ps2D(2,:),'k.'),grid on
    P0 = mean(Ps3D, 2); % center of points
    norm_vec = P0 - Orig; norm_vec = norm_vec / norm(norm_vec, 2);
    if i==1
        figure(1), clf
        plot3(Ps3D(1,:), Ps3D(2,:), Ps3D(3,:), 'r.'), hold on
        plot3(Orig(1), Orig(2), Orig(3), 'r*')
        quiver3(Orig(1), Orig(2), Orig(3), norm_vec(1), norm_vec(2), norm_vec(3),'k'); % z-axis
        quiver3(Orig(1), Orig(2), Orig(3), norm_vec(1), norm_vec(2), norm_vec(3),'b'); % z-axis
        quiver3(Orig(1), Orig(2), Orig(3), norm_vec(1), norm_vec(2), norm_vec(3),'c'); % z-axis
    else
        figure(1), clf
        plot3(Ps3D(1,:), Ps3D(2,:), Ps3D(3,:), 'b.'), hold on
        % plot3(Ps3D(1,:), Ps3D(3,:), Ps3D(2,:), '.', 'color', 0.1+[c(i) 1-c(i)  c(i)]), hold on
        plot3(Orig(1), Orig(2), Orig(3), 'b*')
        quiver3(Orig(1), Orig(2), Orig(3), norm_vec(1), norm_vec(2), norm_vec(3),'b','LineWidth',2);
    end
    % plot3(frame(i).data(1,:), frame(i).data(2,:), frame(i).data(3,:), 'r.')
end
figure(1), axis square 
grid on, hold off
view([27 10])

function b = change_axis(a)
b = a;              % x = x
b(2,:) = -a(3,:);   % y = -z
b(3,:) = a(2,:);    % z = y;



