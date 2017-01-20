function show_fusion_process
% figure(1) is to show camera coordinate systems

% centers = load('data\centers.txt');
ax_seq = [1 2 3];
for i=1:5
     a = load(num2str(i, 'cloud_cluster_%d_0.txt'));
    frame(i).data = change_axis(a);  %#ok<AGROW>
%     frame(i).center = centers(i,ax_seq)';                                    %#ok<AGROW>
end
% c = linspace(0.4, 0.9, 10);
K = [1000, 0,  0; ...
     0,  100, 0; ...
     0,  0,  1]; 
H = 64; W = 300;
for i = 1:5
    Orig = frame(i).data(ax_seq,1)*1000;   %ÊÀ½ç×ø±ê
    Ps3D = frame(i).data(ax_seq,2:end)*1000;
    Ps2D = K*Ps3D; z = Ps3D(3,:);
    Ps2D(1,:) = Ps2D(1,:)./z;
    Ps2D(2,:) = Ps2D(2,:)./z;
%     figure(2), plot(Ps2D(1,:), Ps2D(2,:),'k.'),grid on
    
    us = Ps2D(2,:);  vs = Ps2D(1,:); 
    u_min = min(us); u_max = max(us);
    v_min = min(vs); v_max = max(vs);
    Ps2D_new = Ps2D;
    Ps2D_new(1,:) = Ps2D_new(1,:) - v_min + 10;
    Ps2D_new(2,:) = Ps2D_new(2,:) - u_min + 10;
    D_map = produce_depth(Ps2D_new, H, W);
    imwrite(uint16(D_map), num2str(i,'depth_%d.png'));
    figure(11),imshow(mat2gray(D_map)), colormap(jet)
    
    % compute the fusion weights
    BW = uint8(255*logical(D_map));
    ROI = find_ROI(BW);
    psf = fspecial('averag',7);
    Blur_BW = imfilter(BW,psf,'replicate');
    Blur_BW1 = double( double(Blur_BW) / norm(double(Blur_BW(:)),2) ) * 10;
    figure(10), imshow(BW)
    figure(20), imshow(uint8(mat2gray(Blur_BW1)*255),'border','tight')
    mask = zeros(size(BW)); mask(ROI(1):ROI(2), ROI(3):ROI(4)) = 1;
    Blur_BW1(logical(BW)) = 1;  Blur_BW1(~mask) = 0;
    figure(30), imshow(mat2gray(Blur_BW1))
    imwrite(uint8(mat2gray(Blur_BW1)*255), num2str(i,'weight_%d.png'));
    
    % depth interpolation
    D_in = D_map(ROI(1):ROI(2), ROI(3):ROI(4));
    D_out = depth_interp(D_in);
    
    P0 = mean(Ps3D, 2); % center of points
    norm_vec = P0 - Orig; norm_vec = norm_vec / norm(norm_vec)*300;
    x_axis = Ps3D(:,end-22) - Ps3D(:,end-29); x_axis = x_axis / norm(x_axis)*300;
    y_axis = cross(x_axis, norm_vec); y_axis = -y_axis / norm(y_axis)*200;
    x_ps = [Orig Orig + x_axis];    y_ps = [Orig Orig + y_axis];    z_ps = [Orig Orig + norm_vec];    
    if i==1
        figure(1), clf, hold on
        plot3(x_ps(1,:), x_ps(2,:),x_ps(3,:),'g-','LineWidth',2)  % x-axis
        plot3(y_ps(1,:), y_ps(2,:),y_ps(3,:),'b-','LineWidth',2)  % y-axis
        plot3(z_ps(1,:), z_ps(2,:),z_ps(3,:),'k-','LineWidth',2)  % z-axis
        plot3(Orig(1), Orig(2), Orig(3), 'ro','markersize',8,'Linewidth',2)
        plot3(Ps3D(1,:), Ps3D(2,:), Ps3D(3,:), 'r.'), 
        legend('X axis', 'Y axis', 'Z axis', 'Origin', 'Object points')
   end
    % plot3(frame(i).data(1,:), frame(i).data(2,:), frame(i).data(3,:), 'r.')
end
figure(1), % axis square 
grid on, hold off
view([-6 -36])

function D = produce_depth(Ps2D, H, W)
D = zeros(H, W);
N = size(Ps2D,2); % point number
zs = Ps2D(3,:);
vs = round(Ps2D(1,:));  us = round(Ps2D(2,:));
for i=1:N
    D(us(i),vs(i)) = zs(i);
end

function b = change_axis(a)
b = a;              % x = x
b(2,:) = -a(3,:);   % y = -z
b(3,:) = a(2,:);    % z = y;

function rect = find_ROI(BW)
[H, W] = size(BW);
i_min = 100; i_max = 0;
j_min = 100; j_max = 0;
for i=1:H
    for j=1:W
        if BW(i,j)
            if i<i_min, i_min = i; end;
            if i>i_max, i_max = i; end;
            if j<j_min, j_min = j; end;
            if j>j_max, j_max = j; end;
        end
    end
end
rect = [i_min, i_max, j_min, j_max];

function D_out = depth_interp(D_in)
[H, W] = size(D_in);
D_out = zeros(H, W);
% for i=1:H
%     D_i = D_in(i, :);
%     d_test = D_i(D_i>0);
%     fft_d = dct(D_i, W); fft_d(abs(fft_d)<200) = 0;
%     d_rec = idct(fft_d, W);
%     D_out(i,:) = d_rec;
% end

