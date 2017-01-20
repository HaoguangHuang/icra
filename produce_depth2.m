%只保留Ps2D中，u,v都大于0的点
function D = produce_depth2(Ps2D, H, W)
    D = zeros(H, W);
    N = size(Ps2D,2); % point number
    zs = Ps2D(3,:);
    us = round(Ps2D(1,:));  vs = round(Ps2D(2,:));%取整
    for i=1:N
        if us(i)>0 && vs(i)>0
        D(vs(i),us(i)) = zs(i);   %先行后列    超出范围的时候，matlab自动补全到新的size
        end
    end
    
    %% 统计D中的非零点
if 1
    non0_sum = 0;
    for r = 1:H
        for c = 1:W
            if(D(r,c)~=0) 
                non0_sum = non0_sum + 1;
            end
        end
    end
    display(non0_sum);
    display(size(Ps2D,2));
end
    
end