%ֻ����Ps2D�У�u,v������0�ĵ�
function D = produce_depth2(Ps2D, H, W)
    D = zeros(H, W);
    N = size(Ps2D,2); % point number
    zs = Ps2D(3,:);
    us = round(Ps2D(1,:));  vs = round(Ps2D(2,:));%ȡ��
    for i=1:N
        if us(i)>0 && vs(i)>0
        D(vs(i),us(i)) = zs(i);   %���к���    ������Χ��ʱ��matlab�Զ���ȫ���µ�size
        end
    end
    
    %% ͳ��D�еķ����
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