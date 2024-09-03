% loadName = 'adc_data_ganLow';
loadName = 'adc_data_rectLow';
% loadName = 'adc_data_emptyLow';
load(['D:\fgf\00可能的新方向\000test\20230515\SAR\实验\data\20230729\minData\', ...
    loadName, '.mat'])

% 起始频率 f_c = 77GHz = 77*10^9 Hz
f_c = 77*10^9;
% 频率坡度 K = 60.012 = 60.012*10^12 Hz
K = 60.012*10^12;
num_t_c = 512;
t_c = 10^(-7);
t = 0:t_c:t_c*(num_t_c-1);


% numL_mat = [30];
numL_mat = [31 31];

num_L = sum(numL_mat);

start_Sft_xy = [0, -0.15, 0];

shift_xy = [0, 0.01, 0.0; ...
    0, 0.0, 0.01];

% shift_xy = [ ...
%     0, 0.01, 0.0; ...
%     0, 0.0, 0.01];

% shift_xy = [0, 0.01, 0.00;];
% shift_xy = [0, 0.0, 0.01;];

% shift_xy = [0, shift_unit, 0.0;];

sft_xy_mat = zeros(num_L, 3);
cur_L = 0;

cur_L = cur_L + 1;
sft_xy_mat(1, :) = start_Sft_xy;
for unitStep = 1:numL_mat(1)-1
    cur_L = cur_L + 1;
    sft_xy_mat(cur_L, :) = ...
        sft_xy_mat(cur_L-1, :) + shift_xy(1,:);
    
end

sft_xy_mat(cur_L+1, :) = sft_xy_mat(cur_L, :);
cur_L = cur_L + 1;
for unitStep = 1:numL_mat(2)-1
    cur_L = cur_L + 1;
    sft_xy_mat(cur_L, :) = sft_xy_mat(cur_L-1, :) + shift_xy(2,:);
end

sft_xy_mat(1:31,:) = sft_xy_mat(32:end,:);

% x_range = [0.2 1];
% y_range = [-0.4 0.4];
% z_range = [-0.4 0.4];

x_range = [0.2 1];
y_range = [-0.4 0.4];
z_range = [-0.4 0.4];

% 光速
c = 3 * 10^8;
Lambda = 0.0039;

% 接收天线, 发射天线
% 单个雷达
Rx_xyz_0_cell = {[0, 0, 0; ...
    0, 0.00195, 0; ...
    0, 0.0039, 0; ...
    0, 0.00585, 0 ...
    ]; ...
    };

Tx_xyz_0_cell = {[0, 0, 0; ...
%     0, 0.00195, 0.0149; ...
%     0, 0, 0.0188 ...
    ]; ...
    };

% 天线数
num_ant_Rx = size(Rx_xyz_0_cell{1}, 1);
num_ant_Tx = size(Tx_xyz_0_cell{1}, 1);

num_xAxis = 21;
num_yAxis = 21;
num_zAxis = 21;

% 检测像素的矩阵坐标
xy_I_mat = cell(num_xAxis, num_yAxis, num_zAxis);
for i_xAxis = 1:num_xAxis
    for i_yAxis = 1:num_yAxis
        for i_zAxis = 1:num_zAxis
            xy_I_mat{i_xAxis, i_yAxis, i_zAxis} = [x_range(1)+(i_xAxis-1)*(x_range(2)-x_range(1))/(num_xAxis-1), ...
                y_range(1)+(i_yAxis-1)*(y_range(2)-y_range(1))/(num_yAxis-1), ...
                z_range(1)+(i_zAxis-1)*(z_range(2)-z_range(1))/(num_zAxis-1)];
            %         xy_I_mat{i_xAxis, i_yAxis} = [0+(i_xAxis-1)*8/(num_xAxis-1), -4+(i_yAxis-1)*8/(num_yAxis-1)];
        end
    end
end

radarNum = size(Rx_xyz_0_cell, 1);

s_dc_allRadar_cell = cell(radarNum, 1);

for i_radar = 1:radarNum
    
    Tx_xyz_0_mat = Tx_xyz_0_cell{i_radar};
    
    Rx_xyz_0_mat = Rx_xyz_0_cell{i_radar};
    
%     s_dc_all_cell = product_s_dc_3d_sar_mulAnt_scanPath(xy_target, amp_target, sft_xy_mat, ...
%         Rx_xyz_0_mat, Tx_xyz_0_mat);
%     
%     s_dc_allRadar_cell{i_radar} = s_dc_all_cell;

    s_dc_allRadar_cell{i_radar}{1} = sarData(:,:,:);

end









I_mat = zeros(num_xAxis, num_yAxis, num_zAxis);

tic

for i_radar = 1:radarNum
    
    Tx_xyz_0_mat = Tx_xyz_0_cell{i_radar};
    Rx_xyz_0_mat = Rx_xyz_0_cell{i_radar};
%     s_dc_all_cell = product_s_dc_3d_sar_mulAnt(xy_target, amp_target, staSft_xy, num_L, unitSft_xy, ...
%         Rx_xyz_0_mat, Tx_xyz_0_mat);
    
    Rx_xyz_0 = Rx_xyz_0_mat;
    
    for i_xAxis = 1:num_xAxis
        for i_yAxis = 1:num_yAxis
            for i_zAxis = 1:num_zAxis
                
                Tg_xyz = xy_I_mat{i_xAxis, i_yAxis, i_zAxis}; %假设目标的坐标
                
                for i_d = 1:num_L
%                 for i_d = 1:31
    
                    
                    for i_Tx = 1:num_ant_Tx
                        
                        Tx_xyz_0 = Tx_xyz_0_mat(i_Tx, :);
                        
                        Rx_xyz = Rx_xyz_0 + sft_xy_mat(i_d, :);
                        Rx_Tg = (Rx_xyz(:,1)-Tg_xyz(1,1)).^2 + (Rx_xyz(:,2)-Tg_xyz(1,2)).^2 + (Rx_xyz(:,3)-Tg_xyz(1,3)).^2;
                        Rx_Tg = Rx_Tg.^0.5;
                        
                        Tx_xyz = Tx_xyz_0 + sft_xy_mat(i_d, :);
                        Tx_Tg = (Tx_xyz(:,1)-Tg_xyz(1,1)).^2 + (Tx_xyz(:,2)-Tg_xyz(1,2)).^2+ (Tx_xyz(:,3)-Tg_xyz(1,3)).^2;
                        Tx_Tg = Tx_Tg.^0.5;
                        
                        R_0_all = Rx_Tg + Tx_Tg;
                        f_r_all = R_0_all*K/c;
                        
                        s_dc_all_mat = s_dc_allRadar_cell{i_radar}{i_Tx};
                        
                        s_rc_all = sum(s_dc_all_mat(:,:,i_d)...
                            .*exp(-1j*2*pi*f_r_all*t), 2);
                        
                        cur_I = mean(s_rc_all .* exp(-1j*2*pi*R_0_all/Lambda));
                        
%                         I_mat(i_xAxis, i_yAxis, i_zAxis) = I_mat(i_xAxis, i_yAxis, i_zAxis) + abs(cur_I);
                        I_mat(i_xAxis, i_yAxis, i_zAxis) = I_mat(i_xAxis, i_yAxis, i_zAxis) + (cur_I);
                        
                    end
                    
                end
            end
            
            disp([i_radar]);
            disp([i_xAxis, i_yAxis, i_zAxis]);
            toc
        end
    end
    
end

I_mat = abs(I_mat);

xyI_mat = max(I_mat,[], 3);

xzI_mat_temp = max(I_mat,[],2);
xzI_mat = zeros(num_xAxis, num_zAxis);
for i_xAxis = 1:num_xAxis
    for i_zAxis = 1:num_zAxis
        xzI_mat(i_xAxis, i_zAxis) = xzI_mat_temp(i_xAxis,1,i_zAxis);
    end
end

yzI_mat_temp = max(I_mat,[], 1);
yzI_mat = zeros(num_yAxis, num_zAxis);
for i_yAxis = 1:num_yAxis
    for i_zAxis = 1:num_zAxis
        yzI_mat(i_yAxis, i_zAxis) = yzI_mat_temp(1,i_yAxis,i_zAxis);
    end
end

x_axiz = x_range(1):(x_range(2)-x_range(1))/(num_xAxis-1):x_range(2);
y_axiz = y_range(1):(y_range(2)-y_range(1))/(num_yAxis-1):y_range(2);
z_axiz = z_range(1):(z_range(2)-z_range(1))/(num_zAxis-1):z_range(2);

figure;
set(gcf,'Position',[300,300,900,250]);

subplot(1,3,1);
pcolor(x_axiz,y_axiz, xyI_mat');
xlabel('X轴');
ylabel('Y轴');
shading flat
axis equal;
xlim(x_range);
ylim(y_range);

subplot(1,3,2);
pcolor(x_axiz,z_axiz, xzI_mat');
xlabel('X轴');
ylabel('Z轴');
shading flat
axis equal;
xlim(x_range);
ylim(z_range);

subplot(1,3,3);
pcolor(y_axiz,z_axiz, yzI_mat');
xlabel('Y轴');
ylabel('Z轴');
shading flat
axis equal;
xlim(y_range);
ylim(z_range);

I_mat_dot = I_mat;
I_mat_dot(I_mat<0.7*max(I_mat(:)) ) = 0;

[x_dot_ind, y_dot_ind, z_dot_ind] = ind2sub(size(I_mat_dot), find(I_mat_dot));

I_v = zeros(length(x_dot_ind), 1);
for i_I_v = 1:length(x_dot_ind)
    I_v(i_I_v) = I_mat_dot(x_dot_ind(i_I_v), y_dot_ind(i_I_v), z_dot_ind(i_I_v) );
end

figure;
scatter3(x_axiz(x_dot_ind), y_axiz(y_dot_ind), z_axiz(z_dot_ind), [], I_v, 'filled');
xlim(x_range);
ylim(y_range);
zlim(z_range);
colorbar
xlabel('X轴');
ylabel('Y轴');
zlabel('Z轴');










