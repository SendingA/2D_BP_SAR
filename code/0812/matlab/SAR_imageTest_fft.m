loadPath = 'D:\fgf\00可能的新方向\000test\20230515\SAR\实验\data\20230812\minData\';
loadFile = 'adc_data_3target_1';
load([loadPath, loadFile, '.mat']);

x_range = [0 10.8];
y_range = [-3.75 3.75];

% 起始频率 f_c = 77GHz = 77*10^9 Hz
f_c = 77*10^9;
% 频率坡度 K = 60.012 = 60.012*10^12 Hz
K = 60.012*10^12;
% frame间隔时间 t_f = 20ms = 0.02 s
% chirp采用间隔 t_c = 10^(-7) s
t_c = 10^(-7);
num_t_c = 512;
% chirp数量 num_t_c = 512;
t = 0:t_c:t_c*(num_t_c-1);
% t = 2*t_c:t_c:t_c*(num_t_c+1);

c = 3 * 10^8;
Lambda = 0.0039;



keyDot_Sft_xy_Mat = [ ...
    -0.295, 0; ...
    0, 0; ...
    ];

sub_num_L = 31;

lineNum = size(keyDot_Sft_xy_Mat, 1)-1;
% sub_num_L = num_L/lineNum;
num_L = sub_num_L*lineNum;


unitSft_xy_Mat = zeros(lineNum, 2);
for i_keyDot = 1:lineNum
    unitSft_xy_Mat(i_keyDot, :) = (keyDot_Sft_xy_Mat(i_keyDot+1, :) ...
        - keyDot_Sft_xy_Mat(i_keyDot, :) )/(sub_num_L - 1);
end

for i_keyDot = 1:lineNum
    for unitStep = 1:sub_num_L
        sft_xy_mat(unitStep+(i_keyDot-1)*sub_num_L, :) = ...
            keyDot_Sft_xy_Mat(i_keyDot,:) + unitSft_xy_Mat(i_keyDot, :)*(unitStep-1);
    end
end

% 单个雷达
Rx_xyz_0_cell = {[0, 0; ...
        0, 0.00195; ...
        0, 0.0039; ...
        0, 0.00585...
    ]; ...
    };

% Rx_xyz_0_cell = {[0, -0.00585; ...
%         0, -0.0039; ...
%         0, -0.00195; ...
%         0, -0....
%     ]; ...
%     };

Tx_xyz_0_cell = {[0, 0.011; ...
%         0, 0.0149; ...
%         0, 0.0188 ...
    ]; ...
    };

% 天线数
num_ant_Rx = size(Rx_xyz_0_cell{1}, 1);
num_ant_Tx = size(Tx_xyz_0_cell{1}, 1);

% 格子数
num_xAxis = 108;
num_yAxis = 75;


% 检测像素的矩阵坐标
xy_I_mat = cell(num_xAxis, num_yAxis);
for i_xAxis = 1:num_xAxis
    for i_yAxis = 1:num_yAxis
        xy_I_mat{i_xAxis, i_yAxis} = [ ...
            x_range(1)+(i_xAxis-1)*(x_range(2)-x_range(1))/(num_xAxis-1), ...
            y_range(1)+(i_yAxis-1)*(y_range(2)-y_range(1))/(num_yAxis-1), ...
            ];
    end
end

radarNum = size(Rx_xyz_0_cell, 1);

s_dc_allRadar_cell = cell(radarNum, 1);

for i_radar = 1:radarNum
    
    Tx_xyz_0_mat = Tx_xyz_0_cell{i_radar};
    
    Rx_xyz_0_mat = Rx_xyz_0_cell{i_radar};
    
    %     s_dc_all_cell = product_s_dc_3d_sar_mulAnt_scanPath(xy_target, amp_target, sft_xy_mat, ...
    %         Rx_xyz_0_mat, Tx_xyz_0_mat);
    s_dc_all_cell = cell(1, num_ant_Tx);
    
    s_dc_all_cell{1,1} = sarData;
    
    s_dc_allRadar_cell{i_radar} = s_dc_all_cell;
    
end



I_mat = zeros(num_xAxis, num_yAxis);

tic

for i_radar = 1:radarNum
    
    Tx_xyz_0_mat = Tx_xyz_0_cell{i_radar};
    Rx_xyz_0_mat = Rx_xyz_0_cell{i_radar};
    
    Rx_xyz_0 = Rx_xyz_0_mat;
    
    for i_xAxis = 1:num_xAxis
        for i_yAxis = 1:num_yAxis
            
            Tg_xyz = xy_I_mat{i_xAxis, i_yAxis}; %假设目标的坐标
            
            for i_d = 1:num_L
                
                for i_Tx = 1:num_ant_Tx
                    
                    Tx_xyz_0 = Tx_xyz_0_mat(i_Tx, :);
                    
                    Rx_xyz = Rx_xyz_0 + sft_xy_mat(i_d, :);
                    Rx_Tg = (Rx_xyz(:,1)-Tg_xyz(1,1)).^2 + (Rx_xyz(:,2)-Tg_xyz(1,2)).^2;
                    Rx_Tg = Rx_Tg.^0.5;
                    
                    Tx_xyz = Tx_xyz_0 + sft_xy_mat(i_d, :);
                    Tx_Tg = (Tx_xyz(:,1)-Tg_xyz(1,1)).^2 + (Tx_xyz(:,2)-Tg_xyz(1,2)).^2;
                    Tx_Tg = Tx_Tg.^0.5;
                    
                    R_0_all = Rx_Tg + Tx_Tg;
                    f_r_all = R_0_all*K/c;
                    
                    s_dc_all_mat = s_dc_allRadar_cell{i_radar}{i_Tx};
                    
                    s_rc_all_fft = fft(s_dc_all_mat(:,:,i_d), 512*50, 2 );
                    %                     K = 60.012*10^12;
                    % % frame间隔时间 t_f = 20ms = 0.02 s
                    % % chirp采用间隔 t_c = 10^(-7) s
                    % t_c = 10^(-7);
                    % num_t_c = 512;

                    binR = c/2/(num_t_c*t_c*K)/50;
                    rangeBin_distanceMat = binR*(1:512*50);
                    for i_rx = 1:4
                        [~, indexBin] = min( abs(rangeBin_distanceMat- R_0_all(i_rx) ));
                        fft_indexMat(i_rx) = indexBin;
                    end
                    
                    for i_rx = 1:4
                        s_rc_all(i_rx,1) = s_rc_all_fft(i_rx, fft_indexMat(i_rx) );
                    end
%                     s_rc_all = sum(s_dc_all_mat(:,:,i_d) ...
%                         .*exp(-1j*2*pi*f_r_all*t), 2);
                    
                    cur_I = sum( s_rc_all .* exp(-1j*2*pi*R_0_all/Lambda) );

                    I_mat(i_xAxis, i_yAxis) = I_mat(i_xAxis, i_yAxis) + (cur_I);

                end
            end
            
            disp([i_radar]);
            disp([i_xAxis, i_yAxis]);
            toc
        end
    end
    
end

% I_mat = log(abs(I_mat));

I_mat_display = abs(I_mat);

% I_mat_display = abs(I_mat(5:70,:));

x_axiz = x_range(1):(x_range(2)-x_range(1))/(num_xAxis-1):x_range(2);
y_axiz = y_range(1):(y_range(2)-y_range(1))/(num_yAxis-1):y_range(2);

figure;
set(gcf,'Position',[250,250,650,400]);
pcolor(x_axiz,y_axiz, I_mat_display');
axis equal
xlabel('X轴');
ylabel('Y轴');
shading flat
xlim(x_range);
ylim(y_range);



