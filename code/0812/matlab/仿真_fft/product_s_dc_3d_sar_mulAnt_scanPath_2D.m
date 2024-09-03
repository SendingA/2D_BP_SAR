function s_dc_all_cell = product_s_dc_3d_sar_mulAnt_scanPath_2D(xy_target, amp_target, sft_xy_mat, ...
    Rx_xyz_0_mat, Tx_xyz_0_mat)

%% 雷达参数
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
% 光速
c = 3 * 10^8;

num_ant_Rx = size(Rx_xyz_0_mat, 1);
num_ant_Tx = size(Tx_xyz_0_mat, 1);

targetNum = size(xy_target, 1);

s_dc_all_cell = cell(1, num_ant_Tx);

% s_dc_all_mat = zeros(num_ant_Rx, num_t_c, num_L);

Rx_xyz_0 = Rx_xyz_0_mat;

num_L = size(sft_xy_mat, 1);

rwin = hamming(num_t_c);
rxNum = num_ant_Rx;
rwinMat = zeros(rxNum, length(rwin));
for i_rxNum = 1:rxNum
    rwinMat(i_rxNum, :) = rwin.';
end

for i_Tx = 1:num_ant_Tx
    
    Tx_xyz_0 = Tx_xyz_0_mat(i_Tx, :);
    s_dc_all_mat = zeros(num_ant_Rx, num_t_c, num_L );
    
    for i_target = 1:targetNum
        Tg_xyz = xy_target(i_target,:);
    
        for i_d = 1:num_L
%         sft_xy_mat
        Rx_xyz = Rx_xyz_0 + sft_xy_mat(i_d, :);
        Rx_Tg = (Rx_xyz(:,1)-Tg_xyz(1,1)).^2 + (Rx_xyz(:,2)-Tg_xyz(1,2)).^2;
        Rx_Tg = Rx_Tg.^0.5;
        
        Tx_xyz = Tx_xyz_0 + sft_xy_mat(i_d, :);
        Tx_Tg = (Tx_xyz(:,1)-Tg_xyz(1,1)).^2 + (Tx_xyz(:,2)-Tg_xyz(1,2)).^2;
        Tx_Tg = Tx_Tg.^0.5;
        
        R_0_all = Rx_Tg + Tx_Tg;
        
%         disp(R_0_all);
        
        s_rx_all = exp(1j*(2*pi*f_c*(t-R_0_all/c) + pi*K*(t-R_0_all/c).^2 ));
        s_tx_all = exp(1j*(2*pi*f_c*t + pi*K*t.^2 ));
        
        s_dc_all = amp_target(i_target, :).*s_tx_all.*conj(s_rx_all).*rwinMat;
        
        s_dc_all_mat(:,:, i_d) = s_dc_all_mat(:,:, i_d) + s_dc_all;
        
        end

        disp([i_target, targetNum]);
        
    end

    s_dc_all_cell{i_Tx} = s_dc_all_mat;
    
end

