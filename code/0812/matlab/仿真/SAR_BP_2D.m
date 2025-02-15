
invl = 0.05;
x_range = [0 10.8];
y_range = [-3.75 3.75];

% (3, -2) (3, -1.5) (3.6, -1.5)
target1_temp1 = [ones(length(-2:invl:-1.5), 1)*3 (-2:invl:-1.5)'];
target1_temp2 = [(3:invl:3.6)' ones(length(3:invl:3.6), 1)*(-1.5)];
target1 = [target1_temp1(1:end-1, :); target1_temp2];

target2_temp1 = [ones(length(1.5:invl:2), 1)*3.6 (1.5:invl:2)'];
target2_temp2 = [(3.6:invl:4.6)' ones(length(3.6:invl:4.6), 1)*1.5];
target2 = [target2_temp1(2:end, :); target2_temp2];

target3 = [ones(length(-0.5:invl:0.5), 1)*5.4 (-0.5:invl:0.5)'];


target4_temp1 = [ones(length(3:invl:3.75), 1)*6.3 (3:invl:3.75)'];
target4_temp2 = [(6.3:invl:9)' ones(length(6.3:invl:9), 1)*3];
target4 = [target4_temp1(2:end, :); target4_temp2];

target5 = [ones(length(-3:invl:3), 1)*9 (-3:invl:3)'];
target5(1, :) = [];
target5(end, :) = [];

target6_temp1 = [ones(length(-3.75:invl:-1), 1)*6.3 (-3.75:invl:-1)'];
target6_temp2 = [(6.3:invl:9)' ones(length(6.3:invl:9), 1)*(-3)];
target6 = [target6_temp1(1:end-1, :); target6_temp2];

xy_target = [5 ,2];

amp_target = ones(size(xy_target, 1), 1);

figure;
subplot(2,2,1);
subtitle('模拟点位置')
set(gcf,'Position',[250,250,650,400]);
scatter(xy_target(:,1), xy_target(:,2), 'filled');
axis equal
xlim(x_range);
ylim(y_range);








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

c = 3 * 10^8;
Lambda = 0.0039;

keyDot_Sft_xy_Mat = [ ...
    0.001, 0; ...
    0, 0; ...
    ];

sub_num_L = 2;

lineNum = size(keyDot_Sft_xy_Mat, 1)-1;
% sub_num_L = num_L/lineNum;
num_L = sub_num_L*lineNum;
sft_xy_mat = zeros(num_L, 2);

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
Rx_xyz_0_cell_prod = {[0, 0; ...
        0, 0.00195; ...
        0, 0.0039; ...
        0, 0.00585...
    ]; ...
    };

Tx_xyz_0_cell_prod = {[0, 0.011...
    ]; ...
    };

% 单个雷达
Rx_xyz_0_cell = {[0, 0; ...
        0, 0.00195; ...
        0, 0.0039; ...
        0, 0.00585...
    ]; ...
    };

Tx_xyz_0_cell = {[0, 0.011 ...
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
    
    Tx_xyz_0_mat = Tx_xyz_0_cell_prod{i_radar};
    
    Rx_xyz_0_mat = Rx_xyz_0_cell_prod{i_radar};
    
    %     s_dc_all_cell = product_s_dc_3d_sar_mulAnt_scanPath(xy_target, amp_target, sft_xy_mat, ...
    %         Rx_xyz_0_mat, Tx_xyz_0_mat);
    
    s_dc_all_cell = product_s_dc_3d_sar_mulAnt_scanPath_2D(xy_target, amp_target, sft_xy_mat, ...
        Rx_xyz_0_mat, Tx_xyz_0_mat );
    
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
                    
                    s_rc_all = sum(s_dc_all_mat(:,:,i_d) ...
                        .*exp(-1j*2*pi*f_r_all*t), 2);
                    
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

I_mat = abs(I_mat);

x_axiz = x_range(1):(x_range(2)-x_range(1))/(num_xAxis-1):x_range(2);
y_axiz = y_range(1):(y_range(2)-y_range(1))/(num_yAxis-1):y_range(2);

subplot(2,2,2);
subtitle("BP法计算点的坐标")
set(gcf,'Position',[250,250,650,400]);
pcolor(x_axiz,y_axiz, I_mat');
axis equal
xlabel('X轴');
ylabel('Y轴');
shading flat
xlim(x_range);
ylim(y_range);
%% 
data=s_dc_all_cell{1}(1:4,1:512,10);
sample_rate = 10000e3  ; %采样率
frequency_slope = 60.012e12  ; %斜率
c = 3e8  ; %光速
sample_size = 512  ; %采样个数
chip_p_frame = 1  ; %每个frame内的chip个数
start_frequency = 77e9  ; %起始频率
idle_time = 100 / 1e6;
ramp_end_time = 60 / 1e6;
f_peridic = 24 / 1e3;
c_peridic = idle_time + ramp_end_time  ; %1chirp持续的时间
cc_peridic = c_peridic * 4  ; %相同发射天线chirp间隔时间
dx = 0.65 / 1e3  ; %虚拟天线距离

%对于第a1个chip：
r=(fft(data,512,2));
a=fftshift(fftshift(fft(r,200,1),1),2);


val_Angle=[-1.4985015  -1.48351648 -1.46853147 -1.45354645 -1.43856144 -1.42357642 -1.40859141 -1.39360639 -1.37862138 -1.36363636 -1.34865135 -1.33366633 -1.31868132 -1.3036963  -1.28871129 -1.27372627 -1.25874126 -1.24375624 -1.22877123 -1.21378621 -1.1988012  -1.18381618 -1.16883117 -1.15384615 -1.13886114 -1.12387612 -1.10889111 -1.09390609 -1.07892108 -1.06393606 -1.04895105 -1.03396603 -1.01898102 -1.003996   -0.98901099 -0.97402597 -0.95904096 -0.94405594 -0.92907093 -0.91408591 -0.8991009  -0.88411588 -0.86913087 -0.85414585 -0.83916084 -0.82417582 -0.80919081 -0.79420579 -0.77922078 -0.76423576 -0.74925075 -0.73426573 -0.71928072 -0.7042957 -0.68931069 -0.67432567 -0.65934066 -0.64435564 -0.62937063 -0.61438561 -0.5994006  -0.58441558 -0.56943057 -0.55444555 -0.53946054 -0.52447552 -0.50949051 -0.49450549 -0.47952048 -0.46453546 -0.44955045 -0.43456543 -0.41958042 -0.4045954  -0.38961039 -0.37462537 -0.35964036 -0.34465534 -0.32967033 -0.31468531 -0.2997003  -0.28471528 -0.26973027 -0.25474525 -0.23976024 -0.22477522 -0.20979021 -0.19480519 -0.17982018 -0.16483516 -0.14985015 -0.13486513 -0.11988012 -0.1048951  -0.08991009 -0.07492507 -0.05994006 -0.04495504 -0.02997003 -0.01498501  0.          0.01498501  0.02997003  0.04495504  0.05994006  0.07492507  0.08991009  0.1048951  0.11988012  0.13486513  0.14985015  0.16483516  0.17982018  0.19480519  0.20979021  0.22477522  0.23976024  0.25474525  0.26973027  0.28471528  0.2997003   0.31468531  0.32967033  0.34465534  0.35964036  0.37462537  0.38961039  0.4045954   0.41958042  0.43456543  0.44955045  0.46453546  0.47952048  0.49450549  0.50949051  0.52447552  0.53946054  0.55444555  0.56943057  0.58441558  0.5994006   0.61438561  0.62937063  0.64435564  0.65934066  0.67432567  0.68931069  0.7042957   0.71928072  0.73426573  0.74925075  0.76423576  0.77922078  0.79420579  0.80919081  0.82417582  0.83916084  0.85414585  0.86913087  0.88411588  0.8991009   0.91408591  0.92907093  0.94405594  0.95904096  0.97402597  0.98901099  1.003996  1.01898102  1.03396603  1.04895105  1.06393606  1.07892108  1.09390609  1.10889111  1.12387612  1.13886114  1.15384615  1.16883117  1.18381618  1.1988012   1.21378621  1.22877123  1.24375624  1.25874126  1.27372627  1.28871129  1.3036963   1.31868132  1.33366633  1.34865135  1.36363636  1.37862138  1.39360639  1.40859141  1.42357642  1.43856144  1.45354645  1.46853147  1.48351648];

val_R=[-12.4975005  -12.44868214 -12.39986378 -12.35104542 -12.30222705 -12.25340869 -12.20459033 -12.15577197 -12.10695361 -12.05813525 -12.00931689 -11.96049853 -11.91168016 -11.8628618  -11.81404344 -11.76522508 -11.71640672 -11.66758836 -11.61877    -11.56995163 -11.52113327 -11.47231491 -11.42349655 -11.37467819 -11.32585983 -11.27704147 -11.22822311 -11.17940474 -11.13058638 -11.08176802 -11.03294966 -10.9841313  -10.93531294 -10.88649458 -10.83767621 -10.78885785 -10.74003949 -10.69122113 -10.64240277 -10.59358441 -10.54476605 -10.49594769 -10.44712932 -10.39831096 -10.3494926 -10.30067424 -10.25185588 -10.20303752 -10.15421916 -10.10540079 -10.05658243 -10.00776407  -9.95894571  -9.91012735  -9.86130899  -9.81249063  -9.76367227  -9.7148539   -9.66603554  -9.61721718  -9.56839882  -9.51958046  -9.4707621   -9.42194374  -9.37312537  -9.32430701  -9.27548865  -9.22667029  -9.17785193  -9.12903357  -9.08021521  -9.03139685  -8.98257848  -8.93376012  -8.88494176  -8.8361234   -8.78730504  -8.73848668  -8.68966832  -8.64084996  -8.59203159  -8.54321323  -8.49439487  -8.44557651  -8.39675815  -8.34793979  -8.29912143  -8.25030306  -8.2014847   -8.15266634  -8.10384798  -8.05502962  -8.00621126  -7.9573929   -7.90857454  -7.85975617  -7.81093781  -7.76211945  -7.71330109  -7.66448273  -7.61566437  -7.56684601  -7.51802764  -7.46920928  -7.42039092  -7.37157256  -7.3227542   -7.27393584  -7.22511748  -7.17629912  -7.12748075  -7.07866239  -7.02984403  -6.98102567  -6.93220731  -6.88338895  -6.83457059  -6.78575222  -6.73693386  -6.6881155  -6.63929714  -6.59047878  -6.54166042  -6.49284206  -6.4440237  -6.39520533  -6.34638697  -6.29756861  -6.24875025  -6.19993189  -6.15111353  -6.10229517  -6.0534768   -6.00465844  -5.95584008  -5.90702172  -5.85820336  -5.809385    -5.76056664  -5.71174828  -5.66292991  -5.61411155  -5.56529319  -5.51647483  -5.46765647  -5.41883811  -5.37001975  -5.32120138  -5.27238302  -5.22356466  -5.1747463   -5.12592794  -5.07710958  -5.02829122  -4.97947286  -4.93065449  -4.88183613  -4.83301777  -4.78419941  -4.73538105  -4.68656269  -4.63774433  -4.58892596  -4.5401076   -4.49128924  -4.44247088  -4.39365252  -4.34483416  -4.2960158   -4.24719744  -4.19837907  -4.14956071  -4.10074235  -4.05192399  -4.00310563  -3.95428727  -3.90546891  -3.85665054  -3.80783218  -3.75901382  -3.71019546  -3.6613771   -3.61255874  -3.56374038  -3.51492202  -3.46610365  -3.41728529  -3.36846693  -3.31964857  -3.27083021  -3.22201185  -3.17319349  -3.12437512  -3.07555676  -3.0267384  -2.97792004  -2.92910168  -2.88028332  -2.83146496  -2.7826466  -2.73382823  -2.68500987  -2.63619151  -2.58737315  -2.53855479  -2.48973643  -2.44091807  -2.39209971  -2.34328134  -2.29446298  -2.24564462  -2.19682626  -2.1480079   -2.09918954  -2.05037118  -2.00155281  -1.95273445  -1.90391609  -1.85509773  -1.80627937  -1.75746101  -1.70864265  -1.65982429  -1.61100592  -1.56218756  -1.5133692   -1.46455084  -1.41573248  -1.36691412  -1.31809576  -1.26927739  -1.22045903  -1.17164067  -1.12282231  -1.07400395  -1.02518559  -0.97636723  -0.92754887  -0.8787305   -0.82991214  -0.78109378  -0.73227542  -0.68345706  -0.6346387   -0.58582034  -0.53700197  -0.48818361  -0.43936525  -0.39054689  -0.34172853  -0.29291017  -0.24409181  -0.19527345  -0.14645508  -0.09763672  -0.04881836   0.           0.04881836   0.09763672   0.14645508   0.19527345   0.24409181   0.29291017   0.34172853   0.39054689   0.43936525   0.48818361   0.53700197   0.58582034   0.6346387   0.68345706   0.73227542   0.78109378   0.82991214   0.8787305   0.92754887   0.97636723   1.02518559   1.07400395   1.12282231   1.17164067   1.22045903   1.26927739   1.31809576   1.36691412   1.41573248   1.46455084   1.5133692    1.56218756   1.61100592   1.65982429   1.70864265   1.75746101   1.80627937   1.85509773   1.90391609   1.95273445   2.00155281   2.05037118   2.09918954   2.1480079    2.19682626   2.24564462   2.29446298   2.34328134   2.39209971   2.44091807   2.48973643   2.53855479   2.58737315   2.63619151   2.68500987   2.73382823   2.7826466    2.83146496   2.88028332   2.92910168   2.97792004   3.0267384    3.07555676   3.12437512   3.17319349   3.22201185   3.27083021   3.31964857   3.36846693   3.41728529   3.46610365   3.51492202   3.56374038   3.61255874   3.6613771    3.71019546   3.75901382   3.80783218   3.85665054   3.90546891   3.95428727   4.00310563   4.05192399   4.10074235   4.14956071   4.19837907   4.24719744   4.2960158   4.34483416   4.39365252   4.44247088   4.49128924   4.5401076   4.58892596   4.63774433   4.68656269   4.73538105   4.78419941   4.83301777   4.88183613   4.93065449   4.97947286   5.02829122   5.07710958   5.12592794   5.1747463    5.22356466   5.27238302   5.32120138   5.37001975   5.41883811   5.46765647   5.51647483   5.56529319   5.61411155   5.66292991   5.71174828   5.76056664   5.809385     5.85820336   5.90702172   5.95584008   6.00465844   6.0534768    6.10229517   6.15111353   6.19993189   6.24875025   6.29756861   6.34638697   6.39520533   6.4440237    6.49284206   6.54166042   6.59047878   6.63929714   6.6881155    6.73693386   6.78575222   6.83457059   6.88338895   6.93220731   6.98102567   7.02984403   7.07866239   7.12748075   7.17629912   7.22511748   7.27393584   7.3227542    7.37157256   7.42039092   7.46920928   7.51802764   7.56684601   7.61566437   7.66448273   7.71330109   7.76211945   7.81093781   7.85975617   7.90857454   7.9573929   8.00621126   8.05502962   8.10384798   8.15266634   8.2014847   8.25030306   8.29912143   8.34793979   8.39675815   8.44557651   8.49439487   8.54321323   8.59203159   8.64084996   8.68966832   8.73848668   8.78730504   8.8361234    8.88494176   8.93376012   8.98257848   9.03139685   9.08021521   9.12903357   9.17785193   9.22667029   9.27548865   9.32430701   9.37312537   9.42194374   9.4707621    9.51958046   9.56839882   9.61721718   9.66603554   9.7148539    9.76367227   9.81249063   9.86130899   9.91012735   9.95894571  10.00776407  10.05658243  10.10540079  10.15421916  10.20303752  10.25185588  10.30067424  10.3494926   10.39831096  10.44712932  10.49594769  10.54476605  10.59358441  10.64240277  10.69122113  10.74003949  10.78885785  10.83767621  10.88649458 10.93531294  10.9841313   11.03294966  11.08176802  11.13058638  11.17940474  11.22822311  11.27704147  11.32585983  11.37467819  11.42349655  11.47231491  11.52113327  11.56995163  11.61877 11.66758836  11.71640672  11.76522508  11.81404344  11.8628618 11.91168016  11.96049853  12.00931689  12.05813525  12.10695361 12.15577197  12.20459033  12.25340869  12.30222705  12.35104542 12.39986378  12.44868214];
[rho, theta_dis] = meshgrid(val_R, val_Angle);
[x, y] = pol2cart(theta_dis, rho);
subplot(2,2,3);
%figure;
subtitle("1发四收FFT计算点的坐标")
for i =1:512
    for j=1:200
        if x(j,i)<0
            y(j,i)=0;
            x(j,i)=0;
            a(j,i)=0;
        end
        y(j,i)=-y(j,i);
    end
end

pcolor(x,y,abs(a) );
axis equal;
xlim([0 10.8])
ylim([-3.75 3.75])
xlabel('x方向距离(m)');
ylabel('y方向距离(m)');

shading flat
subplot(2, 2, 4);
%lfigure;
subtitle("1发四收FFT计算点的坐标")
pcolor(x,y,abs(a) );
axis equal;
xlim([min(min(x)) max(max(x))])
ylim([min(min(y)) max(max(y))])
xlabel('x方向距离(m)');
ylabel('y方向距离(m)');
title('t1')
shading flat
%imagesc(y)
%% 
figure
imagesc(abs(a))
