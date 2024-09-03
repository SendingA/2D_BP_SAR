clear;
% close all;

% ����ԭʼ�״�����
totalPath = 'D:\��ҵ\22-������\sar����\code_update\2D-BP-SAR-20230824\code\0812\matlab\data\';
rawData_FrameName = 'adc_data_cs_lr_4_l.bin';

savePath = 'D:\��ҵ\22-������\sar����\code_update\2D-BP-SAR-20230824\code\0812\matlab\data1\';

rawData = read1443Data([totalPath, rawData_FrameName]);

% ѡ������
detectRange = 20;


%% ����
f_c = 77*10^9;
K = 60.012*10^12;
chirp_typeNum = 1;
chirp_num = 1;
sample_num = 512;
frame_num = 3000;

radarNum = 1;
TxNum = 1;
RxNum = 4;

shiftSample2Frame = chirp_num *chirp_typeNum *sample_num;

using_rawData = zeros(RxNum, sample_num, frame_num);

rwin = hamming(sample_num);
rwinMat = zeros(RxNum, length(rwin));

for i_rxNum = 1:RxNum
    rwinMat(i_rxNum, :) = rwin.';
end
for i_frame = 1:frame_num
    using_rawData(:,:,i_frame) = rawData(:, ...
        shiftSample2Frame*(i_frame-1)+1: shiftSample2Frame*i_frame).*rwinMat;
end

% for i_frame = 1:frame_num
%     using_rawData(:,:,i_frame) = rawData(:, ...
%         shiftSample2Frame*(i_frame-1)+1: shiftSample2Frame*i_frame);
% end

dataForFFT = zeros(sample_num, frame_num);
for i_frame = 1:frame_num
   dataForFFT(:, i_frame) = using_rawData(3, :,i_frame).'; 
end

rangeFFT = fft(dataForFFT);

figure;
imagesc(abs(rangeFFT(10:270,:)));

figure;
plot(abs(rangeFFT(detectRange,:)));

% ��������
PointNum = 30; % �״���ڵ�λ��
shiftFrame = 32.35; %ÿ������֡
startFrame1 = 37;  %��ʼ֡

selectFrameMat = ones(1, PointNum);
for i_frame1 = 1:PointNum
    selectFrameMat(i_frame1) = startFrame1 + round(shiftFrame*(i_frame1-1));
end

hold on
scatter(selectFrameMat, abs(rangeFFT(detectRange, selectFrameMat)));

sarData = using_rawData(:,:, selectFrameMat);

save( [savePath, rawData_FrameName(1:end-4), '.mat'],'sarData');





