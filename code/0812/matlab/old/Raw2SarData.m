clear;
% close all;

% 导入原始雷达数据
totalPath = 'D:\fgf\00可能的新方向\000test\20230515\SAR\实验\data\20230729\rawData\';
rawData_FrameName = 'adc_data_emptyLow.bin';

savePath = 'D:\fgf\00可能的新方向\000test\20230515\SAR\实验\data\20230729\minData\';

rawData = read1443Data([totalPath, rawData_FrameName]);

detectRange = 225;

%% 参数
f_c = 77*10^9;
K = 60.012*10^12;
% K = 50.018*10^12;
frame = 1000;
chirp_typeNum = 1;
chirp_num = 1;
sample_num = 512;
frame_num = 1100;

radarNum = 1;
TxNum = 3;
RxNum = 4;

shiftSample2Frame = chirp_num *chirp_typeNum *sample_num;

using_rawData = zeros(RxNum, sample_num, frame_num);

for i_frame = 1:frame_num
    using_rawData(:,:,i_frame) = rawData(:, ...
        shiftSample2Frame*(i_frame-1)+1: shiftSample2Frame*i_frame);
end

dataForFFT = zeros(sample_num, frame_num);
for i_frame = 1:frame_num
   dataForFFT(:, i_frame) = using_rawData(3, :,i_frame).'; 
end

rangeFFT = fft(dataForFFT);

figure;
imagesc(abs(rangeFFT(:,:)));

figure;
plot(abs(rangeFFT(detectRange,:)));

lineNum1 = 31;
lineNum2 = 31;
shiftFrame = 15.4;
startFrame1 = 74;

% selectFrameMat(1) = startF
selectFrameMat = ones(1, lineNum1+lineNum2);
for i_frame1 = 1:lineNum1
    selectFrameMat(i_frame1) = startFrame1 + round(shiftFrame*(i_frame1-1));
end

startFrame2 = selectFrameMat(i_frame1)+ 19;
for i_frame2 = 1:lineNum2
    selectFrameMat(i_frame1+i_frame2) = startFrame2 + round(shiftFrame*(i_frame2-1));
end

hold on
scatter(selectFrameMat, abs(rangeFFT(detectRange, selectFrameMat)));

sarData = using_rawData(:,:, selectFrameMat);

% savePath

save( [savePath, rawData_FrameName(1:end-4), '.mat'],'sarData');








