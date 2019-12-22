clear all
close all
clc
tic


% H = genH(256, 512);
load matrixH.mat H
ind = find(H==1);
[r,c] = ind2sub(size(H),ind);
[rows,cols] = size(H);
h = sparse(H);
n = cols;
k = n-rows;
R = k/n;

trials=[10 100 300 1200 40];


dB = [1, 1.5, 2, 2.5, 3];               % SNR的范围（dB）
SNRpbit = 10.^(dB/10);
No_uncoded = 1./SNRpbit;
No = No_uncoded./R;                     % Eb=1 ebno = Eb/N0 (energy per bit to noise power spectral density ratio)


iter=10;                                %迭代次数
maximum_blockerror=100;                 % 对于每个SNR下达到的最大的帧错误数
FER1=zeros(1,length(dB));
BER1=zeros(1,length(dB));
FER2=zeros(1,length(dB));
BER2=zeros(1,length(dB));
FER3=zeros(1,length(dB));
BER3=zeros(1,length(dB));
FER4=zeros(1,length(dB));
BER4=zeros(1,length(dB));
%%%%%%%%%%%%%%%%%%采用比特翻转译码算法%%%%%%%%%%%%%%%%%%%%%%%
for z=1:length(SNRpbit) % 对SNR循环
    rand('seed',584);
    randn('seed',843);
    biterrors1=0;
    blockerrors1=0;
    block1=0;
    biterrors2=0;
    blockerrors2=0;
    block2=0;
    biterrors3=0;
    blockerrors3=0;
    block3=0;
    biterrors4=0;
    blockerrors4=0;
    block4=0;
    for i = 1:trials(z)
    % while (blockerrors1 < maximum_blockerror)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 对消息进行编码 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        s=round(rand(1, cols-rows));%随机产生一个消息序列
        %使用H矩阵进行LDPC编码
        [u,P,rearranged_cols]=ldpc_encode(s,H);
        %对编码序列进行bpsk调制
        tx_waveform=bpsk(u);
        %将调制好的序列发送如AWGN信道
        N0 = No(z)/2;
        sigma=sqrt(N0);
        rx_waveform=tx_waveform + sigma*randn(1,length(tx_waveform));
        tx=rx_waveform';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%译码%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %比特翻转译码算法
        vhat1 = decodeBitFlip(tx, H, iter);
        uhat1=vhat1;
        [num1, rat1] = biterr(vhat1, u);
        if num1 ~= 0
            blockerrors1 = blockerrors1 + 1;
        end
        biterrors1 = biterrors1 + num1;
        block1 = block1+1;
        
        
        %BP译码算法
        vhat2 = decodeProbDomain(tx, H, N0, iter);
        uhat2 = vhat2;
        [num2, rat2] = biterr(vhat2, u);
        if num2 ~= 0
            blockerrors2 = blockerrors2 + 1;
        end   
        biterrors2 = biterrors2 + num2;
        block2 = block2+1;
        
        %最小和译码算法
        vhat3 = decodeMinSumLogDomain(tx, H, N0, iter);
        uhat3 = vhat3;
        [num3, rat3] = biterr(vhat3, u);
        if num3 ~= 0
            blockerrors3 = blockerrors3 + 1;
        end   
        biterrors3 = biterrors3 + num3;
        block3 = block3+1;
        
        %归一化最小和译码算法
        vhat4 =decodeNormalMinSumLogDomain(tx, H, N0, iter);
        uhat4 = vhat4;
        [num4, rat4] = biterr(vhat4, u);
         if num4 ~= 0
            blockerrors4 = blockerrors4 + 1;
        end   
        biterrors4 = biterrors4 + num4;
        block4 = block4+1;
    end
    BER1(z) = biterrors1/(block1*n);
    FER1(z) = blockerrors1/block1;
    BER2(z) = biterrors2/(block2*n);
    FER2(z) = blockerrors2/block2;
    BER3(z) = biterrors3/(block3*n);
    FER3(z) = blockerrors3/block3;
    BER4(z) = biterrors4/(block4*n);
    FER4(z) = blockerrors4/block4;
end
    
figure(1)
semilogy(dB,BER1,'r-o')
hold on
semilogy(dB,BER2,'g-o')
hold on
semilogy(dB,BER3,'b-o')
hold on
semilogy(dB,BER4,'y-o')
hold on
semilogy(dB,FER1,'r-*')
hold on
semilogy(dB,FER2,'g-*')
hold on
semilogy(dB,FER3,'b-*')
hold on
semilogy(dB,FER4,'y-*')
hold on
title('误比特率/误分组率')
ylabel('误比特率/误分组率')
xlabel('信噪比 (dB)')
grid

save('errors.mat', 'BER1', 'BER2', 'BER3', 'BER4', 'FER1', 'FER2', 'FER3', 'FER4', ...
    'biterrors1', 'biterrors2', 'biterrors3', 'biterrors4', 'block1', 'block2', 'block3', 'block4', ...
    'blockerrors1', 'blockerrors2', 'blockerrors3', 'blockerrors4')


toc