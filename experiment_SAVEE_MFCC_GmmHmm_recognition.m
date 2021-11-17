clc; clear; close all;

%% Load SAVEE MFCC
oldpath = path;
path(oldpath, '.\MAT file');
load('SAVEE_MFCC.mat');



%% Sequence Pre-proccessing
for i = 1:length(MFCC_NeutDC)
	MFCC_NeutDC{i} = MFCC_NeutDC{i}{1}';
end
for i = 1:length(MFCC_NeutJE)
	MFCC_NeutJE{i} = MFCC_NeutJE{i}{1}';
end
for i = 1:length(MFCC_NeutJK)
	MFCC_NeutJK{i} = MFCC_NeutJK{i}{1}';
end
for i = 1:length(MFCC_NeutKL)
	MFCC_NeutKL{i} = MFCC_NeutKL{i}{1}';
end
for i = 1:length(MFCC_HpyDC)
	MFCC_HpyDC{i} = MFCC_HpyDC{i}{1}';
end
for i = 1:length(MFCC_HpyJE)
	MFCC_HpyJE{i} = MFCC_HpyJE{i}{1}';
end
for i = 1:length(MFCC_HpyJK)
	MFCC_HpyJK{i} = MFCC_HpyJK{i}{1}';
end
for i = 1:length(MFCC_HpyKL)
	MFCC_HpyKL{i} = MFCC_HpyKL{i}{1}';
end
for i = 1:length(MFCC_SadDC)
	MFCC_SadDC{i} = MFCC_SadDC{i}{1}';
end
for i = 1:length(MFCC_SadJE)
	MFCC_SadJE{i} = MFCC_SadJE{i}{1}';
end
for i = 1:length(MFCC_SadJK)
	MFCC_SadJK{i} = MFCC_SadJK{i}{1}';
end
for i = 1:length(MFCC_SadKL)
	MFCC_SadKL{i} = MFCC_SadKL{i}{1}';
end
for i = 1:length(MFCC_AgyDC)
	MFCC_AgyDC{i} = MFCC_AgyDC{i}{1}';
end
for i = 1:length(MFCC_AgyJE)
	MFCC_AgyJE{i} = MFCC_AgyJE{i}{1}';
end
for i = 1:length(MFCC_AgyJK)
	MFCC_AgyJK{i} = MFCC_AgyJK{i}{1}';
end
for i = 1:length(MFCC_AgyKL)
	MFCC_AgyKL{i} = MFCC_AgyKL{i}{1}';
end
for i = 1:length(MFCC_FearDC)
	MFCC_FearDC{i} = MFCC_FearDC{i}{1}';
end
for i = 1:length(MFCC_FearJE)
	MFCC_FearJE{i} = MFCC_FearJE{i}{1}';
end
for i = 1:length(MFCC_FearJK)
	MFCC_FearJK{i} = MFCC_FearJK{i}{1}';
end
for i = 1:length(MFCC_FearKL)
	MFCC_FearKL{i} = MFCC_FearKL{i}{1}';
end
for i = 1:length(MFCC_DisgDC)
	MFCC_DisgDC{i} = MFCC_DisgDC{i}{1}';
end
for i = 1:length(MFCC_DisgJE)
	MFCC_DisgJE{i} = MFCC_DisgJE{i}{1}';
end
for i = 1:length(MFCC_DisgJK)
	MFCC_DisgJK{i} = MFCC_DisgJK{i}{1}';
end
for i = 1:length(MFCC_DisgKL)
	MFCC_DisgKL{i} = MFCC_DisgKL{i}{1}';
end
for i = 1:length(MFCC_SurpDC)
	MFCC_SurpDC{i} = MFCC_SurpDC{i}{1}';
end
for i = 1:length(MFCC_SurpJE)
	MFCC_SurpJE{i} = MFCC_SurpJE{i}{1}';
end
for i = 1:length(MFCC_SurpJK)
	MFCC_SurpJK{i} = MFCC_SurpJK{i}{1}';
end
for i = 1:length(MFCC_SurpKL)
	MFCC_SurpKL{i} = MFCC_SurpKL{i}{1}';
end



%% Partition Data
N_Training = 5;
N_Training_Neut = 5;
% N_Training = 6;
% N_Training_Neut = 6;
% N_Training = 7;
% N_Training_Neut = 7;
% N_Training = 14;
% N_Training_Neut = 29;

seq = randperm(size(MFCC_NeutDC, 2));
Train_NeutDC = MFCC_NeutDC(seq(1:N_Training_Neut));
Test_NeutDC = MFCC_NeutDC(seq(N_Training_Neut+1:end));
seq = randperm(size(MFCC_NeutJE, 2));
Train_NeutJE = MFCC_NeutJE(seq(1:N_Training_Neut));
Test_NeutJE = MFCC_NeutJE(seq(N_Training_Neut+1:end));
seq = randperm(size(MFCC_NeutJK, 2));
Train_NeutJK = MFCC_NeutJK(seq(1:N_Training_Neut));
Test_NeutJK = MFCC_NeutJK(seq(N_Training_Neut+1:end));
seq = randperm(size(MFCC_NeutKL, 2));
Train_NeutKL = MFCC_NeutKL(seq(1:N_Training_Neut));
Test_NeutKL = MFCC_NeutKL(seq(N_Training_Neut+1:end));

seq = randperm(size(MFCC_HpyDC, 2));
Train_HpyDC = MFCC_HpyDC(seq(1:N_Training));
Test_HpyDC = MFCC_HpyDC(seq(N_Training+1:end));
seq = randperm(size(MFCC_HpyJE, 2));
Train_HpyJE = MFCC_HpyJE(seq(1:N_Training));
Test_HpyJE = MFCC_HpyJE(seq(N_Training+1:end));
seq = randperm(size(MFCC_HpyJK, 2));
Train_HpyJK = MFCC_HpyJK(seq(1:N_Training));
Test_HpyJK = MFCC_HpyJK(seq(N_Training+1:end));
seq = randperm(size(MFCC_HpyKL, 2));
Train_HpyKL = MFCC_HpyKL(seq(1:N_Training));
Test_HpyKL = MFCC_HpyKL(seq(N_Training+1:end));

seq = randperm(size(MFCC_SadDC, 2));
Train_SadDC = MFCC_SadDC(seq(1:N_Training));
Test_SadDC = MFCC_SadDC(seq(N_Training+1:end));
seq = randperm(size(MFCC_SadJE, 2));
Train_SadJE = MFCC_SadJE(seq(1:N_Training));
Test_SadJE = MFCC_SadJE(seq(N_Training+1:end));
seq = randperm(size(MFCC_SadJK, 2));
Train_SadJK = MFCC_SadJK(seq(1:N_Training));
Test_SadJK = MFCC_SadJK(seq(N_Training+1:end));
seq = randperm(size(MFCC_SadKL, 2));
Train_SadKL = MFCC_SadKL(seq(1:N_Training));
Test_SadKL = MFCC_SadKL(seq(N_Training+1:end));

seq = randperm(size(MFCC_AgyDC, 2));
Train_AgyDC = MFCC_AgyDC(seq(1:N_Training));
Test_AgyDC = MFCC_AgyDC(seq(N_Training+1:end));
seq = randperm(size(MFCC_AgyJE, 2));
Train_AgyJE = MFCC_AgyJE(seq(1:N_Training));
Test_AgyJE = MFCC_AgyJE(seq(N_Training+1:end));
seq = randperm(size(MFCC_AgyJK, 2));
Train_AgyJK = MFCC_AgyJK(seq(1:N_Training));
Test_AgyJK = MFCC_AgyJK(seq(N_Training+1:end));
seq = randperm(size(MFCC_AgyKL, 2));
Train_AgyKL = MFCC_AgyKL(seq(1:N_Training));
Test_AgyKL = MFCC_AgyKL(seq(N_Training+1:end));

seq = randperm(size(MFCC_FearDC, 2));
Train_FearDC = MFCC_FearDC(seq(1:N_Training));
Test_FearDC = MFCC_FearDC(seq(N_Training+1:end));
seq = randperm(size(MFCC_FearJE, 2));
Train_FearJE = MFCC_FearJE(seq(1:N_Training));
Test_FearJE = MFCC_FearJE(seq(N_Training+1:end));
seq = randperm(size(MFCC_FearJK, 2));
Train_FearJK = MFCC_FearJK(seq(1:N_Training));
Test_FearJK = MFCC_FearJK(seq(N_Training+1:end));
seq = randperm(size(MFCC_FearKL, 2));
Train_FearKL = MFCC_FearKL(seq(1:N_Training));
Test_FearKL = MFCC_FearKL(seq(N_Training+1:end));

seq = randperm(size(MFCC_DisgDC, 2));
Train_DisgDC = MFCC_DisgDC(seq(1:N_Training));
Test_DisgDC = MFCC_DisgDC(seq(N_Training+1:end));
seq = randperm(size(MFCC_DisgJE, 2));
Train_DisgJE = MFCC_DisgJE(seq(1:N_Training));
Test_DisgJE = MFCC_DisgJE(seq(N_Training+1:end));
seq = randperm(size(MFCC_DisgJK, 2));
Train_DisgJK = MFCC_DisgJK(seq(1:N_Training));
Test_DisgJK = MFCC_DisgJK(seq(N_Training+1:end));
seq = randperm(size(MFCC_DisgKL, 2));
Train_DisgKL = MFCC_DisgKL(seq(1:N_Training));
Test_DisgKL = MFCC_DisgKL(seq(N_Training+1:end));

seq = randperm(size(MFCC_SurpDC, 2));
Train_SurpDC = MFCC_SurpDC(seq(1:N_Training));
Test_SurpDC = MFCC_SurpDC(seq(N_Training+1:end));
seq = randperm(size(MFCC_SurpJE, 2));
Train_SurpJE = MFCC_SurpJE(seq(1:N_Training));
Test_SurpJE = MFCC_SurpJE(seq(N_Training+1:end));
seq = randperm(size(MFCC_SurpJK, 2));
Train_SurpJK = MFCC_SurpJK(seq(1:N_Training));
Test_SurpJK = MFCC_SurpJK(seq(N_Training+1:end));
seq = randperm(size(MFCC_SurpKL, 2));
Train_SurpKL = MFCC_SurpKL(seq(1:N_Training));
Test_SurpKL = MFCC_SurpKL(seq(N_Training+1:end));
   


%% Compute Mean of Trainset
temp = [];
for i = length(Train_NeutDC)
	temp = [temp Train_NeutDC{i}];
end
mean_NeutDC = mean(temp, 2);
temp = [];
for i = length(Train_NeutJE)
	temp = [temp Train_NeutJE{i}];
end
mean_NeutJE = mean(temp, 2);
temp = [];
for i = length(Train_NeutJK)
	temp = [temp Train_NeutJK{i}];
end
mean_NeutJK = mean(temp, 2);
temp = [];
for i = length(Train_NeutKL)
	temp = [temp Train_NeutKL{i}];
end
mean_NeutKL = mean(temp, 2);

temp = [];
for i = length(Train_HpyDC)
	temp = [temp Train_HpyDC{i}];
end
mean_HpyDC = mean(temp, 2);
temp = [];
for i = length(Train_HpyJE)
	temp = [temp Train_HpyJE{i}];
end
mean_HpyJE = mean(temp, 2);
temp = [];
for i = length(Train_HpyJK)
	temp = [temp Train_HpyJK{i}];
end
mean_HpyJK = mean(temp, 2);
temp = [];
for i = length(Train_HpyKL)
	temp = [temp Train_HpyKL{i}];
end
mean_HpyKL = mean(temp, 2);

temp = [];
for i = length(Train_SadDC)
	temp = [temp Train_SadDC{i}];
end
mean_SadDC = mean(temp, 2);
temp = [];
for i = length(Train_SadJE)
	temp = [temp Train_SadJE{i}];
end
mean_SadJE = mean(temp, 2);
temp = [];
for i = length(Train_SadJK)
	temp = [temp Train_SadJK{i}];
end
mean_SadJK = mean(temp, 2);
temp = [];
for i = length(Train_SadKL)
	temp = [temp Train_SadKL{i}];
end
mean_SadKL = mean(temp, 2);

temp = [];
for i = length(Train_AgyDC)
	temp = [temp Train_AgyDC{i}];
end
mean_AgyDC = mean(temp, 2);
temp = [];
for i = length(Train_AgyJE)
	temp = [temp Train_AgyJE{i}];
end
mean_AgyJE = mean(temp, 2);
temp = [];
for i = length(Train_AgyJK)
	temp = [temp Train_AgyJK{i}];
end
mean_AgyJK = mean(temp, 2);
temp = [];
for i = length(Train_AgyKL)
	temp = [temp Train_AgyKL{i}];
end
mean_AgyKL = mean(temp, 2);

temp = [];
for i = length(Train_FearDC)
	temp = [temp Train_FearDC{i}];
end
mean_FearDC = mean(temp, 2);
temp = [];
for i = length(Train_FearJE)
	temp = [temp Train_FearJE{i}];
end
mean_FearJE = mean(temp, 2);
temp = [];
for i = length(Train_FearJK)
	temp = [temp Train_FearJK{i}];
end
mean_FearJK = mean(temp, 2);
temp = [];
for i = length(Train_FearKL)
	temp = [temp Train_FearKL{i}];
end
mean_FearKL = mean(temp, 2);

temp = [];
for i = length(Train_DisgDC)
	temp = [temp Train_DisgDC{i}];
end
mean_DisgDC = mean(temp, 2);
temp = [];
for i = length(Train_DisgJE)
	temp = [temp Train_DisgJE{i}];
end
mean_DisgJE = mean(temp, 2);
temp = [];
for i = length(Train_DisgJK)
	temp = [temp Train_DisgJK{i}];
end
mean_DisgJK = mean(temp, 2);
temp = [];
for i = length(Train_DisgKL)
	temp = [temp Train_DisgKL{i}];
end
mean_DisgKL = mean(temp, 2);

temp = [];
for i = length(Train_SurpDC)
	temp = [temp Train_SurpDC{i}];
end
mean_SurpDC = mean(temp, 2);
temp = [];
for i = length(Train_SurpJE)
	temp = [temp Train_SurpJE{i}];
end
mean_SurpJE = mean(temp, 2);
temp = [];
for i = length(Train_SurpJK)
	temp = [temp Train_SurpJK{i}];
end
mean_SurpJK = mean(temp, 2);
temp = [];
for i = length(Train_SurpKL)
	temp = [temp Train_SurpKL{i}];
end
mean_SurpKL = mean(temp, 2);



%% Compute Mean of beginning of Trainset
temp = [];
for i = 1:length(Train_NeutDC)
	temp = [temp Train_NeutDC{i}(:, 1)];
end
mean_beginNeutDC = mean(temp, 2);
temp = [];
for i = 1:length(Train_NeutJE)
	temp = [temp Train_NeutJE{i}(:, 1)];
end
mean_beginNeutJE = mean(temp, 2);
temp = [];
for i = 1:length(Train_NeutJK)
	temp = [temp Train_NeutJK{i}(:, 1)];
end
mean_beginNeutJK = mean(temp, 2);
temp = [];
for i = 1:length(Train_NeutKL)
	temp = [temp Train_NeutKL{i}(:, 1)];
end
mean_beginNeutKL = mean(temp, 2);

temp = [];
for i = 1:length(Train_HpyDC)
	temp = [temp Train_HpyDC{i}(:, 1)];
end
mean_beginHpyDC = mean(temp, 2);
temp = [];
for i = 1:length(Train_HpyJE)
	temp = [temp Train_HpyJE{i}(:, 1)];
end
mean_beginHpyJE = mean(temp, 2);
temp = [];
for i = 1:length(Train_HpyJK)
	temp = [temp Train_HpyJK{i}(:, 1)];
end
mean_beginHpyJK = mean(temp, 2);
temp = [];
for i = 1:length(Train_HpyKL)
	temp = [temp Train_HpyKL{i}(:, 1)];
end
mean_beginHpyKL = mean(temp, 2);

temp = [];
for i = 1:length(Train_SadDC)
	temp = [temp Train_SadDC{i}(:, 1)];
end
mean_beginSadDC = mean(temp, 2);
temp = [];
for i = 1:length(Train_SadJE)
	temp = [temp Train_SadJE{i}(:, 1)];
end
mean_beginSadJE = mean(temp, 2);
temp = [];
for i = 1:length(Train_SadJK)
	temp = [temp Train_SadJK{i}(:, 1)];
end
mean_beginSadJK = mean(temp, 2);
temp = [];
for i = 1:length(Train_SadKL)
	temp = [temp Train_SadKL{i}(:, 1)];
end
mean_beginSadKL = mean(temp, 2);

temp = [];
for i = 1:length(Train_AgyDC)
	temp = [temp Train_AgyDC{i}(:, 1)];
end
mean_beginAgyDC = mean(temp, 2);
temp = [];
for i = 1:length(Train_AgyJE)
	temp = [temp Train_AgyJE{i}(:, 1)];
end
mean_beginAgyJE = mean(temp, 2);
temp = [];
for i = 1:length(Train_AgyJK)
	temp = [temp Train_AgyJK{i}(:, 1)];
end
mean_beginAgyJK = mean(temp, 2);
temp = [];
for i = 1:length(Train_AgyKL)
	temp = [temp Train_AgyKL{i}(:, 1)];
end
mean_beginAgyKL = mean(temp, 2);

temp = [];
for i = 1:length(Train_FearDC)
	temp = [temp Train_FearDC{i}(:, 1)];
end
mean_beginFearDC = mean(temp, 2);
temp = [];
for i = 1:length(Train_FearJE)
	temp = [temp Train_FearJE{i}(:, 1)];
end
mean_beginFearJE = mean(temp, 2);
temp = [];
for i = 1:length(Train_FearJK)
	temp = [temp Train_FearJK{i}(:, 1)];
end
mean_beginFearJK = mean(temp, 2);
temp = [];
for i = 1:length(Train_FearKL)
	temp = [temp Train_FearKL{i}(:, 1)];
end
mean_beginFearKL = mean(temp, 2);

temp = [];
for i = 1:length(Train_DisgDC)
	temp = [temp Train_DisgDC{i}(:, 1)];
end
mean_beginDisgDC = mean(temp, 2);
temp = [];
for i = 1:length(Train_DisgJE)
	temp = [temp Train_DisgJE{i}(:, 1)];
end
mean_beginDisgJE = mean(temp, 2);
temp = [];
for i = 1:length(Train_DisgJK)
	temp = [temp Train_DisgJK{i}(:, 1)];
end
mean_beginDisgJK = mean(temp, 2);
temp = [];
for i = 1:length(Train_DisgKL)
	temp = [temp Train_DisgKL{i}(:, 1)];
end
mean_beginDisgKL = mean(temp, 2);

temp = [];
for i = 1:length(Train_SurpDC)
	temp = [temp Train_SurpDC{i}(:, 1)];
end
mean_beginSurpDC = mean(temp, 2);
temp = []; 
for i = 1:length(Train_SurpJE)
	temp = [temp Train_SurpJE{i}(:, 1)];
end
mean_beginSurpJE = mean(temp, 2);
temp = [];
for i = 1:length(Train_SurpJK)
	temp = [temp Train_SurpJK{i}(:, 1)];
end
mean_beginSurpJK = mean(temp, 2);
temp = [];
for i = 1:length(Train_SurpKL)
	temp = [temp Train_SurpKL{i}(:, 1)];
end
mean_beginSurpKL = mean(temp, 2);



%% Model Selection
% Set number of GMM components and number of states   
numStates = 2;
numComp = 2;

%numStates = 4;
%numComp = 4;

%numStates = 8;
%numComp = 8;



%% GMM Initialize
w_ini = repmat(1/numComp, [1 numComp]);

% produce offset array for GMM components
GMM_offset = [];
for comp = 1:numComp
    GMM_offset = [GMM_offset 0.5*(comp)];
end

temp_mean_NeutDC = mean_NeutDC * GMM_offset;
temp_mean_NeutJE = mean_NeutJE * GMM_offset;
temp_mean_NeutJK = mean_NeutJK * GMM_offset;
temp_mean_NeutKL = mean_NeutKL * GMM_offset;
temp_mean_HpyDC = mean_HpyDC * GMM_offset;
temp_mean_HpyJE = mean_HpyJE * GMM_offset;
temp_mean_HpyJK = mean_HpyJK * GMM_offset;
temp_mean_HpyKL = mean_HpyKL * GMM_offset;
temp_mean_SadDC = mean_SadDC * GMM_offset;
temp_mean_SadJE = mean_SadJE * GMM_offset;
temp_mean_SadJK = mean_SadJK * GMM_offset;
temp_mean_SadKL = mean_SadKL * GMM_offset;
temp_mean_AgyDC = mean_AgyDC * GMM_offset;
temp_mean_AgyJE = mean_AgyJE * GMM_offset;
temp_mean_AgyJK = mean_AgyJK * GMM_offset;
temp_mean_AgyKL = mean_AgyKL * GMM_offset;
temp_mean_FearDC = mean_FearDC * GMM_offset;
temp_mean_FearJE = mean_FearJE * GMM_offset;
temp_mean_FearJK = mean_FearJK * GMM_offset;
temp_mean_FearKL = mean_FearKL * GMM_offset;
temp_mean_DisgDC = mean_DisgDC * GMM_offset;
temp_mean_DisgJE = mean_DisgJE * GMM_offset;
temp_mean_DisgJK = mean_DisgJK * GMM_offset;
temp_mean_DisgKL = mean_DisgKL * GMM_offset;
temp_mean_SurpDC = mean_SurpDC * GMM_offset;
temp_mean_SurpJE = mean_SurpJE * GMM_offset;
temp_mean_SurpJK = mean_SurpJK * GMM_offset;
temp_mean_SurpKL = mean_SurpKL * GMM_offset;

cov_ini = repmat(diag(ones(1, numCoef)), [1,1,numComp]);



%% State Initialize
TRANS = repmat(1/numStates, [numStates numStates]);

% obtain distance between feature mean and feature beginning
distNeutDC = abs(mean_beginNeutDC - mean_NeutDC);
distNeutJE = abs(mean_beginNeutJE - mean_NeutJE);
distNeutJK = abs(mean_beginNeutJK - mean_NeutJK);
distNeutKL = abs(mean_beginNeutKL - mean_NeutKL);
distHpyDC = abs(mean_beginHpyDC - mean_HpyDC);
distHpyJE = abs(mean_beginHpyJE - mean_HpyJE);
distHpyJK = abs(mean_beginHpyJK - mean_HpyJK);
distHpyKL = abs(mean_beginHpyKL - mean_HpyKL);
distSadDC = abs(mean_beginSadDC - mean_SadDC);
distSadJE = abs(mean_beginSadJE - mean_SadJE);
distSadJK = abs(mean_beginSadJK - mean_SadJK);
distSadKL = abs(mean_beginSadKL - mean_SadKL);
distAgyDC = abs(mean_beginAgyDC - mean_AgyDC);
distAgyJE = abs(mean_beginAgyJE - mean_AgyJE);
distAgyJK = abs(mean_beginAgyJK - mean_AgyJK);
distAgyKL = abs(mean_beginAgyKL - mean_AgyKL);
distFearDC = abs(mean_beginFearDC - mean_FearDC);
distFearJE = abs(mean_beginFearJE - mean_FearJE);
distFearJK = abs(mean_beginFearJK - mean_FearJK);
distFearKL = abs(mean_beginFearKL - mean_FearKL);
distDisgDC = abs(mean_beginDisgDC - mean_DisgDC);
distDisgJE = abs(mean_beginDisgJE - mean_DisgJE);
distDisgJK = abs(mean_beginDisgJK - mean_DisgJK);
distDisgKL = abs(mean_beginDisgKL - mean_DisgKL);
distSurpDC = abs(mean_beginSurpDC - mean_SurpDC);
distSurpJE = abs(mean_beginSurpJE - mean_SurpJE);
distSurpJK = abs(mean_beginSurpJK - mean_SurpJK);
distSurpKL = abs(mean_beginSurpKL - mean_SurpKL);

mean_NeutDC_ini = zeros(numCoef, numComp, numStates);
for state = 1:numStates
    mean_NeutDC_ini(:, :, state) = temp_mean_NeutDC - (0.5*(state-1))*distNeutDC;
end
mean_NeutJE_ini = zeros(numCoef, numComp, numStates);
for state = 1:numStates
    mean_NeutJE_ini(:, :, state) = temp_mean_NeutJE - (0.5*(state-1))*distNeutJE;
end
mean_NeutJK_ini = zeros(numCoef, numComp, numStates);
for state = 1:numStates
    mean_NeutJK_ini(:, :, state) = temp_mean_NeutJK - (0.5*(state-1))*distNeutJK;
end
mean_NeutKL_ini = zeros(numCoef, numComp, numStates);
for state = 1:numStates
    mean_NeutKL_ini(:, :, state) = temp_mean_NeutKL - (0.5*(state-1))*distNeutKL;
end

mean_HpyDC_ini = zeros(numCoef, numComp, numStates);
for state = 1:numStates
    mean_HpyDC_ini(:, :, state) = temp_mean_HpyDC - (0.5*(state-1))*distHpyDC;
end
mean_HpyJE_ini = zeros(numCoef, numComp, numStates);
for state = 1:numStates
    mean_HpyJE_ini(:, :, state) = temp_mean_HpyJE - (0.5*(state-1))*distHpyJE;
end
mean_HpyJK_ini = zeros(numCoef, numComp, numStates);
for state = 1:numStates
    mean_HpyJK_ini(:, :, state) = temp_mean_HpyJK - (0.5*(state-1))*distHpyJK;
end
mean_HpyKL_ini = zeros(numCoef, numComp, numStates);
for state = 1:numStates
    mean_HpyKL_ini(:, :, state) = temp_mean_HpyKL - (0.5*(state-1))*distHpyKL;
end

mean_SadDC_ini = zeros(numCoef, numComp, numStates);
for state = 1:numStates
    mean_SadDC_ini(:, :, state) = temp_mean_SadDC - (0.5*(state-1))*distSadDC;
end
mean_SadJE_ini = zeros(numCoef, numComp, numStates);
for state = 1:numStates
    mean_SadJE_ini(:, :, state) = temp_mean_SadJE - (0.5*(state-1))*distSadJE;
end
mean_SadJK_ini = zeros(numCoef, numComp, numStates);
for state = 1:numStates
    mean_SadJK_ini(:, :, state) = temp_mean_SadJK - (0.5*(state-1))*distSadJK;
end
mean_SadKL_ini = zeros(numCoef, numComp, numStates);
for state = 1:numStates
    mean_SadKL_ini(:, :, state) = temp_mean_SadKL - (0.5*(state-1))*distSadKL;
end

mean_AgyDC_ini = zeros(numCoef, numComp, numStates);
for state = 1:numStates
    mean_AgyDC_ini(:, :, state) = temp_mean_AgyDC - (0.5*(state-1))*distAgyDC;
end
mean_AgyJE_ini = zeros(numCoef, numComp, numStates);
for state = 1:numStates
    mean_AgyJE_ini(:, :, state) = temp_mean_AgyJE - (0.5*(state-1))*distAgyJE;
end
mean_AgyJK_ini = zeros(numCoef, numComp, numStates);
for state = 1:numStates
    mean_AgyJK_ini(:, :, state) = temp_mean_AgyJK - (0.5*(state-1))*distAgyJK;
end
mean_AgyKL_ini = zeros(numCoef, numComp, numStates);
for state = 1:numStates
    mean_AgyKL_ini(:, :, state) = temp_mean_AgyKL - (0.5*(state-1))*distAgyKL;
end

mean_FearDC_ini = zeros(numCoef, numComp, numStates);
for state = 1:numStates
    mean_FearDC_ini(:, :, state) = temp_mean_FearDC - (0.5*(state-1))*distFearDC;
end
mean_FearJE_ini = zeros(numCoef, numComp, numStates);
for state = 1:numStates
    mean_FearJE_ini(:, :, state) = temp_mean_FearJE - (0.5*(state-1))*distFearJE;
end
mean_FearJK_ini = zeros(numCoef, numComp, numStates);
for state = 1:numStates
    mean_FearJK_ini(:, :, state) = temp_mean_FearJK - (0.5*(state-1))*distFearJK;
end
mean_FearKL_ini = zeros(numCoef, numComp, numStates);
for state = 1:numStates
    mean_FearKL_ini(:, :, state) = temp_mean_FearKL - (0.5*(state-1))*distFearKL;
end

mean_DisgDC_ini = zeros(numCoef, numComp, numStates);
for state = 1:numStates
    mean_DisgDC_ini(:, :, state) = temp_mean_DisgDC - (0.5*(state-1))*distDisgDC;
end
mean_DisgJE_ini = zeros(numCoef, numComp, numStates);
for state = 1:numStates
    mean_DisgJE_ini(:, :, state) = temp_mean_DisgJE - (0.5*(state-1))*distDisgJE;
end
mean_DisgJK_ini = zeros(numCoef, numComp, numStates);
for state = 1:numStates
    mean_DisgJK_ini(:, :, state) = temp_mean_DisgJK - (0.5*(state-1))*distDisgJK;
end
mean_DisgKL_ini = zeros(numCoef, numComp, numStates);
for state = 1:numStates
    mean_DisgKL_ini(:, :, state) = temp_mean_DisgKL - (0.5*(state-1))*distDisgKL;
end

mean_SurpDC_ini = zeros(numCoef, numComp, numStates);
for state = 1:numStates
    mean_SurpDC_ini(:, :, state) = temp_mean_SurpDC - (0.5*(state-1))*distSurpDC;
end
mean_SurpJE_ini = zeros(numCoef, numComp, numStates);
for state = 1:numStates
    mean_SurpJE_ini(:, :, state) = temp_mean_SurpJE - (0.5*(state-1))*distSurpJE;
end
mean_SurpJK_ini = zeros(numCoef, numComp, numStates);
for state = 1:numStates
    mean_SurpJK_ini(:, :, state) = temp_mean_SurpJK - (0.5*(state-1))*distSurpJK;
end
mean_SurpKL_ini = zeros(numCoef, numComp, numStates);
for state = 1:numStates
    mean_SurpKL_ini(:, :, state) = temp_mean_SurpKL - (0.5*(state-1))*distSurpKL;
end

cov_ini = repmat(cov_ini, [1,1,1,numStates]);



%% Combine GMM Parameters
for state = 1:numStates
	EMISGmm_NeutDC(state, :) = {w_ini mean_NeutDC_ini(:,:,state)' cov_ini(:, :, :, state)};
	EMISGmm_NeutJE(state, :) = {w_ini mean_NeutJE_ini(:,:,state)' cov_ini(:, :, :, state)};
	EMISGmm_NeutJK(state, :) = {w_ini mean_NeutJK_ini(:,:,state)' cov_ini(:, :, :, state)};
	EMISGmm_NeutKL(state, :) = {w_ini mean_NeutKL_ini(:,:,state)' cov_ini(:, :, :, state)};
	EMISGmm_HpyDC(state, :) = {w_ini mean_HpyDC_ini(:,:,state)' cov_ini(:, :, :, state)};
	EMISGmm_HpyJE(state, :) = {w_ini mean_HpyJE_ini(:,:,state)' cov_ini(:, :, :, state)};
	EMISGmm_HpyJK(state, :) = {w_ini mean_HpyJK_ini(:,:,state)' cov_ini(:, :, :, state)};
	EMISGmm_HpyKL(state, :) = {w_ini mean_HpyKL_ini(:,:,state)' cov_ini(:, :, :, state)};
	EMISGmm_SadDC(state, :) = {w_ini mean_SadDC_ini(:,:,state)' cov_ini(:, :, :, state)};
	EMISGmm_SadJE(state, :) = {w_ini mean_SadJE_ini(:,:,state)' cov_ini(:, :, :, state)};
	EMISGmm_SadJK(state, :) = {w_ini mean_SadJK_ini(:,:,state)' cov_ini(:, :, :, state)};
	EMISGmm_SadKL(state, :) = {w_ini mean_SadKL_ini(:,:,state)' cov_ini(:, :, :, state)};
	EMISGmm_AgyDC(state, :) = {w_ini mean_AgyDC_ini(:,:,state)' cov_ini(:, :, :, state)};
	EMISGmm_AgyJE(state, :) = {w_ini mean_AgyJE_ini(:,:,state)' cov_ini(:, :, :, state)};
	EMISGmm_AgyJK(state, :) = {w_ini mean_AgyJK_ini(:,:,state)' cov_ini(:, :, :, state)};
	EMISGmm_AgyKL(state, :) = {w_ini mean_AgyKL_ini(:,:,state)' cov_ini(:, :, :, state)};
	EMISGmm_FearDC(state, :) = {w_ini mean_FearDC_ini(:,:,state)' cov_ini(:, :, :, state)};
	EMISGmm_FearJE(state, :) = {w_ini mean_FearJE_ini(:,:,state)' cov_ini(:, :, :, state)};
	EMISGmm_FearJK(state, :) = {w_ini mean_FearJK_ini(:,:,state)' cov_ini(:, :, :, state)};
	EMISGmm_FearKL(state, :) = {w_ini mean_FearKL_ini(:,:,state)' cov_ini(:, :, :, state)};
	EMISGmm_DisgDC(state, :) = {w_ini mean_DisgDC_ini(:,:,state)' cov_ini(:, :, :, state)};
	EMISGmm_DisgJE(state, :) = {w_ini mean_DisgJE_ini(:,:,state)' cov_ini(:, :, :, state)};
	EMISGmm_DisgJK(state, :) = {w_ini mean_DisgJK_ini(:,:,state)' cov_ini(:, :, :, state)};
	EMISGmm_DisgKL(state, :) = {w_ini mean_DisgKL_ini(:,:,state)' cov_ini(:, :, :, state)};
	EMISGmm_SurpDC(state, :) = {w_ini mean_SurpDC_ini(:,:,state)' cov_ini(:, :, :, state)};
	EMISGmm_SurpJE(state, :) = {w_ini mean_SurpJE_ini(:,:,state)' cov_ini(:, :, :, state)};
	EMISGmm_SurpJK(state, :) = {w_ini mean_SurpJK_ini(:,:,state)' cov_ini(:, :, :, state)};
	EMISGmm_SurpKL(state, :) = {w_ini mean_SurpKL_ini(:,:,state)' cov_ini(:, :, :, state)};
end



%% HMM-GMM Training
[TRANS_NeutDC_trained, EMISGmm_NeutDC_trained, logliks_NeutDC] = GmmHmmtrain({[Train_NeutDC{:}]}, TRANS, EMISGmm_NeutDC);
[TRANS_NeutJE_trained, EMISGmm_NeutJE_trained, logliks_NeutJE] = GmmHmmtrain({[Train_NeutJE{:}]}, TRANS, EMISGmm_NeutJE);
[TRANS_NeutJK_trained, EMISGmm_NeutJK_trained, logliks_NeutJK] = GmmHmmtrain({[Train_NeutJK{:}]}, TRANS, EMISGmm_NeutJK);
[TRANS_NeutKL_trained, EMISGmm_NeutKL_trained, logliks_NeutKL] = GmmHmmtrain({[Train_NeutKL{:}]}, TRANS, EMISGmm_NeutKL);
[TRANS_HpyDC_trained, EMISGmm_HpyDC_trained, logliks_HpyDC] = GmmHmmtrain({[Train_HpyDC{:}]}, TRANS, EMISGmm_HpyDC);
[TRANS_HpyJE_trained, EMISGmm_HpyJE_trained, logliks_HpyJE] = GmmHmmtrain({[Train_HpyJE{:}]}, TRANS, EMISGmm_HpyJE);
[TRANS_HpyJK_trained, EMISGmm_HpyJK_trained, logliks_HpyJK] = GmmHmmtrain({[Train_HpyJK{:}]}, TRANS, EMISGmm_HpyJK);
[TRANS_HpyKL_trained, EMISGmm_HpyKL_trained, logliks_HpyKL] = GmmHmmtrain({[Train_HpyKL{:}]}, TRANS, EMISGmm_HpyKL);
[TRANS_SadDC_trained, EMISGmm_SadDC_trained, logliks_SadDC] = GmmHmmtrain({[Train_SadDC{:}]}, TRANS, EMISGmm_SadDC);
[TRANS_SadJE_trained, EMISGmm_SadJE_trained, logliks_SadJE] = GmmHmmtrain({[Train_SadJE{:}]}, TRANS, EMISGmm_SadJE);
[TRANS_SadJK_trained, EMISGmm_SadJK_trained, logliks_SadJK] = GmmHmmtrain({[Train_SadJK{:}]}, TRANS, EMISGmm_SadJK);
[TRANS_SadKL_trained, EMISGmm_SadKL_trained, logliks_SadKL] = GmmHmmtrain({[Train_SadKL{:}]}, TRANS, EMISGmm_SadKL);
[TRANS_AgyDC_trained, EMISGmm_AgyDC_trained, logliks_AgyDC] = GmmHmmtrain({[Train_AgyDC{:}]}, TRANS, EMISGmm_AgyDC);
[TRANS_AgyJE_trained, EMISGmm_AgyJE_trained, logliks_AgyJE] = GmmHmmtrain({[Train_AgyJE{:}]}, TRANS, EMISGmm_AgyJE);
[TRANS_AgyJK_trained, EMISGmm_AgyJK_trained, logliks_AgyJK] = GmmHmmtrain({[Train_AgyJK{:}]}, TRANS, EMISGmm_AgyJK);
[TRANS_AgyKL_trained, EMISGmm_AgyKL_trained, logliks_AgyKL] = GmmHmmtrain({[Train_AgyKL{:}]}, TRANS, EMISGmm_AgyKL);
[TRANS_FearDC_trained, EMISGmm_FearDC_trained, logliks_FearDC] = GmmHmmtrain({[Train_FearDC{:}]}, TRANS, EMISGmm_FearDC);
[TRANS_FearJE_trained, EMISGmm_FearJE_trained, logliks_FearJE] = GmmHmmtrain({[Train_FearJE{:}]}, TRANS, EMISGmm_FearJE);
[TRANS_FearJK_trained, EMISGmm_FearJK_trained, logliks_FearJK] = GmmHmmtrain({[Train_FearJK{:}]}, TRANS, EMISGmm_FearJK);
[TRANS_FearKL_trained, EMISGmm_FearKL_trained, logliks_FearKL] = GmmHmmtrain({[Train_FearKL{:}]}, TRANS, EMISGmm_FearKL);
[TRANS_DisgDC_trained, EMISGmm_DisgDC_trained, logliks_DisgDC] = GmmHmmtrain({[Train_DisgDC{:}]}, TRANS, EMISGmm_DisgDC);
[TRANS_DisgJE_trained, EMISGmm_DisgJE_trained, logliks_DisgJE] = GmmHmmtrain({[Train_DisgJE{:}]}, TRANS, EMISGmm_DisgJE);
[TRANS_DisgJK_trained, EMISGmm_DisgJK_trained, logliks_DisgJK] = GmmHmmtrain({[Train_DisgJK{:}]}, TRANS, EMISGmm_DisgJK);
[TRANS_DisgKL_trained, EMISGmm_DisgKL_trained, logliks_DisgKL] = GmmHmmtrain({[Train_DisgKL{:}]}, TRANS, EMISGmm_DisgKL);
% [TRANS_SurpDC_trained, EMISGmm_SurpDC_trained, logliks_SurpDC] = GmmHmmtrain({[Train_SurpDC{:}]}, TRANS, EMISGmm_SurpDC);
% [TRANS_SurpJE_trained, EMISGmm_SurpJE_trained, logliks_SurpJE] = GmmHmmtrain({[Train_SurpJE{:}]}, TRANS, EMISGmm_SurpJE);
% [TRANS_SurpJK_trained, EMISGmm_SurpJK_trained, logliks_SurpJK] = GmmHmmtrain({[Train_SurpJK{:}]}, TRANS, EMISGmm_SurpJK);
% [TRANS_SurpKL_trained, EMISGmm_SurpKL_trained, logliks_SurpKL] = GmmHmmtrain({[Train_SurpKL{:}]}, TRANS, EMISGmm_SurpKL);




%% Testing
% % for speaker DC
model_NeutDC = {TRANS_NeutDC_trained EMISGmm_NeutDC_trained};
model_HpyDC = {TRANS_HpyDC_trained EMISGmm_HpyDC_trained};
model_SadDC = {TRANS_SadDC_trained EMISGmm_SadDC_trained};
model_AgyDC = {TRANS_AgyDC_trained EMISGmm_AgyDC_trained};
model_FearDC = {TRANS_FearDC_trained EMISGmm_FearDC_trained};
model_DisgDC = {TRANS_DisgDC_trained EMISGmm_DisgDC_trained};
% model_SurpDC = {TRANS_SurpDC_trained EMISGmm_SurpDC_trained};
% models_DC =  {model_NeutDC model_HpyDC model_SadDC model_AgyDC model_FearDC model_DisgDC model_SurpDC};
models_DC =  {model_NeutDC model_HpyDC model_SadDC model_AgyDC model_FearDC model_DisgDC};

[recog_NeutDC, rate_NeutDC, confu_NeutDC] = GmmHmmrecog(Test_NeutDC, models_DC, 1);
[recog_HpyDC, rate_HpyDC, confu_HpyDC] = GmmHmmrecog(Test_HpyDC, models_DC, 2);
[recog_SadDC, rate_SadDC, confu_SadDC] = GmmHmmrecog(Test_SadDC, models_DC, 3);
[recog_AgyDC, rate_AgyDC, confu_AgyDC] = GmmHmmrecog(Test_AgyDC, models_DC, 4);
[recog_FearDC, rate_FearDC, confu_FearDC] = GmmHmmrecog(Test_FearDC, models_DC, 5);
[recog_DisgDC, rate_DisgDC, confu_DisgDC] = GmmHmmrecog(Test_DisgDC, models_DC, 6);


% % for speaker JE
model_NeutJE = {TRANS_NeutJE_trained EMISGmm_NeutJE_trained};
model_HpyJE = {TRANS_HpyJE_trained EMISGmm_HpyJE_trained};
model_SadJE = {TRANS_SadJE_trained EMISGmm_SadJE_trained};
model_AgyJE = {TRANS_AgyJE_trained EMISGmm_AgyJE_trained};
model_FearJE = {TRANS_FearJE_trained EMISGmm_FearJE_trained};
model_DisgJE = {TRANS_DisgJE_trained EMISGmm_DisgJE_trained};
% model_SurpJE = {TRANS_SurpJE_trained EMISGmm_SurpJE_trained};
% % models_JE =  {model_NeutJE model_HpyJE model_SadJE model_AgyJE model_FearJE model_DisgJE model_SurpJE};
models_JE =  {model_NeutJE model_HpyJE model_SadJE model_AgyJE model_FearJE model_DisgJE};

[recog_NeutJE, rate_NeutJE, confu_NeutJE] = GmmHmmrecog(Test_NeutJE, models_JE, 1);
[recog_HpyJE, rate_HpyJE, confu_HpyJE] = GmmHmmrecog(Test_HpyJE, models_JE, 2);
[recog_SadJE, rate_SadJE, confu_SadJE] = GmmHmmrecog(Test_SadJE, models_JE, 3);
[recog_AgyJE, rate_AgyJE, confu_AgyJE] = GmmHmmrecog(Test_AgyJE, models_JE, 4);
[recog_FearJE, rate_FearJE, confu_FearJE] = GmmHmmrecog(Test_FearJE, models_JE, 5);
[recog_DisgJE, rate_DisgJE, confu_DisgJE] = GmmHmmrecog(Test_DisgJE, models_JE, 6);


% % for speaker JK
model_NeutJK = {TRANS_NeutJK_trained EMISGmm_NeutJK_trained};
model_HpyJK = {TRANS_HpyJK_trained EMISGmm_HpyJK_trained};
model_SadJK = {TRANS_SadJK_trained EMISGmm_SadJK_trained};
model_AgyJK = {TRANS_AgyJK_trained EMISGmm_AgyJK_trained};
model_FearJK = {TRANS_FearJK_trained EMISGmm_FearJK_trained};
model_DisgJK = {TRANS_DisgJK_trained EMISGmm_DisgJK_trained};
% model_SurpJK = {TRANS_SurpJK_trained EMISGmm_SurpJK_trained};
% % models_JK =  {model_NeutJK model_HpyJK model_SadJK model_AgyJK model_FearJK model_DisgJK model_SurpJK};
models_JK =  {model_NeutJK model_HpyJK model_SadJK model_AgyJK model_FearJK model_DisgJK};

[recog_NeutJK, rate_NeutJK, confu_NeutJK] = GmmHmmrecog(Test_NeutJK, models_JK, 1);
[recog_HpyJK, rate_HpyJK, confu_HpyJK] = GmmHmmrecog(Test_HpyJK, models_JK, 2);
[recog_SadJK, rate_SadJK, confu_SadJK] = GmmHmmrecog(Test_SadJK, models_JK, 3);
[recog_AgyJK, rate_AgyJK, confu_AgyJK] = GmmHmmrecog(Test_AgyJK, models_JK, 4);
[recog_FearJK, rate_FearJK, confu_FearJK] = GmmHmmrecog(Test_FearJK, models_JK, 5);
[recog_DisgJK, rate_DisgJK, confu_DisgJK] = GmmHmmrecog(Test_DisgJK, models_JK, 6);


% % for speaker KL
model_NeutKL = {TRANS_NeutKL_trained EMISGmm_NeutKL_trained};
model_HpyKL = {TRANS_HpyKL_trained EMISGmm_HpyKL_trained};
model_SadKL = {TRANS_SadKL_trained EMISGmm_SadKL_trained};
model_AgyKL = {TRANS_AgyKL_trained EMISGmm_AgyKL_trained};
model_FearKL = {TRANS_FearKL_trained EMISGmm_FearKL_trained};
model_DisgKL = {TRANS_DisgKL_trained EMISGmm_DisgKL_trained};
% model_SurpKL = {TRANS_SurpKL_trained EMISGmm_SurpKL_trained};
% % models_KL =  {model_NeutKL model_HpyKL model_SadKL model_AgyKL model_FearKL model_DisgKL model_SurpKL};
models_KL =  {model_NeutKL model_HpyKL model_SadKL model_AgyKL model_FearKL model_DisgKL};

[recog_NeutKL, rate_NeutKL, confu_NeutKL] = GmmHmmrecog(Test_NeutKL, models_KL, 1);
[recog_HpyKL, rate_HpyKL, confu_HpyKL] = GmmHmmrecog(Test_HpyKL, models_KL, 2);
[recog_SadKL, rate_SadKL, confu_SadKL] = GmmHmmrecog(Test_SadKL, models_KL, 3);
[recog_AgyKL, rate_AgyKL, confu_AgyKL] = GmmHmmrecog(Test_AgyKL, models_KL, 4);
[recog_FearKL, rate_FearKL, confu_FearKL] = GmmHmmrecog(Test_FearKL, models_KL, 5);
[recog_DisgKL, rate_DisgKL, confu_DisgKL] = GmmHmmrecog(Test_DisgKL, models_KL, 6);


%%
% Repeat experiment code 5 times to have 5 sets of data.
% e.g. DC1, DC2, DC3, DC4, DC5

DC1 = [confu_AgyDC; confu_DisgDC; confu_FearDC; confu_HpyDC; confu_NeutDC; confu_SadDC];
JE1 = [confu_AgyJE; confu_DisgJE; confu_FearJE; confu_HpyJE; confu_NeutJE; confu_SadJE];
JK1 = [confu_AgyJK; confu_DisgJK; confu_FearJK; confu_HpyJK; confu_NeutJK; confu_SadJK];
KL1 = [confu_AgyKL; confu_DisgKL; confu_FearKL; confu_HpyKL; confu_NeutKL; confu_SadKL];



%%
% Estimate the averge accuracy and confusion table.
ansDC = (DC1+DC2+DC3+DC4+DC5)/5;
ansJE = (JE1+JE2+JE3+JE4+JE5)/5;
ansJK = (JK1+JK2+JK3+JK4+JK5)/5;
ansKL = (KL1+KL2+KL3+KL4+KL5)/5;