% The script to extract audio data with speaker dependant from SAVEE
% database, and save as MAT file.

clc; clear; close all;

%% Get filename and path
NeutDC = dir('.\Signals\Neutral\SAVEE\DC');
oldpath = path;
path(oldpath,'.\Signals\Neutral\SAVEE\DC');
for i = 1:size(NeutDC, 1)-2
    [x_NeutDC{i}, fs_NeutDC(i)] = audioread(NeutDC(i+2).name);
end

NeutJE = dir('.\Signals\Neutral\SAVEE\JE');
oldpath = path;
path(oldpath,'.\Signals\Neutral\SAVEE\JE');
for i = 1:size(NeutJE, 1)-2
    [x_NeutJE{i}, fs_NeutJE(i)] = audioread(NeutJE(i+2).name);
end

NeutJK = dir('.\Signals\Neutral\SAVEE\JK');
oldpath = path;
path(oldpath,'.\Signals\Neutral\SAVEE\JK');
for i = 1:size(NeutJK, 1)-2
    [x_NeutJK{i}, fs_NeutJK(i)] = audioread(NeutJK(i+2).name);
end

NeutKL = dir('.\Signals\Neutral\SAVEE\KL');
oldpath = path;
path(oldpath,'.\Signals\Neutral\SAVEE\KL');
for i = 1:size(NeutKL, 1)-2
    [x_NeutKL{i}, fs_NeutKL(i)] = audioread(NeutKL(i+2).name);
end

HpyDC = dir('.\Signals\Happy\SAVEE\DC');
oldpath = path;
path(oldpath,'.\Signals\Happy\SAVEE\DC');
for i = 1:size(HpyDC, 1)-2
    [x_HpyDC{i}, fs_HpyDC(i)] = audioread(HpyDC(i+2).name);
end

HpyJE = dir('.\Signals\Happy\SAVEE\JE');
oldpath = path;
path(oldpath,'.\Signals\Happy\SAVEE\JE');
for i = 1:size(HpyJE, 1)-2
    [x_HpyJE{i}, fs_HpyJE(i)] = audioread(HpyJE(i+2).name);
end

HpyJK = dir('.\Signals\Happy\SAVEE\JK');
oldpath = path;
path(oldpath,'.\Signals\Happy\SAVEE\JK');
for i = 1:size(HpyJK, 1)-2
    [x_HpyJK{i}, fs_HpyJK(i)] = audioread(HpyJK(i+2).name);
end

HpyKL = dir('.\Signals\Happy\SAVEE\KL');
oldpath = path;
path(oldpath,'.\Signals\Happy\SAVEE\KL');
for i = 1:size(HpyKL, 1)-2
    [x_HpyKL{i}, fs_HpyKL(i)] = audioread(HpyKL(i+2).name);
end

SadDC = dir('.\Signals\Sad\SAVEE\DC');
oldpath = path;
path(oldpath,'.\Signals\Sad\SAVEE\DC');
for i = 1:size(SadDC, 1)-2
    [x_SadDC{i}, fs_SadDC(i)] = audioread(SadDC(i+2).name);
end

SadJE = dir('.\Signals\Sad\SAVEE\JE');
oldpath = path;
path(oldpath,'.\Signals\Sad\SAVEE\JE');
for i = 1:size(SadJE, 1)-2
    [x_SadJE{i}, fs_SadJE(i)] = audioread(SadJE(i+2).name);
end

SadJK = dir('.\Signals\Sad\SAVEE\JK');
oldpath = path;
path(oldpath,'.\Signals\Sad\SAVEE\JK');
for i = 1:size(SadJK, 1)-2
    [x_SadJK{i}, fs_SadJK(i)] = audioread(SadJK(i+2).name);
end

SadKL = dir('.\Signals\Sad\SAVEE\KL');
oldpath = path;
path(oldpath,'.\Signals\Sad\SAVEE\KL');
for i = 1:size(SadKL, 1)-2
    [x_SadKL{i}, fs_SadKL(i)] = audioread(SadKL(i+2).name);
end

AgyDC = dir('.\Signals\Angry\SAVEE\DC');
oldpath = path;
path(oldpath,'.\Signals\Angry\SAVEE\DC');
for i = 1:size(AgyDC, 1)-2
    [x_AgyDC{i}, fs_AgyDC(i)] = audioread(AgyDC(i+2).name);
end

AgyJE = dir('.\Signals\Angry\SAVEE\JE');
oldpath = path;
path(oldpath,'.\Signals\Angry\SAVEE\JE');
for i = 1:size(AgyJE, 1)-2
    [x_AgyJE{i}, fs_AgyJE(i)] = audioread(AgyJE(i+2).name);
end

AgyJK = dir('.\Signals\Angry\SAVEE\JK');
oldpath = path;
path(oldpath,'.\Signals\Angry\SAVEE\JK');
for i = 1:size(AgyJK, 1)-2
    [x_AgyJK{i}, fs_AgyJK(i)] = audioread(AgyJK(i+2).name);
end

AgyKL = dir('.\Signals\Angry\SAVEE\KL');
oldpath = path;
path(oldpath,'.\Signals\Angry\SAVEE\KL');
for i = 1:size(AgyKL, 1)-2
    [x_AgyKL{i}, fs_AgyKL(i)] = audioread(AgyKL(i+2).name);
end

FearDC = dir('.\Signals\Fearful\SAVEE\DC');
oldpath = path;
path(oldpath,'.\Signals\Fearful\SAVEE\DC');
for i = 1:size(FearDC, 1)-2
    [x_FearDC{i}, fs_FearDC(i)] = audioread(FearDC(i+2).name);
end

FearJE = dir('.\Signals\Fearful\SAVEE\JE');
oldpath = path;
path(oldpath,'.\Signals\Fearful\SAVEE\JE');
for i = 1:size(FearJE, 1)-2
    [x_FearJE{i}, fs_FearJE(i)] = audioread(FearJE(i+2).name);
end

FearJK = dir('.\Signals\Fearful\SAVEE\JK');
oldpath = path;
path(oldpath,'.\Signals\Fearful\SAVEE\JK');
for i = 1:size(FearJK, 1)-2
    [x_FearJK{i}, fs_FearJK(i)] = audioread(FearJK(i+2).name);
end

FearKL = dir('.\Signals\Fearful\SAVEE\KL');
oldpath = path;
path(oldpath,'.\Signals\Fearful\SAVEE\KL');
for i = 1:size(FearKL, 1)-2
    [x_FearKL{i}, fs_FearKL(i)] = audioread(FearKL(i+2).name);
end

DisgDC = dir('.\Signals\Disgust\SAVEE\DC');
oldpath = path;
path(oldpath,'.\Signals\Disgust\SAVEE\DC');
for i = 1:size(DisgDC, 1)-2
    [x_DisgDC{i}, fs_DisgDC(i)] = audioread(DisgDC(i+2).name);
end

DisgJE = dir('.\Signals\Disgust\SAVEE\JE');
oldpath = path;
path(oldpath,'.\Signals\Disgust\SAVEE\JE');
for i = 1:size(DisgJE, 1)-2
    [x_DisgJE{i}, fs_DisgJE(i)] = audioread(DisgJE(i+2).name);
end

DisgJK = dir('.\Signals\Disgust\SAVEE\JK');
oldpath = path;
path(oldpath,'.\Signals\Disgust\SAVEE\JK');
for i = 1:size(DisgJK, 1)-2
    [x_DisgJK{i}, fs_DisgJK(i)] = audioread(DisgJK(i+2).name);
end

DisgKL = dir('.\Signals\Disgust\SAVEE\KL');
oldpath = path;
path(oldpath,'.\Signals\Disgust\SAVEE\KL');
for i = 1:size(DisgKL, 1)-2
    [x_DisgKL{i}, fs_DisgKL(i)] = audioread(DisgKL(i+2).name);
end

SurpDC = dir('.\Signals\Surprised\SAVEE\DC');
oldpath = path;
path(oldpath,'.\Signals\Surprised\SAVEE\DC');
for i = 1:size(SurpDC, 1)-2
    [x_SurpDC{i}, fs_SurpDC(i)] = audioread(SurpDC(i+2).name);
end

SurpJE = dir('.\Signals\Surprised\SAVEE\JE');
oldpath = path;
path(oldpath,'.\Signals\Surprised\SAVEE\JE');
for i = 1:size(SurpJE, 1)-2
    [x_SurpJE{i}, fs_SurpJE(i)] = audioread(SurpJE(i+2).name);
end

SurpJK = dir('.\Signals\Surprised\SAVEE\JK');
oldpath = path;
path(oldpath,'.\Signals\Surprised\SAVEE\JK');
for i = 1:size(SurpJK, 1)-2
    [x_SurpJK{i}, fs_SurpJK(i)] = audioread(SurpJK(i+2).name);
end

SurpKL = dir('.\Signals\Surprised\SAVEE\KL');
oldpath = path;
path(oldpath,'.\Signals\Surprised\SAVEE\KL');
for i = 1:size(SurpKL, 1)-2
    [x_SurpKL{i}, fs_SurpKL(i)] = audioread(SurpKL(i+2).name);
end



%% Save as MAT File
clear oldpath i;
clear NeutDC NeutJE NeutJK NeutKL ...
      HpyDC HpyJE HpyJK HpyKL ... 
      SadDC SadJE SadJK SadKL ... 
      AgyDC AgyJE AgyJK AgyKL ... 
      FearDC FearJE FearJK FearKL ... 
      DisgDC DisgJE DisgJK DisgKL ... 
      SurpDC SurpJE SurpJK SurpKL;

cd '.\MAT file';
save("SAVEE.mat");
