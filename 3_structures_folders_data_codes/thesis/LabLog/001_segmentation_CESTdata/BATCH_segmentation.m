
% load codes that you need for the batch file (first clone from  https://github.com/cest-sources/CEST_EVAL)
addpath(genpath('C:\thesis\Codes\CEST_EVAL_github\'));  % adds all folders and subfolders to the matlab path for this session.

% load data that you need for the batch file ( you cna find Example_data.mat in the above git)
load('C:\thesis\Data\2022_09_29_CEST_example_data\Example_data.mat')


% now your local codes can be played out

%% lets display some data from the data we loaded,and calculate a binary mask
A=Mz_stack(:,:,1,1);
Segment =A>200;

figure,
subplot(1,2,1), imagesc(A);
subplot(1,2,2), imagesc(Segment);
%%% lets use a function (montage1t, NORM_ZSTACK) from the loaded package
Z=NORM_ZSTACK(Mz_stack,M0_stack,P,Segment); % normalizes Mz_stack using M0_stack, and applies Segment
figure, 
subplot(2,2,1),montage1t(Mz_stack); title('raw data');
subplot(2,2,3),montage1t(Z); title('segmented and normalized data');

% lets plot the 4th dimensio of one pixel
subplot(2,2,2),plot(P.SEQ.w,squeeze( mean(mean(Mz_stack(:,:,1,:)))) );title('mean of raw data');
subplot(2,2,4),plot(P.SEQ.w,squeeze( mean(mean(Z(:,:,1,:)))) ); title('mean of segmented and normalized data');




