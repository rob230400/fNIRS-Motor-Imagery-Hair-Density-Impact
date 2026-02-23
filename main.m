%% Group 21

%% Import data and inizialization
clear all, close all, clc

warning('off')

% Loading data
n_sub = 29; % number of subjects
n_chan = 36; % number of channels
n_sess = 3; % number of sessions
n_rep = 20; % number of repetitions
MI_cnt = cell(n_sub,1); MI_mnt = MI_cnt; MI_mrk = MI_cnt;
for sub=1:n_sub
    MI_cnt{sub,1} = load(sprintf("Data/S%d/MI_cnt.mat", sub)).MI_cnt;
    MI_mnt{sub,1} = load(sprintf("Data/S%d/MI_mnt.mat", sub)).MNT;
    MI_mrk{sub,1} = load(sprintf("Data/S%d/MI_mrk.mat", sub)).MI_mrk;
end
fsamp = 10; % equal for every subject
fNy = fsamp/2; % Nyquist frequency  

clear sub

%% Import and calculation of subjects characteristics
subjects_characteristics = readtable("./Data/subjects_characteristics.xlsx");
% Age
ages = 2016-subjects_characteristics.YearOfBirth;
mean_age = mean(ages); std_age = std(ages);

% Gender
gender = subjects_characteristics.Gender;
n_males = length(find(strcmp(gender,'M')==1)); n_females = length(find(strcmp(gender,'F')==1));

% Handedness
handedness = subjects_characteristics.Handedness;
Lhandedness_sub = find(strcmp(handedness,"L")==1);

% Removing left handed subject
if length(Lhandedness_sub) < 3
    for i=Lhandedness_sub
        MI_cnt(i,:) = [];
        MI_mnt(i,:) = [];
        MI_mrk(i,:) = [];
    end
    n_sub = n_sub - length(Lhandedness_sub); % number of subjects with removal of left handedness subjects
end

% Hair Density
HairDensity = subjects_characteristics.HairDensity;
[HairDensityCount, HairDensityLabel] = histcounts(categorical(HairDensity), 'Categories', unique(HairDensity));
HairDensityLabel([1, 2, 3]) = cellfun(@(x) x, HairDensityLabel([2, 3, 1]), 'UniformOutput', false);
HairDensityLabelMean = cell(1,2); HairDensityLabelMean{1,1} = 'L-M'; HairDensityLabelMean{1,2} = 'H-VH';
HairDensityCount(1,[1 2 3]) = HairDensityCount(1,[2 3 1]);  % Forse non serve
for sub=1:n_sub
    if (strcmp(HairDensity{sub,1},"L")==1 | strcmp(HairDensity{sub,1},"M")==1)
        HairDensity{sub,1} = 'L-M';
    else
        HairDensity{sub,1} = 'H-VH';
    end
end

clear i sub ages mean_age std_age gender n_males n_females handedness Lhandedness_sub HairDensityLabel

%% Calculation of NIRS signal through modified Beer-Lambert law

% Differential path length
B = 5.97; % differential path length factor of the brain tissue
% Absorption coefficient in cm^-1/M
% (values taken from https://omlc.org/spectra/hemoglobin/summary.html)
alpha_ox = [586, 1058];
alpha_h = [1548.52, 691.32];
% Optical distance between emitter and receiver (given by the track)
d = 3; % 3cm (30mm)

% Calculation: the solution of the Beer-Lambert law is [dc_h; dc_ox] = inv(M1)*M2
% Use parfor if available
parfor sub=1:n_sub
    for s=1:n_sess
        for c=1:n_chan
            for t=1:length(MI_cnt{sub,1}{1,s}.x(:,c))
                b = [MI_cnt{sub,1}{1,s}.x(t,c); MI_cnt{sub,1}{1,s}.x(t,c+n_chan)];
                A = [alpha_ox; alpha_h]'*d*B;
                x = 1000*(A\b); % *1000 allows to obtain the signal in mM
                MI_cnt{sub,1}{1,s}.ox(t,c) = x(1);
                MI_cnt{sub,1}{1,s}.h(t,c) = x(2);
            end
        end
    end
end

clear c s sub b A x B d alpha_ox alpha_h

disp("1:  Load DONE")

%% CAR - Spatial Filtering
for sub=1:n_sub
    for s=1:n_sess
        % Calculation of CAR filter
        ox_CAR = mean(MI_cnt{sub,1}{1,s}.ox,2);
        h_CAR = mean(MI_cnt{sub,1}{1,s}.h,2);

        % Application of CAR filter
        MI_cnt{sub,1}{1,s}.ox = MI_cnt{sub,1}{1,s}.ox - ox_CAR;
        MI_cnt{sub,1}{1,s}.h = MI_cnt{sub,1}{1,s}.h - h_CAR;
    end
end

clear s sub ox_CAR h_CAR

disp("2: CAR Filtering DONE")

%% Band-Filter (sixth-order zero-phase Butterworth filter with passband of 0.01–0.1 Hz)
n = 4; Wn = 0.1/fNy;
[bl,al] = butter(n,Wn,'low');
%freqz(bl,al,2^10,fsamp)

n = 4; Wn = 0.01/fNy;
[bh,ah] = butter(n,Wn,'high');
%freqz(bh,ah,2^10,fsamp)

for sub=1:n_sub
    for s=1:n_sess
            MI_cnt{sub,1}{1,s}.ox = filtfilt(bl,al,MI_cnt{sub,1}{1,s}.ox);
            MI_cnt{sub,1}{1,s}.h = filtfilt(bl,al,MI_cnt{sub,1}{1,s}.h);
            MI_cnt{sub,1}{1,s}.ox = filtfilt(bh,ah,MI_cnt{sub,1}{1,s}.ox);
            MI_cnt{sub,1}{1,s}.h = filtfilt(bh,ah,MI_cnt{sub,1}{1,s}.h);
    end
end

clear s sub bl al n Wn bh ah

disp("3:  Band-Pass Filtering DONE")

%% Epoch division and side division

for sub=1:n_sub
    for s=1:n_sess
        onset_event = floor(MI_mrk{sub,1}{1,s}.time/1000*fsamp); % values in samples
        side_event = MI_mrk{sub, 1}{1, s}.event.desc';
        Lcont = 1;
        Rcont = 1;
        for e=1:n_rep
            if side_event(e)==1
                MI_cnt{sub,1}{1,s}.oxL(:,:,Lcont) = MI_cnt{sub,1}{1,s}.ox(onset_event(e)-5*fsamp:onset_event(e)+25*fsamp,:);
                MI_cnt{sub,1}{1,s}.hL(:,:,Lcont) = MI_cnt{sub,1}{1,s}.h(onset_event(e)-5*fsamp:onset_event(e)+25*fsamp,:);
                Lcont = Lcont+1;
            else
                MI_cnt{sub,1}{1,s}.oxR(:,:,Rcont) = MI_cnt{sub,1}{1,s}.ox(onset_event(e)-5*fsamp:onset_event(e)+25*fsamp,:);
                MI_cnt{sub,1}{1,s}.hR(:,:,Rcont) = MI_cnt{sub,1}{1,s}.h(onset_event(e)-5*fsamp:onset_event(e)+25*fsamp,:);
                Rcont = Rcont+1;
            end
        end
    end
end


clear s sub e onset_event side_event Lcont Rcont

disp("4:  Epoch and side division DONE")

%% Baseline correction

% Removal of the mean of the signal between -5s and -2s for every repetition
for sub=1:n_sub
    for s=1:n_sess
        % Calculation of the baseline
        oxL_baseline(1,:,:) = mean(MI_cnt{sub,1}{1,s}.oxL(1:3*fsamp-1,:,:),1);
        oxR_baseline(1,:,:) = mean(MI_cnt{sub,1}{1,s}.oxR(1:3*fsamp-1,:,:),1);
        hL_baseline(1,:,:) = mean(MI_cnt{sub,1}{1,s}.hL(1:3*fsamp-1,:,:),1);
        hR_baseline(1,:,:) = mean(MI_cnt{sub,1}{1,s}.hR(1:3*fsamp-1,:,:),1);

        % Removal of the baseline
        MI_cnt{sub,1}{1,s}.oxL = MI_cnt{sub,1}{1,s}.oxL - oxL_baseline;
        MI_cnt{sub,1}{1,s}.oxR = MI_cnt{sub,1}{1,s}.oxR - oxR_baseline;
        MI_cnt{sub,1}{1,s}.hL = MI_cnt{sub,1}{1,s}.hL - hL_baseline;
        MI_cnt{sub,1}{1,s}.hR = MI_cnt{sub,1}{1,s}.hR - hR_baseline;
    end
end
timevect_epoch = (0:length(MI_cnt{1,1}{1,1}.oxL)-1)/fsamp-5;    

clear s sub oxL_baseline oxR_baseline hL_baseline hR_baseline

disp("5:  Baseline correction DONE")


%% Peak distribution calculation (for t-test) (average between session) 
% Vedere se usare matrici al posto di celle
oxL_sess = cell(n_sub,1);
oxR_sess = cell(n_sub,1);
hL_sess = cell(n_sub,1);
hR_sess = cell(n_sub,1);
oxL_sessavg = cell(n_sub,1);
oxR_sessavg = cell(n_sub,1);
hL_sessavg = cell(n_sub,1);
hR_sessavg = cell(n_sub,1);

for sub=1:n_sub
    oxL_sess{sub,1}=[];
    oxR_sess{sub,1}=[];
    hL_sess{sub,1}=[];
    hR_sess{sub,1}=[];
    oxL_sessavg{sub,1} = [];
    oxR_sessavg{sub,1} = [];
    hL_sessavg{sub,1} = [];
    hR_sessavg{sub,1} = [];
    for s=1:3
        oxL_sess{sub,1} = cat(4, oxL_sess{sub,1}, MI_cnt{sub,1}{1,s}.oxL);
        oxR_sess{sub,1} = cat(4, oxR_sess{sub,1}, MI_cnt{sub,1}{1,s}.oxR);
        hL_sess{sub,1} = cat(4, hL_sess{sub,1}, MI_cnt{sub,1}{1,s}.hL);
        hR_sess{sub,1} = cat(4, hR_sess{sub,1}, MI_cnt{sub,1}{1,s}.hR);
    end
    % Average between sessions
    oxL_sessavg{sub,1} = squeeze(mean(oxL_sess{sub,1},4));
    oxR_sessavg{sub,1} = squeeze(mean(oxR_sess{sub,1},4));
    hL_sessavg{sub,1} = squeeze(mean(hL_sess{sub,1},4));
    hR_sessavg{sub,1} = squeeze(mean(hR_sess{sub,1},4));

    % Peak of signals in every epoch of the first session
    max_oxL{sub,:} = squeeze(max(oxL_sessavg{sub,1}(5*fsamp+1:20*fsamp+1,:,:),[],1))';
    max_oxR{sub,:} = squeeze(max(oxR_sessavg{sub,1}(5*fsamp+1:20*fsamp+1,:,:),[],1))'; 
    max_hL{sub,:} = squeeze(max(hL_sessavg{sub,1}(5*fsamp+1:20*fsamp+1,:,:),[],1))';
    max_hR{sub,:} = squeeze(max(hR_sessavg{sub,1}(5*fsamp+1:20*fsamp+1,:,:),[],1))'; 
end

clear sub oxR_sess oxL_sess hR_sess hL_sess oxR_sessavg oxL_sessavg hR_sessavg hL_sessavg

disp("6:  Peak distribution calculation (for t-test) DONE")

%% Average (mean) and Standard Error Computation
% Notes:
%        Average: epochs --> sessions --> subjects (Grand Average)
%        Size of 3-D Matrixes:  (sample, channel, epoch)
%        Size of 4-D Matrixes:  (sample, channel, session, subject)

% Average between epochs, sessions and channels
for sub=1:n_sub
    for s=1:3
        % Average between epochs
        oxL_meanepochs(:,:,s,sub) = mean(MI_cnt{sub,1}{1,s}.oxL,3);
        oxR_meanepochs(:,:,s,sub) = mean(MI_cnt{sub,1}{1,s}.oxR,3);
        hL_meanepochs(:,:,s,sub) = mean(MI_cnt{sub,1}{1,s}.hL,3);
        hR_meanepochs(:,:,s,sub) = mean(MI_cnt{sub,1}{1,s}.hR,3);

        % Standard error of the epochs
        oxL_stdepochs(:,:,s,sub) = std(MI_cnt{sub,1}{1,s}.oxL,0,3)/sqrt(n_rep/2);
        oxR_stdepochs(:,:,s,sub) = std(MI_cnt{sub,1}{1,s}.oxR,0,3)/sqrt(n_rep/2);
        hL_stdepochs(:,:,s,sub) = std(MI_cnt{sub,1}{1,s}.hL,0,3)/sqrt(n_rep/2);
        hR_stdepochs(:,:,s,sub) = std(MI_cnt{sub,1}{1,s}.hR,0,3)/sqrt(n_rep/2);
    end

    % Average and standard error between sessions (1, 2, 3 and all together)
    for ms = 1:4
        % Average
        if ms<4
            oxL_meansessions{ms,1} = mean(oxL_meanepochs(:,:,ms,:),3); oxR_meansessions{ms,1} = mean(oxR_meanepochs(:,:,ms,:),3);
            hL_meansessions{ms,1} = mean(hL_meanepochs(:,:,ms,:),3); hR_meansessions{ms,1} = mean(hR_meanepochs(:,:,ms,:),3);
        else
            oxL_meansessions{ms,1} = mean(oxL_meanepochs,3); oxR_meansessions{ms,1} = mean(oxR_meanepochs,3);
            hL_meansessions{ms,1} = mean(hL_meanepochs,3); hR_meansessions{ms,1} = mean(hR_meanepochs,3);
        end

        % Standard error
        if ms<4
            oxL_stdsessions{ms,1} = oxL_stdepochs(:,:,s,:); oxR_stdsessions{ms,1} = oxR_stdepochs(:,:,s,:);
            hL_stdsessions{ms,1} = hL_stdepochs(:,:,s,:); hR_stdsessions{ms,1} = hR_stdepochs(:,:,s,:);
        else
            oxL_stdsessions{ms,1} = std(oxL_meanepochs,0,3)/sqrt(n_sess); oxR_stdsessions{ms,1} = std(oxR_meanepochs,0,3)/sqrt(n_sess);
            hL_stdsessions{ms,1} = std(hL_meanepochs,0,3)/sqrt(n_sess); hR_stdsessions{ms,1} = std(hR_meanepochs,0,3)/sqrt(n_sess);
        end
    end
end

% Grand Average
% Signal
for ms = 1:4
    oxL_GrandAverageSig{ms,1} = mean(oxL_meansessions{ms,1},4); oxR_GrandAverageSig{ms,1} = mean(oxR_meansessions{ms,1},4);
    hL_GrandAverageSig{ms,1} = mean(hL_meansessions{ms,1},4); hR_GrandAverageSig{ms,1} = mean(hR_meansessions{ms,1},4);
end
% % Standard Error
for ms = 1:4
    oxL_GrandAverageStd{ms,1} = std(oxL_meansessions{ms,1},0,4)/sqrt(n_sub); oxR_GrandAverageStd{ms,1} = std(oxR_meansessions{ms,1},0,4)/sqrt(n_sub);
    hL_GrandAverageStd{ms,1} = std(hL_meansessions{ms,1},0,4)/sqrt(n_sub); hR_GrandAverageStd{ms,1} = std(hR_meansessions{ms,1},0,4)/sqrt(n_sub);
end

clear sub s ms 

disp("7:  Average and Standard Error Computation DONE")

%% Plot of the Ideal Subject
posChanSP = [2 3 5 6 8:10 12:14 16 17 19 20 22:24 26:28 30 31 33 34]; % channel position for the subplot
ordChanSP = [16 18 28 30 14 17 23 25 29 36 15 22 27 35 13 20 24 26 33 34 19 21 32 31]; % channel number associated to the subplot position
sub = 14; % chosen ideal subject (it's subject 15 of the provided dataset; here it's 14 because subject 13 was removed
sessp = 4; % Session to plot (put 1, 2, 3 or average (4))
posChanSPL = 13:24;
posChanSPR = 25:36;
actYLimL = [-2e-4 8e-4]; % y-axes limit (LMI)
actYLimR = [-4e-3 1e-3]; % y-axes limit (RMI)

titleIdSub1 = "Subject " + num2str(sub+1) + " - ";
titleIdSub2 = {"Session 1 ", "Session 2 ", "Session 3 ", "Averaged Session "};
titleIdSub3 = "- RMI";

chan_count = 1;

figure
sgtitle(titleIdSub1+titleIdSub2{sessp}+titleIdSub3, fontsize=35, FontWeight='bold')
for jj = posChanSP
    subplot(5,7,jj)
    % Mean of O2Hb and HHb
    plot(timevect_epoch, oxR_meansessions{sessp,1}(:,ordChanSP(chan_count),1,sub), Color="#D95319", LineWidth=1.5), hold on
    plot(timevect_epoch, hR_meansessions{sessp,1}(:,ordChanSP(chan_count),1,sub), Color="#0072BD", LineWidth=1.5), hold on
    % Onset and offset of the task
    xline(0, '--k', LineWidth=0.8), hold on, xline(10, '--k', LineWidth=0.8), hold on
    % Standard error bands
    ubound_oxR = oxR_meansessions{sessp,1}(:,ordChanSP(chan_count),1,sub)+oxR_stdsessions{sessp,1}(:,ordChanSP(chan_count),1,sub);
    lbound_oxR = oxR_meansessions{sessp,1}(:,ordChanSP(chan_count),1,sub)-oxR_stdsessions{sessp,1}(:,ordChanSP(chan_count),1,sub);
    ubound_hR = hR_meansessions{sessp,1}(:,ordChanSP(chan_count),1,sub)+hR_stdsessions{sessp,1}(:,ordChanSP(chan_count),1,sub);
    lbound_hR = hR_meansessions{sessp,1}(:,ordChanSP(chan_count),1,sub)-hR_stdsessions{sessp,1}(:,ordChanSP(chan_count),1,sub);
    fill([timevect_epoch, fliplr(timevect_epoch)], [ubound_oxR', fliplr(lbound_oxR')], 'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none'), hold on
    fill([timevect_epoch, fliplr(timevect_epoch)], [ubound_hR', fliplr(lbound_hR')], 'b', 'FaceAlpha', 0.1, 'EdgeColor', 'none'), hold off
    % Channel title and axis labels
    title(MI_mnt{sub,1}.clab(ordChanSP(chan_count))+" ("+num2str(ordChanSP(chan_count))+")")
    xlabel("Time (s)"), ylabel("\DeltaHb (mM)")
    % axis limit
    xlim([timevect_epoch(1) timevect_epoch(end)])
    if ismember(ordChanSP(chan_count), posChanSPL)
        ylim(actYLimL)
    else
        ylim(actYLimR)
    end

    % Update the channel counter
    chan_count = chan_count+1;
end
chan_count = chan_count-1;

% Global legend in the center of the figure
subplot(5,7,18)
plot(0, 1, Color="#D95319", LineWidth=1.5), hold on
plot(10, 1, Color="#0072BD", LineWidth=1.5), hold on
fill([0, fliplr(0)], [2', fliplr(-2')], 'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none'), hold on
fill([10, fliplr(10)], [2', fliplr(-2')], 'b', 'FaceAlpha', 0.1, 'EdgeColor', 'none'), hold on
yline(0, '--k', LineWidth=0.8), hold off
axis("off"), legend("\DeltaO_2Hb", "\DeltaHHb", "\DeltaO_2Hb Standard Error Band", "\DeltaHHb Standard Error Band", "Task Onset and Offset",  Location="best", FontSize=13.5)

disp("8:  Plot of the Ideal Subject DONE")

clear jj ubound_oxR lbound_oxR ubound_hR lbound_hR titleIdSub1 titleIdSub2 titleIdSub3 actYLimR actYLimL chan_count sessp

%% t-test: hair density (HD) - for every channel

% Inizialization of the needed cells
oxL_avHD = cell(1,2); oxR_avHD = cell(1,2);
for hd = 1:2
    oxL_avHD{1, hd} = []; oxR_avHD{1, hd} = [];
end
% Concatenation of subjects
% Every single cell of ox/h L/R _HD will be a 3D matrix
% of size: (10, 36, n° of sub with the cell-column associated HD) 
% with 36 being the n° of channels and 10 being the n° of epochs in the first session
for sub=1:n_sub
    hd = find(strcmp(HairDensityLabelMean, HairDensity{sub, 1}));
    oxL_avHD{1, hd} = cat(3, oxL_avHD{1, hd}, max_oxL{sub,1});
    oxR_avHD{1, hd} = cat(3, oxR_avHD{1, hd}, max_oxR{sub,1}); 
end
% Averaging
for hd=1:2
    oxL_avHD{1, hd} = mean(oxL_avHD{1, hd}, 3);
    oxR_avHD{1, hd} = mean(oxR_avHD{1, hd}, 3);
end


% Channel selection 
% L and R stand for left and right hemispheres
chan_act_L = [15 17 20];
chan_act_R = [29 33 35];

% Mean between channels of the same hemisphere
for hd=1:2
    oxL_avHD{1, hd} = mean(oxL_avHD{1, hd}(:,chan_act_R),2);
    oxR_avHD{1, hd} = mean(oxR_avHD{1, hd}(:,chan_act_L),2);
end

% t-test computation
% Inizialization of elements needed for the t-test computation
alpha = 0.05;

% Note: tTestResults structure:
% - struct with value of hypotesis (1 or 0), p-value and stats of t-test
% - h and p are 36x4, that is, n_chan x number of delta_hemoglobin (in order: oxL, hL, oxR, hR)
tTestResults = struct;

% t-test computation 
[hypothesisR, pvalR, ~, ~] = ttest2(oxR_avHD{1,1}, oxR_avHD{1,2}, Alpha=alpha, Vartype="equal", Tail='right');
[hypothesisL, pvalL, ~, ~] = ttest2(oxL_avHD{1,1}, oxL_avHD{1,2}, Alpha=alpha, Vartype="equal", Tail='right');
% Results storing
tTestResults.h = [hypothesisL, hypothesisR];
tTestResults.p = [pvalL, pvalR];

disp("9:  t-test: hair density (HD) for every channel DONE")

clear hd sub hypothesisL hypothesisR pvalR pvalL


%% Bar diagram of average of max_oxL/R
% Asteriks indicating significant statistical difference were displayed
% using: https://github.com/arsalanfiroozi/sigasterisk

% Elements needed to get a proper figure
screenSize = get(0, 'ScreenSize');
screenWidth = screenSize(3);
screenHeight = screenSize(4);
maxSize = min(screenWidth, screenHeight)-200;

% Data selection
Data = cat(2,oxL_avHD{1,1},oxR_avHD{1,1},oxL_avHD{1,2},oxR_avHD{1,2});

% Error-bar evaluation
error = std(Data)/sqrt(10);
error = reshape(error,2,2);

Data = mean(Data);
Data = reshape(Data,2,2);

% Bar Diagram Plot
figure
H = bar(Data); hold on
set(gca,XTickLabel={'LMI', 'RMI'}, FontSize=30, FontWeight='bold');
set(gcf, Position=[100, 100, maxSize, maxSize]);
set(H(1,1), FaceColor="#77AC30");
set(H(1,2), FaceColor="#7E2F8E");
ylabel("\DeltaO_2Hb (mM)")
legend("L-M", "H-VH", FontSize=30)
sigasterisk(1,2,1,1,"***",[2.7e-04,2.7e-04;2.7e-04,2.7e-04]), hold on
sigasterisk(1,2,2,2,"**",[2.25e-04,2.25e-04;2.25e-04,2.25e-04]), hold on
add_errorbar(error,Data), hold off

disp("10: Bar diagram of average of max_oxL/R DONE")

clear screenSize screenWidth screenHeight maxSize H

%% Topography - Evolution of Grand Average through the sessions
% Columns: 1=LMI, 2=RMI
% Rows: session 1, 2, 3 and averaged sessions (4)

% Data selection
for sess = 1:n_sess+1
    ox_MaxPeakAmp{sess,1} = max(oxL_GrandAverageSig{sess,1}(5*fsamp+1:20*fsamp+1,:));
    ox_MaxPeakAmp{sess,2} = max(oxR_GrandAverageSig{sess,1}(5*fsamp+1:20*fsamp+1,:));
end

% Plot Title definition
TotalTitle1 = sprintf("Topographic distributions of O_2Hb peak amplitudes during motor imagery task of the Grand Average\n");
TotalTitle2 = {"Session 1"; "Session 2"; "Session 3"; "Averaged Sessions"};
TotalTitle3 = sprintf("\nLMI - RMI");
colorbarTitle = sprintf("Peak Amplitude of\nO_2Hb");

% Make multiple subplots to get a bigger plot
figRows = 4;
figCols = 20;
allSubplots = 1:figRows*figCols;
subplotMatrix = reshape(allSubplots, figCols, figRows)';
subplotLMI = subplotMatrix(:, 1:figCols/2); subplotLMI = reshape(subplotLMI', [numel(subplotLMI) 1]);
subplotRMI = subplotMatrix(:, figCols/2+1:figCols); subplotRMI = reshape(subplotRMI', [numel(subplotRMI) 1]);

% Find the minimum and the maximum value of topography between all channels
% and all sessions to make comparable plots
catOx_MaxPeakAmp = cell2mat(reshape(ox_MaxPeakAmp, 1, []));
minColorbar = min(catOx_MaxPeakAmp(:));
maxColorbar = max(catOx_MaxPeakAmp(:));

% Plot
for sess = 1:n_sess
    figure
    sgtitle(TotalTitle2{sess,1}, fontsize=55, FontWeight='bold')
    

    % LMI
    subplot(figRows,figCols,subplotLMI), 
    topOx{sess,1} = fNIRS_Topography_project(ox_MaxPeakAmp{sess,1});
    set(gca, FontSize=30, FontWeight='bold')
    c = colorbar;
    c.Ruler.TickLabelFormat = '%g mM'; % display measure unit
    clim([minColorbar maxColorbar]) % set limits of colorbar
    if sess==1
        title("LMI", FontSize=40)
    end
    
    % RMI
    subplot(figRows,figCols,subplotRMI),
    topOx{sess,2} = fNIRS_Topography_project(ox_MaxPeakAmp{sess,2});
    set(gca, FontSize=20)
    c = colorbar;
    c.Ruler.TickLabelFormat = '%g mM'; % display measure unit
    clim([minColorbar maxColorbar]) % set limits of colorbar
    if sess==1
        title("RMI", FontSize=40)
    end
end

clear sess TotalTitle1 TotalTitle2 TotalTitle3 figRows figCols allSubplots subplotMatrix subplotLMI subplotRMI topOx ax c catOx_MaxPeakAmp minColorbar maxColorbar colorbarTitle

disp("11: Topography DONE")
