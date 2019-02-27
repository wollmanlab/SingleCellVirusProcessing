%% This script loads the raw image stacks, creates SAR, makes drift correction, creates 2d EDoF projections,
%%makes a new MD for these, run a batch PIV analysis for the sequence and
%%do post analysis
BaseStr = regexprep([char(ispc.*'Z:\Images2018\') char(isunix.*'/bigstore/Images2018/')],char(0),'');
Usr = 'Jen';
Project = 'NFkBDynamics';
Dataset = 'H2BiRFP_ICRP_TestTL_20x_2018Nov30';
acquisition = 5;

% define nuclear channel and channel to use for tracking
NucChannel = 'FarRed'; %Channel to segment on
TrackChannel = 'Yellow'; %"smoothness" channel, usually virus
DeathChannel = 'FarRed';
%% Get MD of raw data
acqname = ['acq_' num2str(acquisition)];
fpath = [BaseStr Usr filesep Project filesep Dataset filesep acqname];
MD=Metadata(fpath);
Wells = unique(MD.getSpecificMetadata('Position'));
frames = unique(cell2mat(MD.getSpecificMetadata('frame')));

%% Drift correction, here, we'll calculate the drift in a few of the planes and just take the mean.
%If we really want, we could calculate a distinct drift for each frame
%specifically. That's probably unnecessary and computationally expensive.

% Channels = {'DeepBlue','Yellow','Red'};
% Zindexes = max(unique(cell2mat(MD.getSpecificMetadata('Zindex'))));
% ZsToLoad = 1;
%
% for j=1:numel(Wells)
%     driftXY.dX = 0;
%     driftXY.dY = 0;
%     for i=1:numel(frames)-1;
%         Data = stkread(MD,'Channel',Channels{1}, 'flatfieldcorrection', false, 'frame', [i i+1], 'Position', Wells{j}, 'Zindex', ZsToLoad);
%         imSize = size(Data);
%         Data = reshape(Data,[imSize./[1,1,2],2]);
%         imSize = size(Data);
%
%         dX = [];
%         dY = [];
%
%         for i=1:imSize(3)
%             currRef = Data(:,:,i,1);
%             currCorr = Data(:,:,i,2);
%             %if you don't value your time, feel free to use imregtform or
%             %whatever else you see fit
%             imXcorr = convnfft(currRef - mean(currRef(:)),rot90(currCorr,2)-mean(currCorr(:)),'same');
%             [maxCorrX, maxCorrY] = find(imXcorr == max(imXcorr(:)));
%             dX = [dX maxCorrX-size(currRef,1)/2];
%             dY = [dY maxCorrY-size(currRef,2)/2];
%         end
%         driftXY.dX = [driftXY.dX round(mean(dX))];
%         driftXY.dY = [driftXY.dY round(mean(dY))]
%     end
%     CummulDriftXY.dX = cumsum(driftXY.dX);
%     CummulDriftXY.dY = cumsum(driftXY.dY);
%
%
%
%     %% Add drift to MD
%     Typ = MD.Types;
%     Vals = MD.Values;
%
%     if ~any(strcmp('driftTform',Typ))
%         Typ{end+1}='driftTform'; %Will become a standard in MD.
%     end
%     Ntypes = size(Typ,2);
%     % put the right drift displacements in the right place
%     for i=1:numel(frames)
%         i
%         inds = MD.getIndex({'frame', 'Position'},{i, Wells{j}});
%         for j1=1:numel(inds)
%             Vals{inds(j1),Ntypes} = [1 0 0 , 0 1 0 , CummulDriftXY.dY(i), CummulDriftXY.dX(i), 1];
%         end
%     end
%     MD.Types = Typ;
%     MD.Values = Vals;
% end
% MD.saveMetadataMat;
% files = dir([MD.pth filesep 'Metadata.txt']);
% if ~isempty(files)
%     movefile([MD.pth filesep files(1).name], [MD.pth filesep 'Metadata_BAK.txt'],'f');
% end



% %% New way to make movies
% j=7
% Channels = {'DeepBlue', 'Red'};
% Data = stkread(MD,'sortby','Channel','Channel',Channels, 'flatfieldcorrection', false,'blindflatfield',true, 'Position', Wells{j},'resize', 0.5,'register',true);
%
% stkshow(Data);
% PixelSize = MD.getSpecificMetadataByIndex('PixelSize', 1);
% PixelSize = PixelSize{1};
%
% MIJ.selectWindow('Data');
% MIJ.run('Stack to Hyperstack...', ['order=xytcz channels=' num2str(numel(Channels)) ' slices=1 frames=' num2str(size(Data,3)/numel(Channels)) ' display=Composite']); %Make into xytc hyperstack
% MIJ.run('Properties...', ['channels=' num2str(numel(Channels)) ' slices=1 frames=' num2str(size(Data,3)/numel(Channels)) ' unit=um pixel_width=' num2str(PixelSize) ' pixel_height=' num2str(PixelSize) ' voxel_depth=1.0000']);%set pixel size
% MIJ.run('Scale Bar...', 'width=100 height=8 font=28 color=White background=None location=[Lower Right] bold overlay');%add scale bar
% MIJ.run('Time Stamper', ['starting=0 interval=0.5 x=0 y=50 font=50 decimal=1 anti-aliased or=h overlay']);%add time stamp
%
% f = figure;%wait for user to adjust colors, intensity, etc
% set(f,'Name','Please adjust brightness and contrast','Position',[360 538 240 60],'NumberTitle', 'off', 'Toolbar','none','Menubar','none')
% h = uicontrol('Position',[20 20 200 40],'String','OK',...
%               'Callback','uiresume(gcbf)');
% uiwait(f);
% close(f);
% %% save
% inds = MD.getIndex({'Position'}, {Wells{j}});
% savepath = MD.getImageFilename({'index'}, {inds(1)});
% sepInd = strfind(savepath, filesep);
% savepath = savepath(1:sepInd(end));
% MIJ.run('AVI... ', ['compression=JPEG frame=14 save=' sprintf('%sMovie_%s%s',savepath,Wells{j},'.avi')]);
%
%

%% Make results object

R = MultiPositionSingleCellVirusResults(fpath)
R.PosNames=unique(MD.getSpecificMetadata('Position'));

frames = unique(cell2mat(MD.getSpecificMetadata('frame')));
R.Frames = frames;
R.analysisScript=fullfile([fpath filesep 'SingleCellAnalysis_2018Dec04.m']);%Change this to the right file
R.reportPth = [BaseStr 'Reports' filesep 'Alon' filesep Project filesep Dataset];



% %% Calculate flat fields
% pos = Wells{1};
% Channels = MD.unique('Channel');
% FF = struct;
% for i=1:numel(Channels)
%     img = stkread(MD,'Channel',Channels{i}, 'flatfieldcorrection', false,'blindflatfield',false, 'frame', 1, 'Position', pos,'register',false);
%     FFimg = awt2Dlite(img,8);
%     FF(i).img = squeeze(FFimg(:,:,:,end));
%     FF(i).channel = Channels{i};
% end
%% Add WellsLbl to results, to single cell segmentation

for WellNum=1:numel(Wells)
    pos = Wells{WellNum};
    Welllbl = cell(numel(frames)-1,1);
    
    parfor i=R.Frames(1:end-1)'
        Welllbl{i} = WellsConstructor_v2FRET(fpath, pos,i,NucChannel);
    end
    
    
    R.setWellsLbl(Welllbl,pos)
end
%%
R.saveResults;
%% Link adjecent frames
for WellNum=1:numel(R.PosNames)
    WellNum
    R.Link(R.PosNames(WellNum),TrackChannel)
end

%% Close gaps
for WellNum=1:numel(R.PosNames)
    R.closeGaps(R.PosNames(WellNum),TrackChannel,DeathChannel)
end

R.saveResults

%% Calculate more features about single cell tracks

CyanChannel = 'Cyan';
YellowChannel = 'Yellow';
DeadC = [];
%CyanC = [];
%Ratio = [];
for WellNum=1:numel(R.PosNames)
    
    WellCells = R.getWellsLbl(R.PosNames{WellNum});
    
    %indtrckChnl = find(strcmp(TrackChannel,WellCells{1}.channels));
    indCyanChnl = find(strcmp(CyanChannel,WellCells{1}.channels));
    indYellowChnl = find(strcmp(YellowChannel,WellCells{1}.channels));
    indNucChnl = find(strcmp(NucChannel,WellCells{1}.channels));

    
    CyanIntensities = cellfun(@(x) x.Intensities{indCyanChnl}, WellCells,'UniformOutput', false);
    YellowIntensities = cellfun(@(x) x.Intensities{indYellowChnl}, WellCells,'UniformOutput', false);
    RatioInts = cellfun(@(x) x.Intensities{indCyanChnl}./x.Intensities{indYellowChnl}, WellCells,'UniformOutput', false);
    NucIntensities = cellfun(@(x) x.Intensities{indNucChnl}, WellCells,'UniformOutput', false);

    
    Tracks = R.getTracks(R.PosNames{WellNum});

    if numel(Tracks)
        %Time
        T = arrayfun(@(x) (x.seqOfEvents(1,1):x.seqOfEvents(2,1)), Tracks,'UniformOutput',false);
        [Tracks.('T')] = T{:};
        for i=1:numel(Tracks)
            i
            CyanTrack{i} = zeros(1,numel(Tracks(i).T));
            YellowTrack{i} = zeros(1,numel(Tracks(i).T));
            RatioTrack{i} = zeros(1,numel(Tracks(i).T));
            NucTrack{i} = zeros(1,numel(Tracks(i).T));
            
            for j=1:numel(Tracks(i).T)
                % j
                if Tracks(i).tracksFeatIndxCG(j)
                    CyanTrack{i}(j) = CyanIntensities{Tracks(i).T(j)}(Tracks(i).tracksFeatIndxCG(j))';
                    YellowTrack{i}(j) = YellowIntensities{Tracks(i).T(j)}(Tracks(i).tracksFeatIndxCG(j))';
                    RatioTrack{i}(j) = RatioInts{Tracks(i).T(j)}(Tracks(i).tracksFeatIndxCG(j))';
                    NucTrack{i}(j) = NucIntensities{Tracks(i).T(j)}(Tracks(i).tracksFeatIndxCG(j))';
                else
                    CyanTrack{i}(j) = NaN;
                    YellowTrack{i}(j) = NaN;
                    RatioTrack{i}(j) = NaN;
                    NucTrack{i}(j) = NaN;
                end
            end
        end
        [Tracks.('CyanTrack')] = CyanTrack{:};
        [Tracks.('YellowTrack')] = YellowTrack{:};
        [Tracks.('RatioTrack')] = RatioTrack{:};
        [Tracks.('NucTrack')] = NucTrack{:};
        
        %ThreshVirus = mean(WellCells{1}.Int90Prctile{indVirusChnl})+2*std(WellCells{1}.Int90Prctile{indVirusChnl});
        ThreshDeath = nanmean(WellCells{1}.Int90Prctile{indNucChnl})-1*nanstd(WellCells{1}.Int90Prctile{indNucChnl});
        
    %   for i=1:numel(Tracks)
            %         Tracks(i).Infected = Smoothing(Tracks(i).VirusTrack)>ThreshVirus;
            %         Tracks(i).CellsGetInfected = sum(Tracks(i).Infected)>4;
    %        Tracks(i).Dead = Smoothing(Tracks(i).NucTrack,'neigh', 5)<ThreshDeath;
    %        Tracks(i).CellDies = sum(Tracks(i).Dead)>4 && sum(diff(Tracks(i).Dead))==1;
    %    end
        %
        R.setTracks(Tracks,R.PosNames{WellNum})
    end
    
    %DeadC = [DeadC sum(NucIntensities{:}<ThreshDeath)./numel(YellowIntensities{:})];
    %InfectedC = [InfectedC sum(VirusIntensities{:}>ThreshVirus)./numel(VirusIntensities{:})];
    
    
    
    
end
clearvars DeathTrack  VirusTrack Tracks RatioTrack;
%%
R.saveResults

%%

DeadC = reshape(DeadC, [4 48]);
InfectedC = reshape(InfectedC, [4 48]);

DeadCmn = mean(DeadC); DeadCst = std(DeadC)./sqrt(4);
InfectedCmn = mean(InfectedC); InfectedCst = std(InfectedC)./sqrt(4);

%%

TNF = [100 33.3 11.1 3.7 1.23 0.411 0.137 0.1];

figure()
subplot(2,2,1)
hold on
errorbar(TNF, DeadCmn(1:6:end)*100, DeadCst(1:6:end)*100, '-ko', 'markerfacecolor','k','markersize',9, 'linewidth', 2)
errorbar(TNF, DeadCmn(3:6:end)*100, DeadCst(3:6:end)*100, '-ko', 'markerfacecolor', 'r', 'markersize', 9, 'linewidth', 2)
set(gca, 'fontsize', 13, 'xlim', [0.08 120], 'xscale', 'log')
xlabel('[TNF] (ng/ml)')
ylabel('Dead Cells (%)')
legend('p65V MEF','WT MEF')
title('24h, no virus')
box on
subplot(2,2,2)
hold on
errorbar(TNF, DeadCmn(2:6:end)*100, DeadCst(2:6:end)*100, '-ko', 'markerfacecolor','k','markersize',9, 'linewidth', 2)
errorbar(TNF, DeadCmn(4:6:end)*100, DeadCst(4:6:end)*100, '-ko', 'markerfacecolor', 'r', 'markersize', 9, 'linewidth', 2)
set(gca, 'fontsize', 13, 'xlim', [0.08 120], 'xscale', 'log')
xlabel('[TNF] (ng/ml)')
ylabel('Dead Cells (%)')
%legend('p65V MEF','WT MEF')
title('24h + HSV')
box on
subplot(2,2,3)
hold on
errorbar(TNF, InfectedCmn(1:6:end)*100, InfectedCst(1:6:end)*100, '-ko', 'markerfacecolor','k','markersize',9, 'linewidth', 2)
errorbar(TNF, InfectedCmn(3:6:end)*100, InfectedCst(3:6:end)*100, '-ko', 'markerfacecolor', 'r', 'markersize', 9, 'linewidth', 2)
set(gca, 'fontsize', 13, 'xlim', [0.08 120], 'xscale', 'log')
xlabel('[TNF] (ng/ml)')
ylabel('HSV+ Cells (%)')
%legend('p65V MEF','WT MEF')
title('24h, no virus')
box on
subplot(2,2,4)
hold on
errorbar(TNF, InfectedCmn(2:6:end)*100, InfectedCst(2:6:end)*100, '-ko', 'markerfacecolor','k','markersize',9, 'linewidth', 2)
errorbar(TNF, InfectedCmn(4:6:end)*100, InfectedCst(4:6:end)*100, '-ko', 'markerfacecolor', 'r', 'markersize', 9, 'linewidth', 2)
set(gca, 'fontsize', 13, 'xlim', [0.08 120], 'xscale', 'log')
xlabel('[TNF] (ng/ml)')
ylabel('HSV+ Cells (%)')
%legend('p65V MEF','WT MEF')
title('24h + HSV')
box on

%% Load results object

pth = '/bigstore/Images2018/Jen/NFkBDynamics/ApoptosisModification_2018Oct18/acq_2';
R = MultiPositionSingleCellVirusResults(pth);

%%
%dtinfdeathCell = {};
% nErrDeadCells = {};
% for WellNum=1:numel(R.PosNames)
% %Tracks = R.getTracks(R.PosNames{WellNum});
%
% Jinfdie = intersect(R.TracksThatDie(R.PosNames{WellNum}),R.TracksThatGetInfected(R.PosNames{WellNum}));
% Jerrdie = intersect(R.TracksThatDie(R.PosNames{WellNum}),find(~[Tracks.CellsGetInfected]));
%
% %dtinfdeath = NaN(numel(Jinfdie),1);
% for i=1:numel(Jinfdie)
% dtinfdeath(i) = find(diff([0 Tracks(Jinfdie(i)).Dead])==1,1,'last')-find([Tracks(Jinfdie(i)).Infected],1);
% end
% numel(R.TracksThatGetInfected(R.PosNames{WellNum}));
% mean(dtinfdeath)
% dtinfdeathCell{WellNum} = dtinfdeath;
% nErrDeadCells{WellNum} = numel(Jerrdie);


end

%% Plot for jen
close all
figure('Position',[82 318 1024 280],'color','w')
%subplot(2,4,1)
h = histogram(dtinfdeathCell{30},10,'Normalization','pdf','FaceColor','k')
dim = [.2 .60 .3 .3];
hbx = annotation('textbox',dim,'String',['mean = ' num2str(mean(dtinfdeathCell{30})/2,2) 'h'],'LineStyle','none')
title('No TNF')
set(gca,'xtick',[-100:50:100], 'xticklabels',[-100:50:100]/2)
ylabel('pdf')
xlabel('Time between infection and death(h)')
hbx = annotation('arrow',[0.1, 0.4],[0.5,0.5])
%%
edges = h.BinEdges;
subplot(2,4,2)

h7 = histogram(dtinfdeathCell{29},edges,'Normalization','pdf','FaceColor','k')
dim = [.41 .60 .3 .3];
hbx = annotation('textbox',dim,'String',['median = ' num2str(mean(dtinfdeathCell{29})/2,2) 'h'],'LineStyle','none')
title('0.14ng/ml TNF')
set(gca,'xtick',[-100:50:100], 'xticklabels',[-100:50:100]/2)

subplot(2,4,3)

h6 = histogram(dtinfdeathCell{10},edges,'Normalization','pdf','FaceColor','k')
dim = [.62 .60 .3 .3];
hbx = annotation('textbox',dim,'String',['median = ' num2str(mean(dtinfdeathCell{10})/2,2)  'h'],'LineStyle','none')
title('0.411ng/ml TNF')
set(gca,'xtick',[-100:50:100], 'xticklabels',[-100:50:100]/2)

subplot(2,4,4)

h5 = histogram(dtinfdeathCell{2},edges,'Normalization','pdf','FaceColor','k')
dim = [.83 .60 .3 .3];
hbx = annotation('textbox',dim,'String',['median = ' num2str(mean(dtinfdeathCell{2})/2,2) 'h'],'LineStyle','none')
title('1.23ng/ml TNF')
set(gca,'xtick',[-100:50:100], 'xticklabels',[-100:50:100]/2)

subplot(2,4,5)

h4 = histogram(dtinfdeathCell{25},edges,'Normalization','pdf','FaceColor','k')
dim = [.2 .13 .3 .3];
hbx = annotation('textbox',dim,'String',['median = ' num2str(mean(dtinfdeathCell{25})/2,2) 'h'],'LineStyle','none')
title('3.7ng/ml TNF')
set(gca,'xtick',[-100:50:100], 'xticklabels',[-100:50:100]/2)

subplot(2,4,6)

h3 = histogram(dtinfdeathCell{17},edges,'Normalization','pdf','FaceColor','k')
dim = [.41 .13 .3 .3];
hbx = annotation('textbox',dim,'String',['median = ' num2str(mean(dtinfdeathCell{17})/2,2)  'h'],'LineStyle','none')
title('11.1ng/ml TNF')
set(gca,'xtick',[-100:50:100], 'xticklabels',[-100:50:100]/2)

subplot(2,4,7)

h2 = histogram(dtinfdeathCell{9},edges,'Normalization','pdf','FaceColor','k')
dim = [.62 .13 .3 .3];
hbx = annotation('textbox',dim,'String',['median = ' num2str(mean(dtinfdeathCell{9})/2,2) 'h'],'LineStyle','none')
title('33.3ng/ml TNF')
set(gca,'xtick',[-100:50:100], 'xticklabels',[-100:50:100]/2)

subplot(2,4,8)

h1 = histogram(dtinfdeathCell{1},edges,'Normalization','pdf','FaceColor','k')
dim = [.83 .13 .3 .3];
hbx = annotation('textbox',dim,'String',['median = ' num2str(mean(dtinfdeathCell{1})/2,2) 'h'],'LineStyle','none')
title('100ng/ml TNF')
set(gca,'xtick',[-100:50:100], 'xticklabels',[-100:50:100]/2)
%%
shg
set(gcf, 'PaperPositionMode','auto','InvertHardCopy','off')
print(gcf,'-dpng','-r300',['/bigstore/Images2018/Jen/NFkBDynamics/HSV-1SpreadwTNF_2018Aug21/acq_4/' 'HistoPlot']);
%%

TNF = [100 33.3 11.1 3.7 1.23 0.411 0.137 0.1]
median = [2.5 2.5 1.8 4.5 5 6 7.5 8.5];
mean = [3.5 2.7 2.4 6.9 6.6 6.7 8.6 12];

figure()
hold on
plot (TNF, median, '-ko' ,'markerfacecolor', 'k', 'markersize', 9, 'linewidth',2)
plot(TNF, mean, '-ko', 'markerfacecolor', 'r' , 'markersize', 9, 'linewidth', 2)
set(gca, 'fontsize', 13, 'xscale', 'log')
xlabel('[TNF] (ng/ml)')
ylabel('\Deltat between infection and death (h)')
legend('Median', 'Mean')
box on
title('Effect of TNF on time interval between inf. and death')

%%  Plot error
TNF = [100 33.3 11.1 3.7 1.23 0.411 0.137 0.1];

Frac = reshape([nErrDeadCells{1:24}]./nums(1:24),3,[])


figure()
hold on
plot(TNF, reshape([nErrDeadCells{1:24}]./nums(1:24),[],3), '-o', 'linewidth', 2, 'markersize', 9)
set(gca, 'fontsize', 13, 'xscale', 'log')
xlabel('[TNF] (ng/ml)')
ylabel('Fraction dead, Error')
legend('DMSO', 'LCL-161', 'Z-Vad, Necrostatin-1')
title('Effect of different drugs on death')
box on


%%
Error = [0.3820 0.4411 0.2256 0.1295 0.0737 0.0820 0.0694 0.1082];
median = [2.5 2.5 1.8 4.5 5 6 7.5 8.5];
mean = [3.5 2.7 2.4 6.9 6.6 6.7 8.6 12];

figure()
subplot(1,2,1)
plot(median, Error, '-ko', 'markerfacecolor', 'k', 'markersize', 9)
set(gca, 'fontsize', 13)
ylabel('Error (fraction dead)')
xlabel('(\Deltat between inf and death)')
title('Median \Deltat')
subplot(1,2,2)
plot(mean, Error, '-ko', 'markerfacecolor', 'k', 'markersize', 9)
set(gca, 'fontsize', 13)
ylabel('Error (fraction dead)')
title('Mean \Deltat')
xlabel('(\Deltat between inf and death)')
%%


% %%
% j=1
% pos = R.PosNames{j};
% WellLbl = R.getWellsLbl(pos)
% fpath='/bigstore/Images2018/Jen/NFkBDynamics/HSV-1SpreadwTNF_2018Aug21/acq_4/';
% outputVideo = VideoWriter(sprintf('%sSingleCellScatter_%s%s',fpath,pos,'.avi'));
% outputVideo.FrameRate = 14;
% open(outputVideo)
%
% for i=1:181
% cla
% WellLbl{i}.scatter;
% freezeColors
% hold on
% WellLbl{i}.scatter('channel','nuclei')
% drawnow; shg
% pause(0.01)
%     frame = getframe(gcf);
%     im = frame2im(frame);
%
%     writeVideo(outputVideo,im);
% end
% close(outputVideo)
% %% single tracks
% j=6
% pos = R.PosNames{j};
%
% indinf = R.TracksThatGetInfected(pos)
% i=9
% R.plotCompTrack(pos,indinf(i))
% %%
% set(gcf, 'PaperPositionMode','auto','InvertHardCopy','off')
% print(gcf,'-dpng','-r300',['/bigstore/Images2018/Jen/NFkBDynamics/HSV-1SpreadwTNF_2018Aug21/acq_4/' 'SingleTrack' pos '_track' num2str(indinf(i))]);