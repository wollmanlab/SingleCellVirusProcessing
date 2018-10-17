%% This script loads the raw image stacks, creates SAR, makes drift correction, creates 2d EDoF projections,
%%makes a new MD for these, run a batch PIV analysis for the sequence and
%%do post analysis
BaseStr = regexprep([char(ispc.*'Z:\Images2018\') char(isunix.*'/bigstore/Images2018/')],char(0),'');
Usr = 'Jen';
Project = 'NFkBDynamics';
Dataset = 'HSV-1SpreadwTNF_2018Aug21';
acquisition = 4;
%% Get MD of raw data
acqname = ['acq_' num2str(acquisition)];
fpath = [BaseStr Usr filesep Project filesep Dataset filesep acqname];
MD=Metadata(fpath);
Wells = unique(MD.getSpecificMetadata('Position'));
frames = unique(cell2mat(MD.getSpecificMetadata('frame')));

%% Drift correction, here, we'll calculate the drift in a few of the planes and just take the mean.
%If we really want, we could calculate a distinct drift for each frame
%specifically. That's probably unnecessary and computationally expensive.

Channels = {'DeepBlue', 'Red', 'Yellow', 'FarRed'};
Zindexes = max(unique(cell2mat(MD.getSpecificMetadata('Zindex'))));
ZsToLoad = 1;

for j=1:numel(Wells)
    driftXY.dX = 0;
    driftXY.dY = 0;
    for i=1:numel(frames)-1;
        Data = stkread(MD,'Channel',Channels{1}, 'flatfieldcorrection', false, 'frame', [i i+1], 'Position', Wells{j}, 'Zindex', ZsToLoad);
        imSize = size(Data);
        Data = reshape(Data,[imSize./[1,1,2],2]);
        imSize = size(Data);
        
        dX = [];
        dY = [];
        
        for i=1:imSize(3)
            currRef = Data(:,:,i,1);
            currCorr = Data(:,:,i,2);
            %if you don't value your time, feel free to use imregtform or
            %whatever else you see fit
            imXcorr = convnfft(currRef - mean(currRef(:)),rot90(currCorr,2)-mean(currCorr(:)),'same');
            [maxCorrX, maxCorrY] = find(imXcorr == max(imXcorr(:)));
            dX = [dX maxCorrX-size(currRef,1)/2];
            dY = [dY maxCorrY-size(currRef,2)/2];
        end
        driftXY.dX = [driftXY.dX round(mean(dX))];
        driftXY.dY = [driftXY.dY round(mean(dY))]
    end
    CummulDriftXY.dX = cumsum(driftXY.dX);
    CummulDriftXY.dY = cumsum(driftXY.dY);
    
    
    
    %% Add drift to MD
    Typ = MD.Types;
    Vals = MD.Values;
    
    if ~any(strcmp('driftTform',Typ))
        Typ{end+1}='driftTform'; %Will become a standard in MD.
    end
    Ntypes = size(Typ,2);
    % put the right drift displacements in the right place
    for i=1:numel(frames)
        i
        inds = MD.getIndex({'frame', 'Position'},{i, Wells{j}});
        for j1=1:numel(inds)
            Vals{inds(j1),Ntypes} = [1 0 0 , 0 1 0 , CummulDriftXY.dY(i), CummulDriftXY.dX(i), 1];
        end
    end
    MD.Types = Typ;
    MD.Values = Vals;
end
MD.saveMetadataMat;
files = dir([MD.pth filesep 'Metadata.txt']);
if ~isempty(files)
    movefile([MD.pth filesep files(1).name], [MD.pth filesep 'Metadata_BAK.txt']);
end
MD=Metadata(fpath);



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
R.analysisScript=fullfile([fpath filesep 'AnalysisScriptTemplate.m']);%Change this to the right file
R.reportPth = [BaseStr 'Reports' filesep 'Alon' filesep Project filesep Dataset];


% define nuclear channel and channel to use for tracking
NucChannel = 'DeepBlue'; %Channel to segment on
TrackChannel = 'Red'; %"smoothness" channel, usually virus

%% Calculate flat fields
pos = Wells{1};
Channels = MD.unique('Channel');
FF = struct;
for i=1:numel(Channels)
    img = stkread(MD,'Channel',Channels{i}, 'flatfieldcorrection', false,'blindflatfield',false, 'frame', 1, 'Position', pos,'register',false);
    FFimg = awt2Dlite(img,8);
    FF(i).img = squeeze(FFimg(:,:,:,end));
    FF(i).channel = Channels{i};
end
%% Add WellsLbl to results, to single cell segmentation

for WellNum=1:numel(Wells)
    pos = Wells{WellNum};
    Welllbl = cell(numel(frames),1);
    
    parfor i=R.Frames'
        Welllbl{i} = WellsConstructor(fpath, pos,i,FF,NucChannel);
    end
    R.setWellsLbl(Welllbl,pos)
end

%% Link adjecent frames
for WellNum=1:numel(R.PosNames)
    R.Link(R.PosNames(WellNum),TrackChannel)
end

%% Close gaps
for WellNum=1:numel(R.PosNames)
    R.closeGaps(R.PosNames(WellNum),TrackChannel,NucChannel)
end

%R.saveResults
%% Load Results
%R = MultiPositionSingleCellVirusResults(fpath)
%Wells = R.PosNames;
%frames = R.Frames;


%% Calculate more features about single cell tracks


for WellNum=1:numel(Wells)
    Tracks = R.getTracks(R.PosNames{WellNum});
    
    %Time
    T = arrayfun(@(x) (x.seqOfEvents(1,1):x.seqOfEvents(2,1)), Tracks,'UniformOutput',false);
    [Tracks.('T')] = T{:};
    
    
    WellCells = R.getWellsLbl(R.PosNames{WellNum});
    
    indtrckChnl = find(strcmp(TrackChannel,WellCells{1}.channels));
    indNucChnl = find(strcmp(NucChannel,WellCells{1}.channels));
    
    
    VirusIntensities = cellfun(@(x) x.Int90Prctile{indtrckChnl}, WellCells,'UniformOutput', false);
    NuclearIntensities = cellfun(@(x) x.Int90Prctile{indNucChnl}, WellCells,'UniformOutput', false);
    
    for i=1:numel(Tracks)
        i
        VirusTrack{i} = zeros(1,numel(Tracks(i).T));
        NuclearTrack{i} = zeros(1,numel(Tracks(i).T));
        
        for j=1:numel(Tracks(i).T)
            % j
            if Tracks(i).tracksFeatIndxCG(j)
                VirusTrack{i}(j) = VirusIntensities{Tracks(i).T(j)}(Tracks(i).tracksFeatIndxCG(j))';
                NuclearTrack{i}(j) = NuclearIntensities{Tracks(i).T(j)}(Tracks(i).tracksFeatIndxCG(j))';
            else
                VirusTrack{i}(j) = NaN;
                NuclearTrack{i}(j) = NaN;
            end
        end
    end
    [Tracks.('VirusTrack')] = VirusTrack{:};
    [Tracks.('NuclearTrack')] = NuclearTrack{:};
    
    
    WellCells = R.getWellsLbl(R.PosNames{WellNum});
    ThreshVirus = mean(WellCells{1}.Int90Prctile{indtrckChnl})+2*std(WellCells{1}.Int90Prctile{indtrckChnl});
    ThreshNuc = mean(WellCells{1}.Int90Prctile{indNucChnl})+2*std(WellCells{1}.Int90Prctile{indNucChnl});
    
    for i=1:numel(Tracks)
        Tracks(i).Infected = Smoothing(Tracks(i).VirusTrack)>ThreshVirus;
        Tracks(i).CellsGetInfected = sum(Tracks(i).Infected)>4;
        Tracks(i).Dead = Smoothing(Tracks(i).NuclearTrack,'neigh', 5)>ThreshNuc;
        Tracks(i).CellDies = sum(Tracks(i).Dead)>8  && any(diff(Tracks(i).Dead)==1);% && sum(diff(Tracks(i).Dead))==1;
    end
    
    R.setTracks(Tracks,R.PosNames{WellNum})
end
clearvars NuclearTrack  VirusTrack Tracks;
R.saveResults

%
% %%
% dtinfdeathCell = {};
% for WellNum=1:numel(Wells)
% Tracks = R.getTracks(R.PosNames{WellNum});
%
% Jinfdie = intersect(R.TracksThatDie(R.PosNames{WellNum}),R.TracksThatGetInfected(R.PosNames{WellNum}));
%
% dtinfdeath = NaN(numel(Jinfdie),1);
% for i=1:numel(Jinfdie)
% dtinfdeath(i) = find(diff(Tracks(Jinfdie(i)).Dead)==1,1,'last')-find(Tracks(Jinfdie(i)).Infected,1);
% end
% numel(R.TracksThatGetInfected(R.PosNames{WellNum}));
% mean(dtinfdeath)
% dtinfdeathCell{WellNum} = dtinfdeath;
% end
%
% %% Plot for jen
% close all
% figure('Position',[82 318 1024 280],'color','w')
% subplot(1,4,1)
% h = histogram(dtinfdeathCell{12},10,'Normalization','pdf','FaceColor','k')
% dim = [.15 .63 .3 .3];
% hbx = annotation('textbox',dim,'String',['mean = ' num2str(mean(dtinfdeathCell{12})/2,2) 'h'],'LineStyle','none')
% hbx = annotation('textbox',dim+[0.02 0.05 0 0],'String',['No TNF'],'LineStyle','none')
% set(gca,'xtick',[-100:50:100], 'xticklabels',[-100:50:100]/2)
% ylabel('pdf')
% xlabel('Time between infection and death(h)')
% hbx = annotation('arrow',[0.1, 0.4],[0.05,0.05])
%
% edges = h.BinEdges;
% subplot(1,4,2)
%
% h4 = histogram(dtinfdeathCell{6},edges,'Normalization','pdf','FaceColor','c')
% dim = [.35 .63 .3 .3];
% hbx = annotation('textbox',dim,'String',['mean = ' num2str(mean(dtinfdeathCell{6})/2,2) 'h'],'LineStyle','none')
% hbx = annotation('textbox',dim+[0.02 0.05 0 0],'String',['0.3 ng/ml'],'LineStyle','none')
% set(gca,'xtick',[-100:50:100], 'xticklabels',[-100:50:100]/2)
%
% subplot(1,4,3)
%
% h3 = histogram(dtinfdeathCell{4},edges,'Normalization','pdf','FaceColor','r')
% dim = [.56 .63 .3 .3];
% hbx = annotation('textbox',dim,'String',['mean = ' num2str(mean(dtinfdeathCell{4})/2,2)  'h'],'LineStyle','none')
% hbx = annotation('textbox',dim+[0.02 0.05 0 0],'String',['3 ng/ml'],'LineStyle','none')
% set(gca,'xtick',[-100:50:100], 'xticklabels',[-100:50:100]/2)
%
% subplot(1,4,4)
%
% h2 = histogram(dtinfdeathCell{2},edges,'Normalization','pdf','FaceColor','b')
% dim = [.77 .63 .3 .3];
% hbx = annotation('textbox',dim,'String',['mean = ' num2str(mean(dtinfdeathCell{2})/2,2) 'h'],'LineStyle','none')
% hbx = annotation('textbox',dim+[0.02 0.05 0 0],'String',['30 ng/ml'],'LineStyle','none')
% set(gca,'xtick',[-100:50:100], 'xticklabels',[-100:50:100]/2)
%
% shg
% set(gcf, 'PaperPositionMode','auto','InvertHardCopy','off')
% print(gcf,'-dpng','-r300',['/bigstore/Images2018/Jen/NFkBDynamics/HSV-1SpreadwTNF_2018Aug21/acq_4/' 'HistoPlot']);
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