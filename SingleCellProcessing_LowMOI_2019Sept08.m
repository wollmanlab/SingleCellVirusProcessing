%% This script loads the raw image stacks, creates SAR, makes drift correction, creates 2d EDoF projections,
%%makes a new MD for these, run a batch PIV analysis for the sequence and
%%do post analysis
BaseStr = regexprep([char(ispc.*'Z:\Images2019\') char(isunix.*'/bigstore/Images2019/')],char(0),'');
Usr = 'Jen';
Project = 'NFkBDynamics';
Dataset = 'TNF_Titr_HSV_MOI1_withReplic_Rep3_2019Sep03';
acquisition = 1;
%% Get MD of raw data
acqname = ['acq_' num2str(acquisition)];
fpath = [BaseStr Usr filesep Project filesep Dataset filesep acqname];
MD=Metadata(fpath,[],1);
Wells = unique(MD.getSpecificMetadata('Position'));
frames = unique(cell2mat(MD.getSpecificMetadata('frame')));

cTNF = ([100*(1/2).^[0:8], 0]);

%% Drift correction
for j=1:numel(Wells)
    MD.CalculateDriftCorrection(Wells{j});
end
MD.saveMetadataMat;
MD=Metadata(fpath,[],1);

%% Make results object

R = MultiPositionSingleCellVirusResults(fpath)
R.PosNames=unique(MD.getSpecificMetadata('Position'));

frames = unique(cell2mat(MD.getSpecificMetadata('frame')));
R.Frames = frames;
R.analysisScript=fullfile([fpath filesep 'SingleCellProcessing_LowMOI_2019Sept08.m']);%Change this to the right file
R.reportPth = [BaseStr 'Reports' filesep 'Alon' filesep Project filesep Dataset];

    R.setData('cTNF',cTNF)

% define nuclear channel and channel to use for tracking 
NucChannel = 'DeepBlue'; %Channel to segment on
TrackChannel = 'Cyan'; %"smoothness" channel, usually virus

%% Add WellsLbl to results, to single cell segmentation

for WellNum=1:numel(Wells)
    pos = Wells{WellNum};
    Welllbl = cell(numel(frames),1);
    
    parfor i=R.Frames'
        Welllbl{i} = WellsConstructor_v3(fpath, pos,i,NucChannel,'EstSize',5,'maskflag', false, 'hthresh',0.001);
    end
    R.setWellsLbl(Welllbl,pos)
end
%R.saveResults


%% Correct vir channel because of weird Cyan thing
VirusChannel = 'Cyan';

for i=1:numel(R.PosNames);
    WellCells = R.getWellsLbl(R.PosNames{i});
    WellCellsCorr = WellCells;
    indVirusChnl = find(strcmp(VirusChannel,WellCells{1}.channels));

    for j=1:numel(WellCells)
        Virus90Intensities = WellCells{j}.Int90Prctile{indVirusChnl}./prctile(WellCells{j}.Int90Prctile{indVirusChnl},1);
        WellCellsCorr{j}.Int90Prctile{indVirusChnl} = Virus90Intensities;
        VirusIntensities = WellCells{j}.Intensities{indVirusChnl}./prctile(WellCells{j}.Intensities{indVirusChnl},1);
        WellCellsCorr{j}.Intensities{indVirusChnl} = VirusIntensities;
    end
    
    
    R.setWellsLbl(WellCellsCorr, R.PosNames{i});

end
R.saveResults

 %% Link adjecent frames
for WellNum=1:numel(R.PosNames)  
    R.Link(R.PosNames(WellNum),TrackChannel, 'intMultiplier', 50)
    WellNum
end

%% Close gaps
for WellNum=1:numel(R.PosNames)
    R.closeGaps(R.PosNames(WellNum),TrackChannel,NucChannel,'maxAmpDiff', 1.5)
end
%%

R.saveResults





%% Load Results
fpath = '/bigstore/Images2019/Jen/NFkBDynamics/TNF_DRUGS_HSVMOI1_2019Jun20/acq_1';
R = MultiPositionSingleCellVirusResults(fpath)
Wells = R.PosNames;
frames = R.Frames;

%% Calculate more features about single cell tracks
DeathChannel = 'FarRed';
VirusChannel = 'Cyan';
NucChannel = 'DeepBlue'; 

for WellNum=1:numel(Wells)
    Tracks = R.getTracks(R.PosNames{WellNum});
    
    %Time
    T = arrayfun(@(x) (x.seqOfEvents(1,1):x.seqOfEvents(2,1)), Tracks,'UniformOutput',false);
    [Tracks.('T')] = T{:};
    
    
    WellCells = R.getWellsLbl(R.PosNames{WellNum});
    
    indDeathChnl = find(strcmp(DeathChannel,WellCells{1}.channels));
    indVirusChnl = find(strcmp(VirusChannel,WellCells{1}.channels));
    indNucChnl = find(strcmp(NucChannel,WellCells{1}.channels));
    
    %correct Vir channel
    VirusIntensities = cellfun(@(x) x.Int90Prctile{indVirusChnl}./prctile(x.Int90Prctile{indVirusChnl},0.1), WellCells,'UniformOutput', false);
    
    NuclearIntensities = cellfun(@(x) x.Int90Prctile{indNucChnl}, WellCells,'UniformOutput', false);
    DeathIntensities = cellfun(@(x) x.Int90Prctile{indDeathChnl}, WellCells,'UniformOutput', false);
    
    for i=1:numel(Tracks)
        i
        VirusTrack{i} = zeros(1,numel(Tracks(i).T));
        DeathTrack{i} = zeros(1,numel(Tracks(i).T));
        NucTrack{i} = zeros(1,numel(Tracks(i).T));
        
        for j=1:numel(Tracks(i).T)
            % j
            if Tracks(i).tracksFeatIndxCG(j)
                VirusTrack{i}(j) = VirusIntensities{Tracks(i).T(j)}(Tracks(i).tracksFeatIndxCG(j))';
                NucTrack{i}(j) = NuclearIntensities{Tracks(i).T(j)}(Tracks(i).tracksFeatIndxCG(j))';
                DeathTrack{i}(j) = DeathIntensities{Tracks(i).T(j)}(Tracks(i).tracksFeatIndxCG(j))';

            else
                VirusTrack{i}(j) = NaN;
                NucTrack{i}(j) = NaN;
                DeathTrack{i}(j) = NaN;
            end
        end
    end
    [Tracks.('VirusTrack')] = VirusTrack{:};
    [Tracks.('DeathTrack')] = DeathTrack{:};
    [Tracks.('NucTrack')] = NucTrack{:};
    
    
    %WellCells = R.getWellsLbl(R.PosNames{35});
    %ThreshVirus = mean(WellCells{50}.Int90Prctile{indVirusChnl})+std(WellCells{50}.Int90Prctile{indVirusChnl});
    %ThreshVirus = 2;%triThreshImg(WellCells{100}.Ints90('Cyan'));
    A = cellfun(@(y) cellfun(@(x) prctile(x.Int90Prctile{indVirusChnl}, 99.99), y),R.WellLbls(1:30),'UniformOutput', false);
    ThreshVirus = max(cat(2,A{:}),[],2)';
    
    %WellCells = R.getWellsLbl(R.PosNames{WellNum});    
    %ThreshDeath = mean(WellCells{1}.Int90Prctile{indDeathChnl})+std(WellCells{1}.Int90Prctile{indDeathChnl});
    %ThreshDeath = triThreshImg(WellCells{100}.Ints90('FarRed')); 
    B = cellfun(@(y) cellfun(@(x) prctile(x.Int90Prctile{indDeathChnl}, 95), y),R.WellLbls(10:10:30),'UniformOutput', false);
    ThreshDeath = max(cat(2,B{:}),[],2)';
    
    for i=1:numel(Tracks)
        Tracks(i).Infected = Smoothing(Tracks(i).VirusTrack)>ThreshVirus(Tracks(i).T);
        %Tracks(i).CellsGetInfected = sum(Tracks(i).Infected)>2;
        Tracks(i).indWhenCellGetsInfected = find(conv(double(Tracks(i).Infected),[1 1 1 1 1], 'same')>=4,1,'first')-1;
        
        if ~isempty(Tracks(i).indWhenCellGetsInfected)
            Tracks(i).Infected(1:Tracks(i).indWhenCellGetsInfected-1)=0;
            Tracks(i).Infected(Tracks(i).indWhenCellGetsInfected:end)=1;
        else
            Tracks(i).Infected(:)=0;
        end
        Tracks(i).CellsGetInfected = any(Tracks(i).Infected);% && sum(diff(Tracks(i).Dead))==1;

        
        Tracks(i).Dead = (single(Smoothing(Tracks(i).DeathTrack,'neigh', 5)>ThreshDeath(Tracks(i).T))  + single(cumsum(Tracks(i).NucTrack>0.9)>10))>0;
        Tracks(i).DeadSum = cumsum(Tracks(i).Dead)>3;
        Tracks(i).indWhenCellDies = find(conv(double(Tracks(i).Dead),[1 1 1 1 1], 'same')>=4,1,'first')-1;
        if ~isempty(Tracks(i).indWhenCellDies)
            Tracks(i).Dead(1:Tracks(i).indWhenCellDies-1)=0;
            Tracks(i).Dead(Tracks(i).indWhenCellDies:end)=1;
        else
            Tracks(i).Dead(:)=0;
        end
        Tracks(i).CellDies = any(Tracks(i).Dead);% && sum(diff(Tracks(i).Dead))==1;

        
    end
    
    R.setTracks(Tracks,R.PosNames{WellNum})
end
clearvars DeathTrack  VirusTrack Tracks;
%R.saveResults



















%% measure live cells over time
W1 = R.getWellsLbl(R.PosNames{j});
i=40
scatter(log(W1{i}.Ints90('DeepBlue')),log(W1{i}.Ints90('FarRed')));pause(0.1)
rect = getrect;
xv = [rect(1) rect(1) rect(1)+rect(3) rect(1)+rect(3)]
yv = [rect(2) rect(2)+rect(4) rect(2)+rect(4) rect(2)]
for j=1:numel(R.PosNames)
    pos = R.PosNames{j}
    W1 = R.getWellsLbl(pos);
%     i=40
%     scatter(log(W1{i}.Ints90('DeepBlue')),log(W1{i}.Ints90('FarRed')));pause(0.1)
%     set(gca,'XTick',log([0.01 0.05 0.1 0.25 0.5 0.75 1]), 'XTickLabel',[0.01 0.05 0.1 0.25 0.5 0.75 1]...
%         ,'XLim',log([0.06 0.7]), 'YLim',log([0.008 0.05]))
%     rect = getrect;
%     xv = [rect(1) rect(1) rect(1)+rect(3) rect(1)+rect(3)]
%     yv = [rect(2) rect(2)+rect(4) rect(2)+rect(4) rect(2)]
    
    fracLive = [];
    numLive = [];
    for i=1:145;
        cla
        xq = log(W1{i}.Ints90('DeepBlue'));
        yq = log(W1{i}.Ints90('FarRed'));
        in = inpolygon(xq,yq,xv, yv);
        %plot(xq(in),yq(in),'bo') % points inside
        %hold on
        %plot(xq(~in),yq(~in),'r+') % points outside
        %pause(0.01);
        shg;
        fracLive = [fracLive sum(in)/W1{i}.num];
        numLive = [numLive, sum(in)];
    end
    
    plot(fracLive); shg
    pause(0.1)
    R.setData('fracLive',fracLive,pos)
    R.setData('numLive',numLive,pos)
end



%% measure confusion matrix over time
j=40
W1 = R.getWellsLbl(R.PosNames{j});
i=200
[Cent, Lbls] = quadrantPlot(log(W1{i}.Ints90('Cyan')),single(log(W1{i}.Ints90('DeepBlue'))));

for j=1:numel(R.PosNames)
    pos = R.PosNames{j}
    W1 = R.getWellsLbl(pos);
    
    
    TN = [];
    FN = [];
    TP = [];
    FP = [];
    for i=1:numel(W1);
        cla
        x1 = log(W1{i}.Ints90('Cyan'));
        x2 = single(log(W1{i}.Ints90('DeepBlue')));
        [~, ~, ~,QuadPer] = quadrantPlot(x1,x2, 'crosshair', Cent,'show',false);
        TN = [TN QuadPer(4)];
        FN = [FN QuadPer(2)];
        TP = [TP QuadPer(1)];
        FP = [FP QuadPer(3)];
    end
    
    
    R.setData('TN',TN,pos)
    R.setData('FN',FN,pos)
    R.setData('TP',TP,pos)
    R.setData('FP',FP,pos)
    
end

%%
figure('color','w','Position',[100,100, 450, 450])
ax2 = axes('Position', [0.15 0.15 0.7 0.7])

LinearRange = 0.4

tp = 48*3+1
hold on
TP = R.getData('TP');
TP = cat(1,TP{:});


FP = R.getData('FP');
FP = cat(1,FP{:});


TN = R.getData('TN');
TN = cat(1,TN{:});


FN = R.getData('FN');
FN = cat(1,FN{:});

totalCells = TP+FP+FN+TN;
tc = reshape(totalCells(:,tp),6,[])'

TP = reshape(TP(:,tp),6,[])';
h1 = ploterr(asinh(cTNF/LinearRange),mean(TP(:,4:6)')./mean(mean(tc(:,4:6)')),[],std(TP(:,4:6)')./mean(mean(tc(:,4:6)')))
for i=1:length(h1) h1(i).Color = deadinfColor; end

FP = reshape(FP(:,tp),6,[])';
h2 = ploterr(asinh(cTNF/LinearRange),mean(FP(:,4:6)')./mean(mean(tc(:,4:6)')),[],std(FP(:,4:6)')./mean(mean(tc(:,4:6)')))
for i=1:length(h2) h2(i).Color = deathColor; end

TN = reshape(TN(:,tp),6,[])';
h3 = ploterr(asinh(cTNF/LinearRange),mean(TN(:,4:6)')./mean(mean(tc(:,4:6)')),[],std(TN(:,4:6)')./mean(mean(tc(:,4:6)')))
for i=1:length(h3) h3(i).Color = liveColor; end

FN = reshape(FN(:,tp),6,[])';
h4 = ploterr(asinh(cTNF/LinearRange),mean(FN(:,4:6)')./mean(mean(tc(:,4:6)')),[],std(FN(:,4:6)')./mean(mean(tc(:,4:6)')))
for i=1:length(h4) h4(i).Color = virusColor; end

set(gca,'YLim',[-1,100]/100,'YTick',[0 1])

a1 = sort(- (1:9)'*10.^(1:2));
a2 =(1:9)'*10.^(1:5);
a = [sort(a1(:)); 0; a2(:)]/1000;
emptyStr = {''; ''; '';''; ''; '';''; ''};
ax2.XTickLabel = [emptyStr; ''; emptyStr; {''}; '0'; {''}; emptyStr; '  .1'; emptyStr; ' 1'; emptyStr; ' 10'; emptyStr; '    100'; emptyStr];
ax2.XTick = asinh(a/LinearRange);

ax2.XLim=asinh([-0.1 100]/LinearRange)
ax2.YLim=[-1 100]/100
ax2.YTick = 0:0.2:1;
ax2.YTickLabel = 100*[0:0.2:1];
ax2.Box = 'on';

xlabel('[TNF] ng/ml')
ylabel('Percent of cells')

hleg = legend([h3(2) h4(2) h1(2) h2(2)],'   TN - Healthy','   FN - Infected','   TP - Dead infected','   FP - Dead uninfected');
hleg.Box='off';
hleg.Position(1) = hleg.Position(1);