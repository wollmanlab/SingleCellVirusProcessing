%%
set(0,'DefaultTextInterpreter', 'tex')
set(0, 'DefaultAxesFontName', 'Arial')
set(0, 'DefaultAxesFontSize', 20)
set(0, 'DefaultUIControlFontName', 'Arial')
set(0,'defaulttextfontname','Arial');
set(0,'defaulttextfontsize',22);
set(groot,'defaultFigureColor','w')
set(groot,'defaultAxesColor','w')
set(groot,'DefaultLineMarkerSize',3)
set(groot,'defaultAxesTickLength',[0.03 0.01])
set(groot,'defaultLineLineWidth',2)


%% define colors
liveColor = [64,224,208]/255; %cyan
virusColor = [148, 0, 211]/255; %purple
deathColor = [255,165,0]/255; %orange
deadinfColor = [0.4 0.4 0.4]; % pink

cTNF = ([100*(1/2).^[0:8], 0]);

%cmTNF = @viridis;
cmTNF = @(x) makeColorMap([0,100,0]/255,[152,251,152]/255,x);

cmMAC = @gray;

%% define savepaths
savpath = '/bigstore/GeneralStorage/Alon/Figures/DecisionPaper2019/Figures020421/'

if ~isdir(savpath)
    mkdir(savpath)
    mkdir([savpath 'epss'])
    mkdir([savpath 'pngs'])
    mkdir([savpath 'movies'])
end


%% Load mac spread stuff
BaseStr = regexprep([char(ispc.*'Z:\Images2020\') char(isunix.*'/bigstore/Images2020/')],char(0),'');
Usr = 'Jen';
Project = 'NFkBDynamics';
Dataset = 'MacTitr_HSVSpread_wRepl_2020Jan14';
acquisition = 2;

% Get MD of raw data
acqname = ['acq_' num2str(acquisition)];
fpath = [BaseStr Usr filesep Project filesep Dataset filesep acqname];
MD=Metadata(fpath,[],1);
Wells = unique(MD.getSpecificMetadata('Position'));
frames = unique(cell2mat(MD.getSpecificMetadata('frame')));

% Load results object
MD=Metadata(fpath,[],1);

R = MultiPositionSingleCellVirusResults(fpath)


%% Total Frac infected at 48h
figure('color','w','Position',[100,100, 450, 450])
ax = axes('Position', [0.2, 0.15, 0.7, 0.7])

cmMAC = @gray;

tzeva = parula(10)
allTP = R.getData('TP')
allFN = R.getData('FN')
rang = [1:9 11:19];

totalInfected = (cat(1, allTP{rang})'+cat(1, allFN{rang})');
dt = 1/2;
t = 48./dt;
macs = [50*(1/2).^[0:8]];
LinearRange = 0.1;

infAtT = reshape(totalInfected(t,:),[],2)';
h1 = ploterr(asinh(macs/LinearRange),mean(infAtT),[],std(infAtT),'-o','logy');
for jj=1:numel(h1)
    h1(jj).MarkerFaceColor='k';
    h1(jj).Color='k';    
    h1(jj).MarkerSize=7;
end

%stupid matlab
BarTicklength(h1(2), 0)

shg

hold on
rang = [21:29 31:39];
totalInfected = (cat(1, allTP{rang})'+cat(1, allFN{rang})');

infAtT = reshape(totalInfected(t,:),[],2)';
h2 = ploterr(asinh(macs/LinearRange),mean(infAtT),[],std(infAtT),'-o','logy');
for jj=1:numel(h2)
    h2(jj).MarkerFaceColor='b';
    h2(jj).Color='b';    
    h2(jj).MarkerSize=7;
end
BarTicklength(h2(2), 0)
shg

hold on
rang = [41:49 51:59];
totalInfected = (cat(1, allTP{rang})'+cat(1, allFN{rang})');

infAtT = reshape(totalInfected(t,:),[],2)';
h3 = ploterr(asinh(macs/LinearRange),mean(infAtT),[],std(infAtT),'-o','logy');
for jj=1:numel(h3)
    h3(jj).MarkerFaceColor='r';
    h3(jj).Color='r';    
    h3(jj).MarkerSize=7;
end
BarTicklength(h3(2), 0)
shg



a1 = sort(- (1:9)'*10.^(1:2));
a2 =(1:9)'*10.^(1:5);
a = [sort(a1(:)); 0; a2(:)]/1000;
emptyStr = {''; ''; '';''; ''; '';''; ''};
XTickLabel = [emptyStr; ' '; emptyStr; {''}; '0'; {''}; emptyStr; ' 0.1'; emptyStr; '1'; emptyStr; '10'; emptyStr; '100'; emptyStr];
XTicks = asinh(a/LinearRange);
set(gca, 'xtick', XTicks, 'xticklabel', XTickLabel)

ax.YTick = (1/4).^([3:-1:0]);
ax.YTickLabel = num2str(100*(1/4).^([3:-1:0])','%3.0f');
ax.YLim = [0.005, 1];

ax.XLim = asinh([0.1, 100]/LinearRange)
ylabel('Infected cells (%)')
xlabel('Macrophage density (%)')

hleg = legend([h1(1),h2(1),h3(1)], 'Activated WT', 'Naive WT', 'Activated TNF KO')
hleg.Box='off'
hleg.Position = [0.33    0.2370    0.2267    0.1578]


%% Save
figname = 'Fig1_InfectedVsMacsat48h'
set(gcf, 'PaperPositionMode','auto','color','w','InvertHardcopy','off')
print(gcf,'-depsc','-r600',[savpath '/epss/' figname]);
print(gcf,'-dpng','-r600',[savpath '/pngs/' figname]);




%% Load Ki67 stuff
pth = '/bigstore/Images2019/Jen/NFkBDynamics/HSV_TNF_Ki67_2019Nov25/acq_2';
R = MultiPositionSingleCellVirusResults(pth);
%% plot with errors
fracProlif = R.getData('fracProlif');
close all

figure('color','w','Position',[100,100, 450, 450])
ax = axes('Position', [0.2, 0.15, 0.7, 0.7])

fracMean = squeeze(mean(reshape(fracProlif(1:160),4,10,4),1));
fracStd = squeeze(std(reshape(fracProlif(1:160),4,10,4),[],1));
tzeva = lines(4);

%uninfected
LinearRange = 0.4;
J = 1:2
h1 = ploterr(asinh(cTNF/LinearRange), mean(fracMean(:,J),2),[], createDistanceMatrix([0 0], fracStd(:,J)),'o');
for jj=1:numel(h1)
    h1(jj).Color = liveColor;
    h1(jj).MarkerFaceColor = liveColor;
    h1(jj).MarkerSize=7;
end

hold on

b = annotation('textbox',[0.22    0.1506    0.689    0.1078],'String', ['Uninfected doubling time=' num2str(log(2)/mean(mean(fracMean(:,J),2)),2) '\pm' num2str(sqrt(10)\log(2)*std(mean(fracMean(:,J),2))/mean(mean(fracMean(:,J),2))^2,1) 'h'])
b.EdgeColor = 'none'
b.FontSize=12;

pRateUninf = mean(mean(fracMean(:,J),2));

%infected
J = 3:4
h2 = ploterr(asinh(cTNF/LinearRange), mean(fracMean(:,J),2),[], createDistanceMatrix([0 0], fracStd(:,J)),'o');
for jj=1:numel(h2)
    h2(jj).Color = virusColor;
    h2(jj).MarkerFaceColor = virusColor;
    h2(jj).MarkerSize=7;
end
hold on
%b = annotation('textbox',[0.22    0.1506    0.689    0.1078],'String', ['Uninfected doubling time=' num2str(log(2)/mean(mean(fracMean(:,J),2)),2) '\pm' num2str(log(2)*std(mean(fracMean(:,J),2))/mean(mean(fracMean(:,J),2))^2,1) 'h'])

b = annotation('textbox',[0.22    0.1106    0.689    0.1078],'String', ['Infected doubling time=' num2str(log(2)/mean(mean(fracMean(:,J),2)),2) '\pm' num2str(sqrt(10)\log(2)*std(mean(fracMean(:,J),2))/mean(mean(fracMean(:,J),2))^2,2) 'h'])
b.EdgeColor = 'none'
b.FontSize=12

pRateInf = mean(mean(fracMean(:,J),2));


shg
ylabel('Percent dividing cells')
xlabel('[TNF\alpha] ng/ml')
a1 = sort(- (1:9)'*10.^(1:2));
a2 =(1:9)'*10.^(1:5);
a = [sort(a1(:)); 0; a2(:)]/1000;
emptyStr = {''; ''; '';''; ''; '';''; ''};
XTickLabel = [emptyStr; ''; emptyStr; {''}; '0'; {''}; emptyStr; ' 0.1'; emptyStr; '    1'; emptyStr; '   10'; emptyStr; '    100'; emptyStr];
XTicks = asinh(a/LinearRange);
set(gca, 'xtick', XTicks, 'xticklabel', XTickLabel)

ax.XLim = asinh([-0.1, 150]/LinearRange)
ax.YTick = -0.01:0.01:0.1;
ax.YTickLabel = 100*[-0.01:0.01:0.1];

hl = legend([h1(2), h2(2)],'Uninfected', 'Infected')

hl.Box = 'off'
hl.Position = [0.217    0.8506    0.3089    0.1078];

%% Save
figname = 'Fig2S_ProliferationRates'
set(gcf, 'PaperPositionMode','auto','color','w','InvertHardcopy','off')
print(gcf,'-depsc','-r600',[savpath '/epss/' figname]);
print(gcf,'-dpng','-r600',[savpath '/pngs/' figname]);




%% Load results for infectivity stuff
fpath = '/bigstore/Images2019/Jen/NFkBDynamics/Infectivity_CoCulture_2019Nov01/acq_2';
R = MultiPositionSingleCellVirusResults(fpath)

nnvl = R.getData('nnvl');
f = R.getData('probInf');

%% 
figure('color','w','Position',[100,100, 450, 450])
ax = axes('Position', [0.2, 0.15, 0.7, 0.7])

dt = 1/3;

X = cat(2,nnvl{:})';
Y = dt\cat(2,f{:})';

Y = Y(X~=0);
X = X(X~=0);
%J = find((X>prctile(X, 0.1)) .* (X<prctile(X, 99.9)));

%X = X(J);
%Y = Y(J);

[Xb, Yb, stdXb, stdYb, steXb, steYb, Xt, Yt] = BinData_v3(log(X), Y, 25);

h = ploterr(exp(Xb(1:end-1)), Yb(1:end-1), [], steYb(1:end-1),'o','logx');
BarTicklength(h(2),0)

for jj=1:numel(h)
    h(jj).Color = 'k';
    h(jj).MarkerFaceColor = 'k';
    h(jj).MarkerSize=7;
end

shg
hold on

[BETA,resnorm,residual,exitflag,output,lambda,jacobian] = lsqcurvefit(@HillFunction,[1 100 0], exp(double(Xb([1 3:end-1]))), Yb([1 3:end-1]),[0 0 0],[inf inf 0 ])
%BETA= nlinfitWeight2(exp(double(Xb(1:end-1))), Yb(1:end-1), @HillFunctionWithCoef,[1 .15 0 1], steYb(1:end-1)+0.000001, [],[0 0 0 0.5],[inf 0.3 inf 3])
%BETA = [1.5 1000, 0]
xx = 0.5:1:1000;
yy = HillFunction(BETA, xx);
plot((xx), yy,'k')
set(gca,'XScale','log', 'XLim',[1 1200],'YLim',[-0.01 0.8],'XTick',[1 10 100 1000])
xlabel('\nu - Nearest neighbors viral load (a.u.)')
ylabel('IR - Infection rate (h^{-1})')
ci = nlparci(BETA,residual,'jacobian',jacobian);

text(1.5, 0.6, '$IR(\nu) = 1.5\cdot\frac{1}{1+\frac{1000}{\nu}}$','interpreter','latex' );

y=Yb([1 3:end-1]);
yfit=HillFunction(BETA, exp(double(Xb([1 3:end-1]))));
SStot = sum((y-mean(y)).^2);                            % Total Sum-Of-Squares
SSres = sum((y(:)-yfit(:)).^2);                         % Residual Sum-Of-Squares
Rsq = 1-SSres/SStot;                                    % R^2
%% Save
figname = 'Fig3S_InfectionRateVsLoad'
set(gcf, 'PaperPositionMode','auto','color','w','InvertHardcopy','off')
print(gcf,'-depsc','-r600',[savpath '/epss/' figname]);
print(gcf,'-dpng','-r600',[savpath '/pngs/' figname]);

%% 
figure('color','w','Position',[100,100, 450, 450])
ax = axes('Position', [0.2, 0.15, 0.7, 0.7])
%% More careful stats
dt = 1/3;

X = cat(2,nnvl{:})';
Y = dt\cat(2,f{:})';

Y = Y(X~=0);
X = X(X~=0);

%remove outliers
J = find((X>prctile(X, 0.1)) .* (X<prctile(X, 99.9)));

X = X(J);
Y = Y(J);

nbins=25;
[Xb, Yb, ~, ~, ~, ~, Xt, Yt] = BinData_v3(log(X), Y, nbins);

S = cell(1,nbins)
% calculate mean and sem with n=100 by repeated sampling
n=500;
for i=1:numel(Xt);
    sX=[];
    sY=[];
    for j=1:100;
        allX = Xt{i};
        allY = Yt{i};
        indSamp = datasample(1:numel(allX), n);
        sX=[sX, nanmean(allX(indSamp))];
        sY=[sY, nanmean(allY(indSamp))];
    end
    S{i}=[mean(sX), std(sX),mean(sY), std(sY)];
end

S = cat(1,S{:});

h = ploterr(exp(S(:,1)),S(:,3), [], S(:,4),'o','logx');

%h = ploterr(exp(Xb(1:end-1)), Yb(1:end-1), [], steYb(1:end-1),'o','logx');
BarTicklength(h(2),0)

for jj=1:numel(h)
    h(jj).Color = 'k';
    h(jj).MarkerFaceColor = 'k';
    h(jj).MarkerSize=7;
end

shg
hold on

%[BETA,resnorm,residual,exitflag,output,lambda,jacobian] = lsqcurvefit(@HillFunction,[1 100 0], exp(double(Xb([1 3:end-1]))), Yb([1 3:end-1]),[0 0 0],[inf inf 0 ])
[BETA,resnorm,residual,exitflag,output,lambda,jacobian] = lsqcurvefit(@HillFunction,[1 100 0], exp(double(S(:,1))), double(S(:,3)),[0 0 0],[3 inf 0.013 ])

%BETA= nlinfitWeight2(exp(double(Xb(1:end-1))), Yb(1:end-1), @HillFunction,[1 .15 0], steYb(1:end-1)+0.000001, [],[0 0 0],[inf 0.3 inf])
%BETA = [1.5 1000, 0]
xx = 0.5:1:1000;
yy = HillFunction(BETA, xx);
plot((xx), yy,'k')
set(gca,'XScale','log', 'XLim',[6 600],'YLim',[-0.01 0.6],'XTick',[1 10 100 1000])
xlabel('\nu - Nearest neighbors viral load (a.u.)')
ylabel('IR - Infection rate (h^{-1})')
ci = nlparci(BETA,residual,'jacobian',jacobian);

text(15, 0.4, '$IR(\nu) = 1.5\cdot\frac{1}{1+\frac{1000}{\nu}}$','interpreter','latex' );

y=Yb([1 3:end-1]);
yfit=HillFunction(BETA, exp(double(Xb([1 3:end-1]))));
SStot = sum((y-mean(y)).^2);                            % Total Sum-Of-Squares
SSres = sum((y(:)-yfit(:)).^2);                         % Residual Sum-Of-Squares
Rsq = 1-SSres/SStot;                                    % R^2
%% Save
figname = 'Fig3S_InfectionRateVsLoad'
set(gcf, 'PaperPositionMode','auto','color','w','InvertHardcopy','off')
print(gcf,'-depsc','-r600',[savpath '/epss/' figname]);
print(gcf,'-dpng','-r600',[savpath '/pngs/' figname]);



%% Load high MOI
BaseStr = regexprep([char(ispc.*'Z:\Images2019\') char(isunix.*'/bigstore/Images2019/')],char(0),'');
Usr = 'Jen';
Project = 'NFkBDynamics';
Dataset = 'TNFTitr_HighMOI_Dec122019_2019Dec12';
acquisition = 2;

% Get MD of raw data
acqname = ['acq_' num2str(acquisition)];
fpath = [BaseStr Usr filesep Project filesep Dataset filesep acqname];
MD=Metadata(fpath,[],1);
Wells = unique(MD.getSpecificMetadata('Position'));
frames = unique(cell2mat(MD.getSpecificMetadata('frame')));
dt = 1/3;%h

cTNF = ([100*(1/2).^[0:8], 0]);
% Load Results
R = MultiPositionSingleCellVirusResults(fpath)
Wells = R.PosNames;
frames = R.Frames;
t = [R.Frames-1]'.*dt;




%% Single track infection-death dt
close all

figure('color','w','Position',[100,100, 450, 450])

ax = axes('Position', [0.15, 0.73, 0.7, 0.18])

j=31
pos = R.PosNames{j};
Tracks1 = R.getTracks(pos)
cla
inds = find(logical(arrayfun(@(x) x.CellDies.*x.CellsGetInfected, Tracks1)));
i=38
%plot(Tracks1(inds(i)).T.*dt,[Tracks1(inds(i)).VirusTrack-virMean]./virStd,'g',Tracks1(inds(i)).T.*dt,[Tracks1(inds(i)).NucTrack-nucMean]./nucStd,'r',Tracks1(inds(i)).T.*dt,[Tracks1(inds(i)).DeathTrack-deathMean]./deathStd,'b','LineWidth',2); shg

yyaxis(ax,'left')
v = Smoothing((Tracks1(inds(i)).VirusTrack - min(Tracks1(inds(i)).VirusTrack))./(max(Tracks1(inds(i)).VirusTrack) - min(Tracks1(inds(i)).VirusTrack)),'neigh',5)
h(1) = plot(Tracks1(inds(i)).T.*dt,v,'Color',virusColor); shg
hold on
Iinf = Tracks1(inds(i)).indWhenCellGetsInfected;
plot([Iinf Iinf].*dt,[-2,3],'--','Color',virusColor)
set(gca,'YLim',[0,1], 'YTick', [0, 1],'ycolor',virusColor,'XLim', [0,48.3])
%ylabel('HSV-1 (a.u.)')

yyaxis(ax,'right')
d = Smoothing((Tracks1(inds(i)).DeathTrack - min(Tracks1(inds(i)).DeathTrack))./(max(Tracks1(inds(i)).DeathTrack) - min(Tracks1(inds(i)).DeathTrack)),'neigh',5)
h(2)=plot(Tracks1(inds(i)).T.*dt,d,'Color',deathColor); shg
hold on
Id = Tracks1(inds(i)).indWhenCellDies;
plot([Id Id].*dt,[-2,3],'--','Color',deathColor)
set(gca,'YLim',[0,1], 'YTick', [0, 1],'ycolor',deathColor)
%ylabel('Cytotox red(a.u.)')
%xlabel('t(h)')

%add shading
ptch = patch([Iinf Iinf Id Id]*dt, [0 1 1 0], [0.8 0.8 0.8])
ptch.EdgeColor = 'none';
ptch.FaceAlpha = 0.5;
t = title('100 ng/ml TNF')
ax = axes('Position', [0.15, 0.43, 0.7, 0.18])


j=36
pos = R.PosNames{j};
Tracks1 = R.getTracks(pos)
cla
inds = find(logical(arrayfun(@(x) x.CellDies.*x.CellsGetInfected, Tracks1)));
i=44%13
%plot(Tracks1(inds(i)).T.*dt,[Tracks1(inds(i)).VirusTrack-virMean]./virStd,'g',Tracks1(inds(i)).T.*dt,[Tracks1(inds(i)).NucTrack-nucMean]./nucStd,'r',Tracks1(inds(i)).T.*dt,[Tracks1(inds(i)).DeathTrack-deathMean]./deathStd,'b','LineWidth',2); shg

yyaxis(ax,'left')
v = Smoothing((Tracks1(inds(i)).VirusTrack - min(Tracks1(inds(i)).VirusTrack))./(max(Tracks1(inds(i)).VirusTrack) - min(Tracks1(inds(i)).VirusTrack)),'neigh',5)
h(1) = plot(Tracks1(inds(i)).T.*dt,v,'Color',virusColor); shg
hold on
Iinf = Tracks1(inds(i)).indWhenCellGetsInfected;
plot([Iinf Iinf].*dt,[-2,3],'--','Color',virusColor);
set(gca,'YLim',[0,1], 'YTick', [0, 1],'ycolor',virusColor,'XLim', [0,48.3])
ylabel('HSV-1 (a.u.)')

yyaxis(ax,'right')
d = Smoothing((Tracks1(inds(i)).DeathTrack - min(Tracks1(inds(i)).DeathTrack))./(max(Tracks1(inds(i)).DeathTrack) - min(Tracks1(inds(i)).DeathTrack)),'neigh',5)
h(2)=plot(Tracks1(inds(i)).T.*dt,d,'Color',deathColor); shg
hold on
Id = Tracks1(inds(i)).indWhenCellDies;
plot([Id Id].*dt,[-2,3],'--','Color',deathColor);
set(gca,'YLim',[0,1], 'YTick', [0, 1],'ycolor',deathColor)
ylabel('Death (a.u.)')
%xlabel('t(h)')

%add shading
ptch = patch([Iinf Iinf Id Id]*dt, [0 1 1 0], [0.8 0.8 0.8])
ptch.EdgeColor = 'none';
ptch.FaceAlpha = 0.5;
t = title('3 ng/ml TNF')

ax = axes('Position', [0.15, 0.13, 0.7, 0.18])


j=40
pos = R.PosNames{j};
Tracks1 = R.getTracks(pos)
cla
inds = find(logical(arrayfun(@(x) x.CellDies.*x.CellsGetInfected, Tracks1)));
i=16 %4
%plot(Tracks1(inds(i)).T.*dt,[Tracks1(inds(i)).VirusTrack-virMean]./virStd,'g',Tracks1(inds(i)).T.*dt,[Tracks1(inds(i)).NucTrack-nucMean]./nucStd,'r',Tracks1(inds(i)).T.*dt,[Tracks1(inds(i)).DeathTrack-deathMean]./deathStd,'b','LineWidth',2); shg

yyaxis(ax,'left')
v = Smoothing((Tracks1(inds(i)).VirusTrack - min(Tracks1(inds(i)).VirusTrack))./(max(Tracks1(inds(i)).VirusTrack) - min(Tracks1(inds(i)).VirusTrack)),'neigh',5)
h(1) = plot(Tracks1(inds(i)).T.*dt,v,'Color',virusColor); shg
hold on
Iinf = Tracks1(inds(i)).indWhenCellGetsInfected;
plot([Iinf Iinf].*dt,[-2,3],'--','Color',virusColor);
set(gca,'YLim',[0,1], 'YTick', [0, 1],'ycolor',virusColor,'XLim', [0,48.3])
%ylabel('HSV-1 (a.u.)')

yyaxis(ax,'right')
d = Smoothing((Tracks1(inds(i)).DeathTrack - min(Tracks1(inds(i)).DeathTrack))./(max(Tracks1(inds(i)).DeathTrack) - min(Tracks1(inds(i)).DeathTrack)),'neigh',5)
h(2)=plot(Tracks1(inds(i)).T.*dt,d,'Color',deathColor); shg
hold on
Id = Tracks1(inds(i)).indWhenCellDies;
plot([Id Id].*dt,[-2,3],'--','Color',deathColor);

set(gca,'YLim',[0,1], 'YTick', [0, 1],'ycolor',deathColor)
%ylabel('Cytotox red(a.u.)')

%add shading
ptch = patch([Iinf Iinf Id Id]*dt, [0 1 1 0], [0.8 0.8 0.8])
ptch.EdgeColor = 'none';
ptch.FaceAlpha = 0.5;


xlabel('Time (hours)')
t = title('0 ng/ml TNF')

%%
figname = 'Fig2_SingleCelldT'
set(gcf, 'PaperPositionMode','auto','color','w','InvertHardcopy','off','Renderer', 'painters')
print(gcf,'-depsc','-r600',[savpath '/epss/' figname]);
print(gcf,'-dpng','-r600',[savpath '/pngs/' figname]);


%% plot decays +/- virus for different [TNF]
close all
figure('color','w','Position',[100,100, 600, 300])
axes('Position', [0.1, 0.14, 0.37, 0.75])
t = R.Frames'/3;

a =(1:9)'*10.^(-2:0);
a = a(:);
emptyStr = {''; ''; '';''; ''; '';''; ''};
XTickLabel = ['.01'; emptyStr ;'.1'; emptyStr; '1'; emptyStr; '10'];
XTicks = log(a);

range = [31:40]
tzeva = (cmTNF(10));
A = t'-1/3;
for i=1:numel(range);
    pos = R.PosNames{range(i)};
    %deathMat = (R.putInMat(R.getTracks(R.PosNames{i}),'DeadSum'));
    %deathMat = R.putInMat(R.getTracks(pos),'Dead');
    %survivorFrac = sum(deathMat==0)./(sum(deathMat==1)+sum(deathMat==0));
    survivorFrac = R.getData('fracLive',pos);
    
    h=semilogy(t, (survivorFrac),'o');
    h.Color = tzeva(i,:);
    hold all
    
    BETA = R.getData('ExpFitParam',pos);
    y = Exponent(BETA,t);
    h1 = semilogy(t, y)
    h1.Color = tzeva(i,:);
    h1.LineWidth=1.5
    A=[A,survivorFrac',y'];
end
xlabel('time(h)')
ylabel('% live cells')
set(gca,'XLim',[0,max(t)],'YLim',[0.01,1],'YTick',[0.1 0.25 0.5 0.75 1],'YTickLabel',100*[0.1 0.25 0.5 0.75 1],'TickLength', [0.03 0.01], 'fontsize', 17)
title('+HSV-I')

axes('Position', [0.48, 0.14, 0.37, 0.75])
t = R.Frames'/3;
betais = []
range = [11:20]
tzeva = (cmTNF(10));
for i=1:numel(range);
    pos = R.PosNames{range(i)};
    %deathMat = (R.putInMat(R.getTracks(R.PosNames{i}),'DeadSum'));
    %deathMat = R.putInMat(R.getTracks(pos),'Dead');
    %survivorFrac = sum(deathMat==0)./(sum(deathMat==1)+sum(deathMat==0));
    
    survivorFrac = R.getData('fracLive',pos);
    h=semilogy(t, (survivorFrac),'o');
    h.Color = tzeva(i,:);
    hold all
    
    BETA = R.getData('ExpFitParam',pos);
    y = Exponent(BETA,t);
    h1 = semilogy(t, y)
    h1.Color = tzeva(i,:);
    h1.LineWidth=1.5
    A=[A,survivorFrac',y'];
end
xlabel('time(h)')
set(gca,'XLim',[0,max(t)],'YLim',[0.01,1],'YTick',[0.1 0.25 0.5 0.75 1],'YTickLabel',100*[],'TickLength', [0.03 0.01], 'fontsize', 17)
title('-HSV-I')

cb = colorbar('Position', [0.87 0.1400 0.05 0.7500]);
cb.Ticks = asinh(fliplr([100*(1/2).^[0:8], 0]))./asinh(100);
cb.TickLabels= num2str(fliplr([100*(1/2).^[0:8], 0])','%3.1f');
cb.Label.String = '[TNF\alpha] (ng/ml)'
cb.Label.FontSize=17
cb.Label.Position(1) = 0.2
cb.Label.Color='w'
colormap(flipud(cmTNF(10)))
%%
figname = 'Fig2_Sup_AllCellsDecay'
set(gcf, 'PaperPositionMode','auto','color','w','InvertHardcopy','off')
print(gcf,'-depsc','-r600',[savpath '/epss/' figname]);
print(gcf,'-dpng','-r600',[savpath '/pngs/' figname]);

%% Hists of death times withFit
close all
figure('color','w','Position',[100,100, 600, 300])

ax = axes('Position', [0.12, 0.2, 0.35, 0.7])

inds1 = 32;
edges = 1:1:50
for i=1:numel(inds1)
    pos = R.PosNames{inds1(i)};
    Tracks = R.getTracks(pos)
    
    inds = find(logical(arrayfun(@(x) x.CellDies.*x.CellsGetInfected, Tracks)));
    
    TDies = arrayfun(@(x) x.T(x.indWhenCellDies)*dt, Tracks,'UniformOutput',false);
    Tinf = arrayfun(@(x) x.T(x.indWhenCellGetsInfected)*dt, Tracks,'UniformOutput',false);
    
    
    h = histogram(cat(2,TDies{inds})-cat(2,Tinf{inds}),edges); hold on
    hold on
    yToFit = h.Values;
    xToFit = (h.BinEdges(1:end-1)+h.BinEdges(2:end))/2;
    BETA0 = [max(yToFit), 1/10, 0];
    BETA = lsqcurvefit(@(beta, xdata) Exponent(beta, xdata), BETA0, xToFit, yToFit)
    x = 0:0.1:150;
    y = Exponent(BETA,x);
    h = plot(x,y);
    h.LineWidth=3;
    R.getData('ExpFitParam',R.PosNames{inds1})
end
ax.XTick = dt*[0:30:150];
ax.XTickLabel = dt*[0:30:150];
ax.XLim=[0,50]
ax.YTick = 0:25:300;
ax.YLim = [0, 50]

xlabel('\Deltat_{infection\rightarrow death}(h)')
ylabel('# of Cells')


ax = axes('Position', [0.6, 0.2, 0.35, 0.7])
t = R.Frames'/3;

a =(1:9)'*10.^(-2:0);
a = a(:);
emptyStr = {''; ''; '';''; ''; '';''; ''};
XTickLabel = ['.01'; emptyStr ;'.1'; emptyStr; '1'; emptyStr; '10'];
XTicks = log(a);

range = inds1
tzeva = flipud(parula(numel(range)));
for i=1:numel(range);
    pos = R.PosNames{range(i)};
    
    survivorFrac = R.getData('fracLive',pos);
    h=semilogy(t, (survivorFrac),'o');
    hold all
        
    %interp from fit
    %BETA2Interp = polyval(p,BETA(2)+pRateInf);
    
    y = Exponent([4 BETA(2)+pRateInf],t);
    h1 = semilogy(t, y)
    h1.LineWidth=3;
    
    
end
xlabel('time(h)')
ylabel('% live cells')
ax.XTick = dt*[0:30:150];
ax.XTickLabel = dt*[0:30:150];
set(gca,'XLim',[0,50],'YLim',[0.01,1],'YTick',[0.1 0.5  1],'YTickLabel',100*[0.1 0.5 1],'TickLength', [0.03 0.01])

%%
figname = 'Fig2_HistAndDecayWithFit'
set(gcf, 'PaperPositionMode','auto','color','w','InvertHardcopy','off','Renderer', 'painters')
print(gcf,'-depsc','-r600',[savpath '/epss/' figname]);
print(gcf,'-dpng','-r600',[savpath '/pngs/' figname]);



%% Hists of death times withFit
close all
figure('color','w','Position',[100,100, 900, 600])

ax = axes('Position', [0.12, 0.2, 0.35, 0.7])

inds1 = [41 42:49];
edges = 0:1:48;
A=[];
for i=1:numel(inds1)
    ax = subplot(3,3,i);
    pos = R.PosNames{inds1(i)};
    Tracks = R.getTracks(pos);
    
    inds = find(logical(arrayfun(@(x) x.CellDies.*x.CellsGetInfected, Tracks)));
    
    TDies = arrayfun(@(x) x.T(x.indWhenCellDies)*dt, Tracks,'UniformOutput',false);
    Tinf = arrayfun(@(x) x.T(x.indWhenCellGetsInfected)*dt, Tracks,'UniformOutput',false);
    
    
    h = histogram(cat(2,TDies{inds})-cat(2,Tinf{inds}),edges);
    
    hold on
    yToFit = h.Values;
    xToFit = (h.BinEdges(1:end-1)+h.BinEdges(2:end))/2;
    
    
        
    BETA0 = [max(yToFit), 1/10, 0];
    
    if i==1
         BETA = lsqcurvefit(@(beta, xdata) Exponent(beta, xdata), BETA0, xToFit(3:end), yToFit(3:end))   ;    
    else
        BETA = lsqcurvefit(@(beta, xdata) Exponent(beta, xdata), BETA0, xToFit(2:end), yToFit(2:end));
    end
    BETA(2)
    x = 0:0.1:150;
    y = Exponent(BETA,x);
    
    A=[A,xToFit',yToFit',Exponent(BETA,xToFit)'];
    
    h = plot(x,y);
    h.LineWidth=3;
    %R.getData('ExpFitParam',R.PosNames{inds1(i)});
    ax.XTick = dt*[0:36:150];
    ax.XTickLabel = dt*[0:36:150];
    ax.XLim=[0,50];
    ax.YScale='log';
    ax.YTick = [0.1 1 10 100];
    ax.YLim = [0.5, 100];
    title(['TNF=' num2str(cTNF(i),'%.2f') 'ng/ml']);
end


xl = xlabel('\Deltat_{infection\rightarrow death}(h)');
yl = ylabel('# of Cells');
xl.Position(1)=-110;
%xl.Position(2)=0.015;
yl.Position(1)=-150;
yl.Position(2)=50;
a1 = annotation('arrow',[0.1, 0.5], [0.06, 0.06]);
a2 = annotation('arrow',[0.06, 0.06], [0.1, 0.5]);


%%
figname = 'Fig_ExpTails'
set(gcf, 'PaperPositionMode','auto','color','w','InvertHardcopy','off','Renderer', 'painters')
print(gcf,'-depsc','-r600',[savpath '/epss/' figname]);
print(gcf,'-dpng','-r600',[savpath '/pngs/' figname]);


%% Hists of death times withFit
close all
figure('color','w','Position',[100,100, 800, 300])

%ax = axes('Position', [0.1, 0.2, 0.28, 0.7])


BetasFromGlobal = []
BetasFromSingle = []

inds1 = 31:60;
edges = 1:1:50
tzeva = flipud(cmTNF(10))
tzeva  = [tzeva; tzeva; tzeva];
for i=1:numel(inds1)
    pos = R.PosNames{inds1(i)};
    Tracks = R.getTracks(pos)
    
    inds = find(logical(arrayfun(@(x) x.CellDies.*x.CellsGetInfected, Tracks)));

    TDies = arrayfun(@(x) x.T(x.indWhenCellDies)*dt, Tracks,'UniformOutput',false);
    Tinf = arrayfun(@(x) x.T(x.indWhenCellGetsInfected)*dt, Tracks,'UniformOutput',false);

    [Values, BinEdges] = histcounts(cat(2,TDies{inds})-cat(2,Tinf{inds}),edges);

    %h.FaceColor = tzeva(i,:);
    %hold on
    yToFit = Values;
    xToFit = (BinEdges(1:end-1)+BinEdges(2:end))/2;
    BETA0 = [max(yToFit), 1/10, 0];
    
    %starting fit from 3rd bin
    BETA = lsqcurvefit(@(beta, xdata) Exponent(beta, xdata), BETA0, xToFit, yToFit);
    x = 0:0.1:150;
    %y = Exponent(BETA,x);
    %h = plot(x,y)
    %h.Color = tzeva(i,:);
    B = R.getData('ExpFitParam',R.PosNames{inds1(i)})
    BetasFromGlobal = [BetasFromGlobal B(2)+pRateInf]; %correction for prolif
    BetasFromSingle = [BetasFromSingle BETA(2)];
    %h.LineWidth=3;
    %ax.YLim=[0, 100];
end
% ax.XTick = dt*[0:30:150];
% ax.XTickLabel = dt*[0:30:150];
% ax.XLim=[0,50]
% ax.YTick = 0:25:300;
% ax.YLim = [0, 100]

%xlabel('death time(h)');
%xlabel('\Deltat_{infection\rightarrow death}(h)')
%ylabel('# of Cells');

ax = axes('Position', [0.1, 0.2, 0.28, 0.7])



scatter(BetasFromSingle,BetasFromGlobal, [], tzeva, 'filled');
hold on
x = [0, 0.19];

[p, S] = polyfit(BetasFromSingle,BetasFromGlobal,1);

[y_fit,delta] = polyval(p,x,S);

plot(x,y_fit,'Color',[220,20,60]/255);

%plot(x,y_fit+2*delta,'r--',x,y_fit-2*delta,'r--')

box on
ax.YLim = [0, 0.18];
ax.XLim = [0, 0.16];
ax.XTick = 0:0.04:0.2;
ax.YTick = 0:0.04:0.2;

Err = polyparci(p,S, .68);
clear range; ErrR = range(Err)/2;
text(0.02, 0.15, ['Slope=' num2str(p(1),3) '\pm' num2str(ErrR(1),1) ], 'FontSize', 14,'Color',[220,20,60]/255);
xlabel('death rates from 1/\Deltat')
ylabel('death rates from global')

cb = colorbar('Position', [0.41 0.200 0.0324 0.700]);
cb.Ticks = asinh(fliplr([100*(1/2).^[0:8], 0]))./asinh(100);
cb.TickLabels= num2str(fliplr([100*(1/2).^[0:8], 0])','%3.1f');
cb.Label.String = '[TNF\alpha] (ng/ml)'
cb.Label.FontSize=17
cb.Label.Position(1) = 0.2
cb.Label.Color='w'
colormap(cmTNF(10))
%%
figname = 'Fig2_CompareAllGlobalvdT'
set(gcf, 'PaperPositionMode','auto','color','w','InvertHardcopy','off')
print(gcf,'-depsc','-r600',[savpath '/epss/' figname]);
print(gcf,'-dpng','-r600',[savpath '/pngs/' figname]);


%%  Rate vs TNF infected
allParams = R.getData('ExpFitParam');
allParams = cat(1,allParams{:})
allParams = allParams(:,2);
allParams = reshape(allParams,[],6);

LinearRange = [0.25]

figure('color','w','Position',[100,100, 900, 300])

ax = axes('Position', [0.12, 0.2, 0.25, 0.7])

x = cTNF'
y = mean(allParams(:,4:6),2)+pRateInf;

sigmay = std(allParams(:,4:6),[],2);

h = ploterr(asinh(cTNF/LinearRange),y,[], sigmay,'ko','hhy',0,'hhx',0);

for jj=1:numel(h)
    h(jj).Color = 'k';
    h(jj).MarkerFaceColor = 'k';
    h(jj).MarkerSize=7;
end

ylabel({'Infected death' 'rate - \beta_{{i}}'})
xlabel('')


a1 = sort(- (1:9)'*10.^(1:2));
a2 =(1:9)'*10.^(1:5);
a = [sort(a1(:)); 0; a2(:)]/1000;
emptyStr = {''; ''; '';''; ''; '';''; ''};
XTickLabel = [emptyStr; ''; emptyStr; {''}; '0'; {''}; emptyStr; ' 0.1'; emptyStr; ' 1'; emptyStr; ' 10'; emptyStr; '    100'; emptyStr];
XTicks = asinh(a/LinearRange);
set(gca, 'xtick', XTicks, 'xticklabel', XTickLabel)

ax.XLim = asinh([-0.1, 150]/LinearRange)

ax.YTick = 0:0.04:0.2;



% ax = axes('Position', [0.15, 0.14, 0.2, 0.39])
% 
% 
% range = [1:30]
% tToMeasure = 36;%h
% survAtT = [];
% for i=1:numel(range);
%     pos = R.PosNames{range(i)};
%     
%     survivorFrac = R.getData('fracLive',pos);
%     
%     survAtT = [survAtT survivorFrac(tToMeasure/dt)];
% end
% 
% falsePosAtT =  1-reshape(survAtT,[],3);
% 
% tzeva = flipud(parula(numel(cTNF)));
% 
% 
% h = ploterr(asinh(cTNF/LinearRange), mean(falsePosAtT(:,:),2), [], std(falsePosAtT(:,:),[],2),'o','hhy',0,'hhx',0);
% hold on
% 
% for jj=1:numel(h)
%     h(jj).Color = 'k';
%     h(jj).MarkerFaceColor = 'k';
% end
% 
% 
% set(gca, 'xtick', XTicks, 'xticklabel', XTickLabel)
% 
 xlabel('[TNF] ng/ml')
% ylabel({'False positive' 'fraction'})
% shg
% ax.XLim = asinh([-0.1, 150]/LinearRange)
% ax.YLim = [-0.03, 0.8]
% 


%%
figname = 'Fig2_deathRatevsTNF'
set(gcf, 'PaperPositionMode','auto','color','w','InvertHardcopy','off')
print(gcf,'-depsc','-r600',[savpath '/epss/' figname]);
print(gcf,'-dpng','-r600',[savpath '/pngs/' figname]);


%%  Rate vs TNF infected
allParams = R.getData('ExpFitParam');
allParams = cat(1,allParams{:})
allParams = allParams(:,2);
allParams = reshape(allParams,[],6);

LinearRange = [0.25]

figure('color','w','Position',[100,100, 900, 300])

% ax = axes('Position', [0.12, 0.2, 0.25, 0.7])
% 
% x = cTNF'
% y = mean(allParams(:,4:6),2);
% 
% sigmay = std(allParams(:,4:6),[],2);
% 
% h = ploterr(asinh(cTNF/LinearRange),y,[], sigmay,'ko','hhy',0,'hhx',0);
% 
% for jj=1:numel(h)
%     h(jj).Color = 'k';
%     h(jj).MarkerFaceColor = 'k';
% end
% 
% ylabel({'Infected death' 'rate - \beta_{{i}}'})
% xlabel('')
% 
% 
% a1 = sort(- (1:9)'*10.^(1:2));
% a2 =(1:9)'*10.^(1:5);
% a = [sort(a1(:)); 0; a2(:)]/1000;
% emptyStr = {''; ''; '';''; ''; '';''; ''};
% XTickLabel = [emptyStr; ''; emptyStr; {''}; '0'; {''}; emptyStr; ' 0.1'; emptyStr; ' 1'; emptyStr; ' 10'; emptyStr; '    100'; emptyStr];
% XTicks = asinh(a/LinearRange);
% set(gca, 'xtick', XTicks, 'xticklabel', XTickLabel)
% 
% ax.XLim = asinh([-0.1, 150]/LinearRange)
% 
% ax.YTick = 0:0.04:0.2;
% 


ax = axes('Position',  [0.12, 0.2, 0.25, 0.7])


range = [1:30]
tToMeasure = 24;%h
survAtT = [];
for i=1:numel(range);
    pos = R.PosNames{range(i)};
    
    survivorFrac = R.getData('fracLive',pos);
    
    survAtT = [survAtT survivorFrac(tToMeasure/dt)];
end

falsePosAtT =  1-reshape(survAtT,[],3);

tzeva = flipud(parula(numel(cTNF)));


h = ploterr(asinh(cTNF/LinearRange), mean(falsePosAtT(:,:),2), [], std(falsePosAtT(:,:),[],2),'o','hhy',0,'hhx',0);
hold on

for jj=1:numel(h)
    h(jj).Color = 'k';
    h(jj).MarkerFaceColor = 'k';
    h(jj).MarkerSize=7;
end


set(gca, 'xtick', XTicks, 'xticklabel', XTickLabel,  'ytick', 0:0.2:1, 'yticklabel', 100*[0:0.2:1])

 xlabel('[TNF] ng/ml')
ylabel({'False positive' 'percent'})
shg
ax.XLim = asinh([-0.1, 150]/LinearRange)
ax.YLim = [-0.03, 0.6]

a1 = sort(- (1:9)'*10.^(1:2));
a2 =(1:9)'*10.^(1:5);
a = [sort(a1(:)); 0; a2(:)]/1000;
emptyStr = {''; ''; '';''; ''; '';''; ''};
XTickLabel = [emptyStr; ''; emptyStr; {''}; '0'; {''}; emptyStr; ' 0.1'; emptyStr; ' 1'; emptyStr; ' 10'; emptyStr; '  100'; emptyStr];
XTicks = asinh(a/LinearRange);
set(gca, 'xtick', XTicks, 'xticklabel', XTickLabel)



%%
figname = 'Fig2_falsePosvsTNF'
set(gcf, 'PaperPositionMode','auto','color','w','InvertHardcopy','off')
print(gcf,'-depsc','-r600',[savpath '/epss/' figname]);
print(gcf,'-dpng','-r600',[savpath '/pngs/' figname]);


%%  S/A tradeoff rate vs frac at 36h
allParams = R.getData('ExpFitParam');
allParams = cat(1,allParams{:})
allParams = allParams(:,2);
allParams = reshape(allParams,[],6);

LinearRange = [0.25]

figure('color','w','Position',[100,100, 900, 300])

% ax = axes('Position', [0.15, 0.56, 0.2, 0.39])
% 
% x = cTNF'
% y = mean(allParams(:,4:6),2);
% 
% sigmay = std(allParams(:,4:6),[],2);
% 
% h = ploterr(asinh(cTNF/LinearRange),y,[], sigmay,'ko','hhy',0,'hhx',0);
% 
% for jj=1:numel(h)
%     h(jj).Color = 'k';
%     h(jj).MarkerFaceColor = 'k';
% end
% 
% ylabel({'Infected death' 'rate - \beta_{{i}}'})
% xlabel('')
% 
% 
% a1 = sort(- (1:9)'*10.^(1:2));
% a2 =(1:9)'*10.^(1:5);
% a = [sort(a1(:)); 0; a2(:)]/1000;
% emptyStr = {''; ''; '';''; ''; '';''; ''};
% XTickLabel = [emptyStr; ''; emptyStr; {''}; '0'; {''}; emptyStr; ' 0.1'; emptyStr; ' 1'; emptyStr; ' 10'; emptyStr; '    100'; emptyStr];
% XTicks = asinh(a/LinearRange);
% set(gca, 'xtick', XTicks, 'xticklabel', [])
% 
% ax.XLim = asinh([-0.1, 150]/LinearRange)
% 
% ax.YTick = 0:0.04:0.2;
% 
% 
% 
% ax = axes('Position', [0.15, 0.14, 0.2, 0.39])
% 
% 
% range = [1:30]
% tToMeasure = 36;%h
% survAtT = [];
% for i=1:numel(range);
%     pos = R.PosNames{range(i)};
%     
%     survivorFrac = R.getData('fracLive',pos);
%     
%     survAtT = [survAtT survivorFrac(tToMeasure/dt)];
% end
% 
% falsePosAtT =  1-reshape(survAtT,[],3);
% 
% tzeva = flipud(parula(numel(cTNF)));
% 
% 
% h = ploterr(asinh(cTNF/LinearRange), mean(falsePosAtT(:,:),2), [], std(falsePosAtT(:,:),[],2),'o','hhy',0,'hhx',0);
% hold on
% 
% for jj=1:numel(h)
%     h(jj).Color = 'k';
%     h(jj).MarkerFaceColor = 'k';
% end
% 
% 
% set(gca, 'xtick', XTicks, 'xticklabel', XTickLabel)
% 
% xlabel('[TNF] ng/ml')
% ylabel({'False positive' 'fraction'})
% shg
% ax.XLim = asinh([-0.1, 150]/LinearRange)
% ax.YLim = [-0.03, 0.8]
% 

ax = axes('Position', [0.45, 0.2, 0.25, 0.7])

tzeva = flipud(cmTNF(10));%flipud(cmTNF(numel(cTNF)));

for i=1:numel(cTNF)
    
    h = ploterr(mean(allParams(i,4:6),2), mean(falsePosAtT(i,:),2), std(allParams(i,4:6),[],2), std(falsePosAtT(i,:),[],2),'o','hhy',0,'hhx',0);
    hold on
    
    for jj=1:numel(h)
        h(jj).Color = tzeva(i,:);
        h(jj).MarkerFaceColor = tzeva(i,:);
        h(jj).MarkerSize=10;
    end
    
end
xlabel('Infected death rate - \beta_{{i}}')
ylabel('False positive fraction')
shg
ax.XLim = [0, 0.17]
ax.YLim = [-0.05, 0.6]

ax.XTick = 0:0.04:0.2;

cb = colorbar('Position', [0.72 0.1500 0.03 0.800]);
cb.Ticks = asinh(fliplr([100*(1/2).^[0:8], 0]))./asinh(100);
cb.TickLabels= num2str(fliplr([100*(1/2).^[0:8], 0])','%3.1f');
cb.Label.String = '[TNF\alpha] (ng/ml)'
cb.Label.FontSize=17
cb.Label.Position(1) = 0.2
cb.Label.Color='w'
colormap(cmTNF(10))

%%
figname = 'Fig2_SAToff'
set(gcf, 'PaperPositionMode','auto','color','w','InvertHardcopy','off')
print(gcf,'-depsc','-r600',[savpath '/epss/' figname]);
print(gcf,'-dpng','-r600',[savpath '/pngs/' figname]);


%% Plot S/A tradeoff rates
close all
figure('color','w','Position',[100,100, 300, 300])
axes('Position', [0.27, 0.2, 0.65, 0.65])
allParams = R.getData('ExpFitParam');
allParams = cat(1,allParams{:})

allParams = allParams(:,2);
allParams = reshape(allParams,[],6);

tzeva = flipud(cmTNF(10));%flipud(cmTNF(numel(cTNF)));

for i=1:numel(cTNF)
    h = ploterr(mean(allParams(i,4:6)+pRateInf,2), mean(allParams(i,1:3)+pRateUninf,2), std(allParams(i,4:6),[],2), std(allParams(i,1:3),[],2),'o');
    hold on
    for jj=1:numel(h)
        h(jj).Color = tzeva(i,:);
        h(jj).MarkerFaceColor = tzeva(i,:);
        h(jj).LineWidth = 2;
        h(jj).MarkerSize = 7;
    end
    
end
set(gca,'YLim',[0.02, 0.08],'XLim',[0, 0.18], 'YTick', [0:0.02:0.1], 'XTick', [0:0.05:0.2])
xlabel('Infected death rate - \beta_{{i}}')
ylabel('Bystander death rate - \beta_{{b}}')
%%
figname = 'Fig3S_SAToffRates'
set(gcf, 'PaperPositionMode','auto','color','w','InvertHardcopy','off')
print(gcf,'-depsc','-r600',[savpath '/epss/' figname]);
print(gcf,'-dpng','-r600',[savpath '/pngs/' figname]);

%% Plot both uninfected and infected cell death rates

allParams = R.getData('ExpFitParam');
allParams = cat(1,allParams{:})
allParams = allParams(:,2);
allParams = reshape(allParams,[],6);
figure('color','w','Position',[100,100, 450, 450])

ax = axes('Position', [0.2, 0.15, 0.7, 0.7])
LinearRange = [0.1]
LinearRangeY = [0.3]

hold on
x = cTNF'
y = asinh((mean(allParams(:,4:6),2))/LinearRangeY); %Infected
y2 = asinh((mean(allParams(:,1:3),2))/LinearRangeY); %Uninfected
%Error calculation: asinh' = 1/sqrt(1+x^2)
sigmay = std(allParams(:,4:6),[],2)./sqrt(mean(allParams(:,4:6)./LinearRangeY,2).^2+1)./LinearRangeY;
sigmay2 = std(allParams(:,1:3),[],2)./sqrt(mean(allParams(:,1:3)./LinearRangeY,2).^2+1)./LinearRangeY;
h = errorbar(asinh(cTNF/LinearRange)',y', sigmay','-ko', 'linewidth', 2, 'markersize', 7, 'markerfacecolor', 'k');
h2 = errorbar(asinh(cTNF/LinearRange)',y2',sigmay2','-ko', 'linewidth', 2, 'markersize', 7, 'markerfacecolor', 'b');
a1 = sort(- (1:9)'*10.^(1:2));
a2 =(1:9)'*10.^(1:5);
a = [sort(a1(:)); 0; a2(:)]/1000;
emptyStr = {''; ''; '';''; ''; '';''; ''};
XTickLabel = [emptyStr; ''; emptyStr; {''}; '0'; {''}; emptyStr; ' 0.1'; emptyStr; ' 1'; emptyStr; ' 10'; emptyStr; '    100'; emptyStr];
XTicks = asinh(a/LinearRange);

a1 = sort(- (1:9)'*10.^(1:2));
a2 =(1:9)'*10.^(1:5);
a = [sort(a1(:)); 0; a2(:)]/10000;
emptyStr = {''; ''; '';''; ''; '';''; ''};
YTickLabel = [emptyStr; ' '; emptyStr; {''}; '0'; {''}; emptyStr; ' 0.01'; emptyStr; ' .1'; emptyStr; ' 1'; emptyStr; '    10'; emptyStr];

YTicks = asinh(a/LinearRangeY);

set(gca, 'xtick', XTicks, 'xticklabel', XTickLabel,  'ytick', YTicks, 'yticklabel', YTickLabel, 'xlim', [-0.5 8])
%hl = legend('Infected cells', 'Uninfected cells')
%hl.Box = 'off'
%hl.Position(1) = hl.Position(1)-0.3;
ylabel('Infected death rate - \beta_i')
xlabel('[TNF] ng/ml')
box on
%title('Uninfected and Infected death rates')

%%
figname = 'Fig2_RatesWVirus'
set(gcf, 'PaperPositionMode','auto','color','w','InvertHardcopy','off')
print(gcf,'-depsc','-r600',[savpath '/epss/' figname]);
print(gcf,'-dpng','-r600',[savpath '/pngs/' figname]);



%% ROC at TP with errorbars...
figure('color','w','Position',[100,100, 550, 450])
ax = axes('Position', [0.15, 0.15, 0.7*45/55, 0.7])

fracLive = R.getData('fracLive')
numLive = R.getData('numLive')

rangunInf = 1:30;
rangInf = 30+rangunInf;


%tCells = cat(1, numLive{:})'./cat(1, fracLive{:})';

%totalCells = sum(reshape(tCells,145,10,[]),3);
%totalCellsInf = sum(reshape(tCells(:,rangInf),145,10,[]),3);
%totalCellsUninf = sum(reshape(tCells(:,rangunInf),145,10,[]),3);

%totalTPCells = tCells(:,rangInf).*(1-cat(1, fracLive{rangInf})');
%totalTPCells = sum(reshape(totalTPCells,145,10,[]),3);
%totalFPCells = tCells(:,rangunInf).*(1-cat(1, fracLive{rangunInf})');
%totalFPCells = sum(reshape(totalFPCells,145,10,[]),3);

fracTP = reshape(1-cat(1, fracLive{rangInf})',145,10,[]);
fracFP = reshape(1-cat(1, fracLive{rangunInf})',145,10,[]);

%TP = totalTPCells./totalCellsInf;

%FP = totalFPCells./totalCellsUninf;
meanFP = mean(fracFP,3);
meanTP = mean(fracTP,3);
stdFP = std(fracFP,[],3);
stdTP = std(fracTP,[],3);
tp = 20*3;

h = ploterr(meanFP(tp,:),meanTP(tp,:),stdFP(tp,:),stdTP(tp,:));

for i=1:numel(h)
    h(i).Color = 'k'
end
    BarTicklength(h(2),0)
    BarTicklengthY(h(3),0)
tzeva = cmTNF(size(meanTP,2))

%Prec = TP./(TP+FP);
hold on
for i=1:(size(meanTP,2))
    h = plot(meanFP(tp,i),meanTP(tp,i),'o')
   hold on ;
    h.MarkerFaceColor = tzeva(i,:);
    h.MarkerSize = 10
    h.Color = tzeva(i,:)
end
shg
plot([0,1], [0,1],'--k')
set(gca,'xlim',[0,1],'ylim',[0,1],'XTick',0:0.2:1,'YTick',0:0.2:1)
xlabel('False positive rate')
ylabel('True positive rate')

cb = colorbar('Position', [0.8 0.1500 0.06 0.700]);
cb.Ticks = asinh(fliplr([100*(1/2).^[0:8], 0]))./asinh(100);
cb.TickLabels= num2str(fliplr([100*(1/2).^[0:8], 0])','%3.1f');
cb.Label.String = '[TNF\alpha] (ng/ml)'
cb.Label.FontSize=17
cb.Label.Position(1) = 0.2
cb.Label.Color='w'
colormap(flipud(cmTNF(20)))
%%
figname = 'Fig2_ROC_Curve'
set(gcf, 'PaperPositionMode','auto','color','w','InvertHardcopy','off')
print(gcf,'-depsc','-r600',[savpath '/epss/' figname]);
print(gcf,'-dpng','-r600',[savpath '/pngs/' figname]);



%% Plot both uninfected and infected cell death rates

allParams = R.getData('ExpFitParam');
allParams = cat(1,allParams{:})
allParams = allParams(:,2);
allParams = reshape(allParams,[],6);
figure('color','w','Position',[100,100, 450, 450])

ax = axes('Position', [0.2, 0.15, 0.7, 0.7])
LinearRange = [0.1]
LinearRangeY = [0.05]

hold on
x = cTNF'
y = asinh((mean(allParams(:,4:6),2)+pRateInf)/LinearRangeY); %Infected
y2 = asinh((mean(allParams(:,1:3),2)+pRateUninf)/LinearRangeY); %Uninfected
%Error calculation: asinh' = 1/sqrt(1+x^2)
sigmay = std(allParams(:,4:6),[],2)./sqrt(mean(allParams(:,4:6)./LinearRangeY,2).^2+1)./LinearRangeY;
sigmay2 = std(allParams(:,1:3),[],2)./sqrt(mean(allParams(:,1:3)./LinearRangeY,2).^2+1)./LinearRangeY;
h = errorbar(asinh(cTNF/LinearRange)',y', sigmay','-o', 'linewidth', 2, 'markersize', 7,'color',virusColor, 'markerfacecolor', virusColor,'markeredgecolor','none');
h2 = errorbar(asinh(cTNF/LinearRange)',y2',sigmay2','-o', 'linewidth', 2, 'markersize', 7,'color',liveColor, 'markerfacecolor', liveColor,'markeredgecolor','none');
a1 = sort(- (1:9)'*10.^(1:2));
a2 =(1:9)'*10.^(1:5);
a = [sort(a1(:)); 0; a2(:)]/1000;
emptyStr = {''; ''; '';''; ''; '';''; ''};
XTickLabel = [emptyStr; ''; emptyStr; {''}; '0'; {''}; emptyStr; ' 0.1'; emptyStr; ' 1'; emptyStr; ' 10'; emptyStr; '    100'; emptyStr];
XTicks = asinh(a/LinearRange);

a1 = sort(- (1:9)'*10.^(1:2));
a2 =(1:9)'*10.^(1:5);
a = [sort(a1(:)); 0; a2(:)]/10000;
emptyStr = {''; ''; '';''; ''; '';''; ''};
YTickLabel = [emptyStr; ' '; emptyStr; {''}; '0'; {''}; emptyStr; ' 0.01'; emptyStr; ' .1'; emptyStr; ' 1'; emptyStr; '    10'; emptyStr];

YTicks = asinh(a/LinearRangeY);

set(gca, 'fontsize', 17, 'xtick', XTicks, 'xticklabel', XTickLabel,  'ytick', YTicks, 'yticklabel', YTickLabel, 'xlim', [-0.5 8])
hl = legend('Infected cells', 'Uninfected cells')
hl.Box = 'off'
hl.Position(1) = hl.Position(1)-0.3;
ylabel('Death rate - \beta')
xlabel('[TNF] ng/ml')
box on
%title('Uninfected and Infected death rates')

%%
figname = 'Fig3S_DeathRatesW_WOVirus'
set(gcf, 'PaperPositionMode','auto','color','w','InvertHardcopy','off')
print(gcf,'-depsc','-r600',[savpath '/epss/' figname]);
print(gcf,'-dpng','-r600',[savpath '/pngs/' figname]);


%% Viral growth
close all
figure('color','w','Position',[100,100, 300, 300])
t = R.Frames'/3;

ax = axes('Position', [0.24, 0.18, 0.7, 0.7])

indsInfNoTNF = [40, 50, 60]
%plot(t,cellfun(@(x) mean(x.Ints90('Cyan')), R.WellLbls{20}))
%hold on
viralLoadToFit = 1000*[cellfun(@(x) mean(x.Ints90('Red')), R.WellLbls{50})];
viralLoaderrToFit = 1000*[cellfun(@(x) std(x.Ints90('Red')/sqrt(x.num)), R.WellLbls{50})];

h = ploterr(t,viralLoadToFit,[],viralLoaderrToFit,'o')

for jj=1:numel(h)
    h(jj).Color = virusColor;
    h(jj).MarkerFaceColor = virusColor;
    h(jj).MarkerSize = 7;
end

hold on;
%
set(gca,'Xlim', [0 20], 'Ylim', [0 200])

J = 20:50;
[p, S] = polyfit(t(J),viralLoadToFit(J)',1);

[y_fit,delta] = polyval(p,t(J),S);
plot(t(J),y_fit,'Color','k');

Err = polyparci(p,S, .68);
clear range; ErrR = range(Err)/2;
text(1, 10, ['Slope=' num2str(p(1),2) 'u/h'],'Color','k');



y=viralLoadToFit(J)';
yfit=y_fit;
SStot = sum((y-mean(y)).^2);                            % Total Sum-Of-Squares
SSres = sum((y(:)-yfit(:)).^2);                         % Residual Sum-Of-Squares
Rsq = 1-SSres/SStot;                                    % R^2


ylabel('Average viral load (a.u.)')
xlabel('Time (hours)')

%%
figname = 'Fig3S_ViralGrowth'
set(gcf, 'PaperPositionMode','auto','color','w','InvertHardcopy','off')
print(gcf,'-depsc','-r600',[savpath '/epss/' figname]);
print(gcf,'-dpng','-r600',[savpath '/pngs/' figname]);


%%

close all
figure('color','w','Position',[100,100, 450, 450])

ax = axes('Position', [0.24, 0.18, 0.7, 0.7])

tzeva = cmTNF(10);
hs=[];
for j=1:10
    jj=[j, j+10, j+20]+30
pos = R.PosNames(jj);
Tracks1 = R.getTracks(pos)


A = arrayfun(@(x) x.VirusTrack(x.indWhenCellDies), Tracks1,'UniformOutput', false)%
%h = histogram(cat(1,A{cellfun(@(x) ~isempty(x), A)}),edges);
%h.FaceColor=tzeva(j,:);
A = 1000.*cat(1,A{cellfun(@(x) ~isempty(x), A)});
edges=[min(A):25:500];

h = histcounts(A,edges);
hs=[hs, h'];
centers=(edges(2:end)+edges(1:end-1))/2;
%h1 = plot(centers,asinh((h./sum(h))/0.01));
h1 = semilogy(centers,h./sum(h),'-o');

h1.Color=tzeva(j,:);

xlabel('Viral load at death(a.u.)')
ylabel('probability')
hold on
shg
end

%%
close all
figure('color','w','Position',[100 158 1111 400])

%ax = axes('Position', [0.24, 0.18, 0.7, 0.7])

tzeva = cmTNF(10);
hs=[];
sp=cell(10,1);
for j=1:10
    sp{j}= axes('Position',[0.1+rem(j-1,5)*0.16,0.6-0.4*floor((j-1)/5) ,0.14, 0.35]);
    jj=[j, j+10, j+20]+30
pos = R.PosNames(jj);
Tracks1 = R.getTracks(pos)

Tracks1 = Tracks1(arrayfun(@(x) any(x.Infected),Tracks1))
A = arrayfun(@(x) x.VirusTrack(x.indWhenCellDies), Tracks1,'UniformOutput', false);
%h = histogram(cat(1,A{cellfun(@(x) ~isempty(x), A)}),50,'Normalization', 'probability');
B = arrayfun(@(x) x.T(x.indWhenCellDies), Tracks1,'UniformOutput', false);

A = 1000*cat(1,A{cellfun(@(x) ~isempty(x), A)});
B = cat(1,B{cellfun(@(x) ~isempty(x), B)});
%A=A(B>20);
%B=B(B>20);

scatter(B/3,asinh(A),5,tzeva(j,:),'filled')
[xx,yy, stdXb, stdYb, steXb, steYb,xt] = BinData_v1(B,asinh(A),20)
% J=cellfun(@numel,xt)>20;
%  h1=ploterr(yy(J)/3,xx(J),steYb(J)/3,steXb(J));
%  h1(1).Color=tzeva(j,:);
%  h1(2).Color=tzeva(j,:);
%  h1(3).Color=tzeva(j,:);
hold on
 h1=plot(xx/3,yy,'-');
 h1(1).Color='k';
 h1(1).MarkerSize=5;
 h1(1).MarkerFaceColor=tzeva(j,:);
set(gca,'ylim',asinh([40,300]),'xlim',[0,48])
set(gca,'YTickLabel','','YTick',asinh([50, 100:100:500]),'XTick', 0:12:48 ,'XTickLabel', '' )
if j==1
    set(gca,'YTickLabel',[50, 100:100:500],'YTick',asinh([50, 100:100:500]),'XTick','')

end
if floor((j-1)/5)
    set(gca,'XTick', 0:12:48 ,'XTickLabel', 0:12:48)
end
if j==6
    yl = ylabel('Viral load at death(a.u.)')
    xl = xlabel('time of death (h)')
    set(gca,'YTickLabel',[50, 100:100:500],'YTick',asinh([50, 100:100:500]),'XTick', 0:12:48 ,'XTickLabel', 0:12:48)
end

shg
end
axes(sp{6})
%ylabel('Viral load at death(a.u.)')
%xlabel('time of death (h)')
cb = colorbar('Position', [0.91 0.1500 0.03 0.800]);
cb.Ticks = asinh(fliplr([100*(1/2).^[0:8], 0]))./asinh(100);
cb.TickLabels= num2str(fliplr([100*(1/2).^[0:8], 0])','%3.1f');
cb.Label.String = '[TNF\alpha] (ng/ml)'
cb.Label.FontSize=17
cb.Label.Position(1) = 0.2
cb.Label.Color='w'
colormap(flipud(cmTNF(10)))

a1 = annotation('arrow',[0.1, 0.5], [0.1, 0.1]);
a2 = annotation('arrow',[0.05, 0.05], [0.2, 0.6]);
yl.Position(1) = -20;
yl.Position(2)=6.5
xl.Position(1)=50
xl.Position(2)=3.6
%%
figname = 'ViralLoadAtDeathSingleCells'
set(gcf, 'PaperPositionMode','auto','color','w','InvertHardcopy','off')
print(gcf,'-depsc','-r600',[savpath '/epss/' figname]);
print(gcf,'-dpng','-r600',[savpath '/pngs/' figname]);
%%
figure('color','w','Position',[100,100, 450, 450])

ax = axes('Position', [0.24, 0.18, 0.7, 0.7])

tzeva = cmTNF(10);
hs=[];
for j=1:10
    jj=[j, j+10, j+20]+30
pos = R.PosNames(jj);
Tracks1 = R.getTracks(pos)
Tracks1 = Tracks1(arrayfun(@(x) any(x.Infected),Tracks1))

A = arrayfun(@(x) x.VirusTrack(x.indWhenCellDies), Tracks1,'UniformOutput', false);
%h = histogram(cat(1,A{cellfun(@(x) ~isempty(x), A)}),50,'Normalization', 'probability');
B = arrayfun(@(x) x.T(x.indWhenCellDies), Tracks1,'UniformOutput', false);

A = 1000*cat(1,A{cellfun(@(x) ~isempty(x), A)});
B = cat(1,B{cellfun(@(x) ~isempty(x), B)});
%A=A(B>20);
%B=B(B>20);

%scatter(B/3,asinh(A),5,tzeva(j,:))
%hold on
[xx,yy, stdXb, stdYb, steXb, steYb,xt] = BinData_v1(B,(A),20)
%  J=cellfun(@numel,xt)>0;
%   h1=ploterr(xx(J)/3,yy(J),steXb(J)/3,steYb(J),'o');
%   h1(1).Color=tzeva(j,:);
%   h1(2).Color=tzeva(j,:);
%   h1(3).Color=tzeva(j,:);

 h1=plot(xx/3,yy,'-o');
 h1(1).Color=tzeva(j,:);
 h1(1).MarkerSize=5;
 h1(1).MarkerFaceColor=tzeva(j,:);
ylabel('Viral load at death(a.u.)')
xlabel('Time of death (h)')
set(gca,'YTickLabel',[0:50:500],'YTick',([0:50:500]))
hold on
shg
end
set(gca,'ylim',[0,250],'xlim',[0,48])
%%
figname = 'ViralLoadAtDeathMeans'
set(gcf, 'PaperPositionMode','auto','color','w','InvertHardcopy','off')
print(gcf,'-depsc','-r600',[savpath '/epss/' figname]);
print(gcf,'-dpng','-r600',[savpath '/pngs/' figname]);
%%
%%
figure('color','w','Position',[100,100, 450, 450])

ax = axes('Position', [0.24, 0.18, 0.7, 0.7])

tzeva = cmTNF(10);
hs=[];
for j=1:10
    jj=[j, j+10, j+20]+30
pos = R.PosNames(jj);
Tracks1 = R.getTracks(pos)
Tracks1 = Tracks1(arrayfun(@(x) any(x.Infected),Tracks1))

A = arrayfun(@(x) x.VirusTrack(x.indWhenCellDies), Tracks1,'UniformOutput', false);
%h = histogram(cat(1,A{cellfun(@(x) ~isempty(x), A)}),50,'Normalization', 'probability');
%B = arrayfun(@(x) x.T(x.indWhenCellDies), Tracks1,'UniformOutput', false);

A = 1000*cat(1,A{cellfun(@(x) ~isempty(x), A)});
%B = cat(1,B{cellfun(@(x) ~isempty(x), B)});
%A=A(B>20);
%B=B(B>20);

%scatter(B/3,asinh(A),5,tzeva(j,:))
%hold on
%[xx,yy, stdXb, stdYb, steXb, steYb,xt] = BinData_v1(B,(A),20)
%  J=cellfun(@numel,xt)>0;
%   h1=ploterr(xx(J)/3,yy(J),steXb(J)/3,steYb(J),'o');
%   h1(1).Color=tzeva(j,:);
%   h1(2).Color=tzeva(j,:);
%   h1(3).Color=tzeva(j,:);

 [h1, edges] =histcounts(A, 0:15:450);
 centers = 2\(edges(2:end)+edges(1:end-1));
 h=plot(centers, h1./sum(h1))
 h(1).Color=tzeva(j,:);
 %h1(1).MarkerSize=5;
 %h1(1).FaceColor=tzeva(j,:);
ylabel('probability')
xlabel('Viral load at death (a.u.)')
set(gca,'XTickLabel',[0:100:500],'XTick',([0:100:500]))
hold on
shg
end
set(gca,'xlim',[0,450])
%%
figname = 'hists_viralloadatdeath'
set(gcf, 'PaperPositionMode','auto','color','w','InvertHardcopy','off')
print(gcf,'-depsc','-r600',[savpath '/epss/' figname]);
print(gcf,'-dpng','-r600',[savpath '/pngs/' figname]);
%%
figure('color','w','Position',[100,100, 450, 450])

tzeva = (cmTNF(10));

ax = axes('Position', [0.24, 0.18, 0.7, 0.7])
ps=[]
Errs=[]

YYs = []
for j=1:10
    jj=[j, j+10, j+20]+30
pos = R.PosNames(jj);
Tracks1 = R.getTracks(pos)
Tracks1 = Tracks1(arrayfun(@(x) any(x.Infected),Tracks1))

VT = arrayfun(@(x) x.VirusTrack(find((~x.Dead))),Tracks1,'UniformOutput',false);%(find((x.Dead)))
T = arrayfun(@(x) x.T(find((~x.Dead))),Tracks1,'UniformOutput',false);%(find((x.Dead)))

VT = 1000*cat(2,VT{cellfun(@(x) ~isempty(x), VT)});
T = cat(2,T{cellfun(@(x) ~isempty(x), T)});


[xx,yy, stdXb, stdYb, steXb, steYb,xt] = BinData_v1(T,VT,30)
 J=cellfun(@numel,xt)>0;
  h1=ploterr(xx(J)/3,yy(J),steXb(J)/3,steYb(J),'o');
  h1(1).Color=tzeva(j,:);
  h1(2).Color=tzeva(j,:);
  h1(3).Color=tzeva(j,:);
  hold on
  
  J = 6:11;
  [p, S] = polyfit(xx(J)/3,yy(J),1);
  
  [y_fit,delta] = polyval(p,xx(J)/3,S);
  plot(xx(J)/3,y_fit,'Color',tzeva(j,:));
  ps=[ps;p];
  Err = polyparci(p,S, .95);
  Errs=[Errs; Err(:,1)']
  clear range; ErrR = range(Err)/2;
  %text(1, 10, ['Slope=' num2str(p(1),2) 'u/h'],'Color','k');
  ylabel('Average viral load(a.u.)')
xlabel('time (h)')

    YYs = [YYs, yy(1:12)', steYb(1:12)', nan(12,1)];
end
set(gca,'ylim',[40,250],'xlim',[0,20])

ax2 = axes('Position', [0.26, 0.61, 0.25, 0.25])
LinearRange=0.1
h1= ploterr(asinh(cTNF/LinearRange),ps(:,1),[],{Errs(:,1), Errs(:,2)})

BarTicklength(h1(2),0)
for i=1:length(h1) 
    h1(i).Color = 'k';
    h1(i).MarkerSize=7
    h1(i).MarkerFaceColor = 'k';
    
end

set(gca,'ylim',[0,20],'xlim',[0,asinh(100/LinearRange)])

ax2.YAxisLocation='right'

a1 = sort(- (1:9)'*10.^(1:2));
a2 =(1:9)'*10.^(1:5);
a = [sort(a1(:)); 0; a2(:)]/1000;
emptyStr = {''; ''; '';''; ''; '';''; ''};
XTickLabel = [emptyStr; ''; emptyStr; {''}; '0'; {''}; emptyStr; ' 0.1'; emptyStr; ' 1'; emptyStr; ' 10'; emptyStr; '  100'; emptyStr];
XTicks = asinh(a/LinearRange);
set(gca, 'xtick', XTicks, 'xticklabel', XTickLabel,'FontSize', 12)
xl = xlabel('[TNF] ng/ml')
xl.FontSize=14
yl = ylabel({'Viral growth rate' '(a.u./h)'})

%%
figname = 'Fig3S_ViralGrowthForDifferentTNF'
set(gcf, 'PaperPositionMode','auto','color','w','InvertHardcopy','off')
print(gcf,'-depsc','-r600',[savpath '/epss/' figname]);
print(gcf,'-dpng','-r600',[savpath '/pngs/' figname]);
%% Plot both uninfected and infected cell death rates

allParams = R.getData('ExpFitParam');
allParams = cat(1,allParams{:})
allParams = allParams(:,2);
allParams = reshape(allParams,[],6);
figure('color','w','Position',[100,100, 450, 450])

ax = axes('Position', [0.2, 0.15, 0.7, 0.7])
LinearRange = [0.1]
LinearRangeY = [0.05]

infDeathRatesExp = mean(allParams(:,4:6),2);
uninfDeathRatesExp = mean(allParams(:,1:3),2);

hold on
x = cTNF'
y = asinh(infDeathRatesExp/LinearRangeY); %Infected
y2 = asinh(uninfDeathRatesExp/LinearRangeY); %Uninfected
%Error calculation: asinh' = 1/sqrt(1+x^2)
sigmay = std(allParams(:,4:6),[],2)./sqrt(mean(allParams(:,4:6)./LinearRangeY,2).^2+1)./LinearRangeY;
sigmay2 = std(allParams(:,1:3),[],2)./sqrt(mean(allParams(:,1:3)./LinearRangeY,2).^2+1)./LinearRangeY;
h = errorbar(asinh(cTNF/LinearRange)',y', sigmay','-o', 'linewidth', 2, 'markersize', 7,'color',virusColor, 'markerfacecolor', virusColor,'markeredgecolor','none');
h2 = errorbar(asinh(cTNF/LinearRange)',y2',sigmay2','-o', 'linewidth', 2, 'markersize', 7,'color',liveColor, 'markerfacecolor', liveColor,'markeredgecolor','none');

TNF = [0 fliplr(100*(1/sqrt(2)).^[0:17])];

infDeathRate = max(interp1(cTNF,filtfilt([1,1,1],3,infDeathRatesExp),TNF),0);
basalDeathRate = max(interp1(cTNF,filtfilt([1,1,1],3,uninfDeathRatesExp),TNF),0);

%h3 = plot(asinh(TNF/LinearRange)', asinh(infDeathRate/LinearRangeY),'-o','color',virusColor, 'linewidth', 2, 'markersize', 7)
%h4 = plot(asinh(TNF/LinearRange)', asinh(basalDeathRate/LinearRangeY),'-o','color',liveColor, 'linewidth', 2, 'markersize', 7)





a1 = sort(- (1:9)'*10.^(1:2));
a2 =(1:9)'*10.^(1:5);
a = [sort(a1(:)); 0; a2(:)]/1000;
emptyStr = {''; ''; '';''; ''; '';''; ''};
XTickLabel = [emptyStr; ''; emptyStr; {''}; '0'; {''}; emptyStr; ' 0.1'; emptyStr; ' 1'; emptyStr; ' 10'; emptyStr; '    100'; emptyStr];
XTicks = asinh(a/LinearRange);







a1 = sort(- (1:9)'*10.^(1:2));
a2 =(1:9)'*10.^(1:5);
a = [sort(a1(:)); 0; a2(:)]/10000;
emptyStr = {''; ''; '';''; ''; '';''; ''};
YTickLabel = [emptyStr; ' '; emptyStr; {''}; '0'; {''}; emptyStr; ' 0.01'; emptyStr; ' .1'; emptyStr; ' 1'; emptyStr; '    10'; emptyStr];

YTicks = asinh(a/LinearRangeY);

set(gca, 'fontsize', 17, 'xtick', XTicks, 'xticklabel', XTickLabel,  'ytick', YTicks, 'yticklabel', YTickLabel, 'xlim', [-0.5 8])
hl = legend('Infected cells', 'Uninfected cells')
hl.Box = 'off'
hl.Position(1) = hl.Position(1)-0.3;
ylabel('Net Death rate:    \beta-k_{p}')
xlabel('[TNF] ng/ml')
box on
%title('Uninfected and Infected death rates')


%%
figname = 'Fig3S_NetDeathRates'
set(gcf, 'PaperPositionMode','auto','color','w','InvertHardcopy','off')
print(gcf,'-depsc','-r600',[savpath '/epss/' figname]);
print(gcf,'-dpng','-r600',[savpath '/pngs/' figname]);





%% This is how you calculate rates:
infDeathRatesExp = mean(allParams(:,4:6),2);
uninfDeathRatesExp = mean(allParams(:,1:3),2);

cTNF = ([100*(1/2).^[0:8], 0]);

infDeathRate = interp1(cTNF,infDeathRatesExp,TNF);
basalDeathRate = interp1(cTNF,uninfDeathRatesExp,TNF);

%% Make triangular grid
A = triangleGrid([-10 -10 10 10], [0,0], 1);
nCells = size(A,1);
%find all nearest neighnors and create sparse matrix
[IDX, D] = rangesearch(A, A,1.01);
IDX = cellfun(@(x,y) x(y>0), IDX,D,'uniformoutput',false);
pointsToMatch = arrayfun(@(x,y) repmat(x,y,1),1:numel(IDX), cellfun(@numel,IDX)','UniformOutput', false);
nnMatrix = sparse(cat(1,pointsToMatch{:}), cat(2,IDX{:})',1);

pointSize = 40;
lineSize = 2;
%% Init cells, small frac infected
%fracInf = 0.01;
%infected = rand(nCells,1)<fracInf;
infected = zeros(nCells,1);
infected(116)=1;
clear h

x= [~infected, infected, zeros(nCells,1), zeros(nCells,1)]';
x = x(:);

figure('color','w')
ax1 = axes('position',[0,0,1,1]);
%voronoi(A(:,1), A(:,2)); hold on;
h(1) = scatter(A(logical(x(1:4:end)),1),A(logical(x(1:4:end)),2),pointSize,liveColor,'filled')
h(1).MarkerEdgeColor = liveColor;
hold on;
h(2) = scatter(A(logical(x(2:4:end)),1),A(logical(x(2:4:end)),2),pointSize,virusColor,'filled')
h(2).MarkerEdgeColor = virusColor;
h(3) = scatter(A(logical(x(3:4:end)),1),A(logical(x(3:4:end)),2),pointSize,deathColor,'filled')
h(3).MarkerEdgeColor = deathColor;
h(4) = scatter(A(logical(x(4:4:end)),1),A(logical(x(4:4:end)),2),pointSize,virusColor,'filled')
h(4).MarkerEdgeColor = deathColor;
axis equal
set(ax1,'xcolor','w','ycolor','w','xlim', [-10, 10],'ylim', [-10, 10])
shg
%hl = legend(h,'Healthy','Infected','False positive','Dead post infection')
%hl.Location = 'northeastoutside';
%hl.Box = 'off';

figure(gcf)
%set(gcf, 'PaperPositionMode','auto','color','w','InvertHardcopy','off')
%print(gcf,'-dpng','-r300',['/Users/Alonyan/Desktop/' 'Initial']);



%% Make Figure model spread
% Load high MOI
BaseStr = regexprep([char(ispc.*'Z:\Images2019\') char(isunix.*'/bigstore/Images2019/')],char(0),'');
Usr = 'Jen';
Project = 'NFkBDynamics';
Dataset = 'TNFTitr_HighMOI_Dec122019_2019Dec12';
acquisition = 2;

% Get MD of raw data
acqname = ['acq_' num2str(acquisition)];
fpath = [BaseStr Usr filesep Project filesep Dataset filesep acqname];
MD=Metadata(fpath,[],1);
Wells = unique(MD.getSpecificMetadata('Position'));
frames = unique(cell2mat(MD.getSpecificMetadata('frame')));
dt = 1/3;%h

cTNF = ([100*(1/2).^[0:8], 0]);
% Load Results
R = MultiPositionSingleCellVirusResults(fpath)
Wells = R.PosNames;
frames = R.Frames;
t = [R.Frames-1]'.*dt;
%
TNF = R.getData('modelTNF')
Stats = R.getData('modelSTATs')
xs = R.getData('modelResultsxs')

A = triangleGrid([-7 -7 7 7], [0,0], 1);
nCells = size(A,1);

B = triangleGrid([-9 -9 9 9], [0,0], 1);

% 
% figure('color','w','Position',[100, 100, 1200, 400])
% rang = [1 9 19];
% for i=1:numel(rang)
%     xOut = mean(xs{rang(i)});
%     %ax1 = subplot(2,10,i)
%     ax1 = axes('position',[0.33*(i-1),0.1,0.3,0.8]);
%     %h = scatter(A(:,1),A(:,2),150,'k' );%xOut(2:4:end)
%     %h.MarkerEdgeAlpha=0.1
%     h1 = voronoi(B(:,1),B(:,2));
%     h1(1).Color = [0 0 0 0.2];
%     h1(2).Color = [0 0 0 0.2];
%     hold on
%     vCM = makeColorMap([1,1,1],liveColor,1001);
%     h = scatter(A(:,1),A(:,2),270*(xOut(2:4:end)+xOut(1:4:end)+xOut(4:4:end))+0.0001,vCM(interp1(0:0.001:1,0:1000,(xOut(1:4:end))./max(xOut(1:4:end)+0.00001),'nearest')+1,:) ,'filled');%xOut(2:4:end)
%     h.MarkerFaceAlpha = 1
%     hold on;
%     dCM = makeColorMap([1,1,1], [0 0 0],1001);
%     h = scatter(A(:,1),A(:,2),270*(xOut(2:4:end)+xOut(4:4:end))+0.0001,dCM(interp1(0:0.001:1,0:1000,(xOut(4:4:end))./max(xOut(4:4:end)+0.00001),'nearest')+1,:) ,'filled');%xOut(2:4:end)
%     h.MarkerFaceAlpha = 1
% 
%     dCM = makeColorMap([1,1,1], virusColor,1001);
%     h = scatter(A(:,1),A(:,2),270*xOut(2:4:end)+0.0001,dCM(interp1(0:0.001:1,0:1000,xOut(2:4:end)./max(xOut(2:4:end)+0.00001),'nearest')+1,:) ,'filled');%xOut(2:4:end)
%     h.MarkerFaceAlpha = 1
%     
%     
%     axis equal
%     set(ax1,'xcolor','none','ycolor','none','XLim',[-10,10],'YLim',[-10,10]);
%     %an = annotation('textbox',[0.03+0.33*(i-1),0.88,0.3,0.1],'String',{['TNF=' num2str(TNF(rang(i))) 'ng/ml']},'LineStyle','none');
% end

%%
figure('color','w','Position',[100, 100, 1200, 400])
rang = [1 3 5  19];
for i=1:numel(rang)
    xOut = mean(xs{rang(i)});
    %ax1 = subplot(2,10,i)
    ax1 = axes('position',[0.31*(i-1),0.1,0.3,0.8]);
    %h = scatter(A(:,1),A(:,2),150,'k' );%xOut(2:4:end)
    %h.MarkerEdgeAlpha=0.1
    h1 = voronoi(B(:,1),B(:,2));
    h1(1).Color = 'none';
    h1(2).Color = [0 0 0 0.2];
    hold on
    vCM = makeColorMap([1,1,1],virusColor,1001);
    h = scatter(A(:,1),A(:,2),220*(xOut(2:4:end)+xOut(4:4:end)).^2+0.000001,vCM(interp1(0:0.001:1,0:1000,(xOut(2:4:end))./max(xOut(2:4:end)+0.00001),'nearest')+1,:) ,'filled');%xOut(2:4:end)
    h.MarkerFaceAlpha = 1
    hold on;
    dCM = makeColorMap([1,1,1], [0 0 0],1001);
    h = scatter(A(:,1),A(:,2),220*xOut(4:4:end).^2+0.000001,dCM(interp1(0:0.001:1,0:1000,xOut(4:4:end)./max(xOut(4:4:end)+0.00001),'nearest')+1,:) ,'filled');%xOut(2:4:end)
    h.MarkerFaceAlpha = 1
% 
%     dCM = makeColorMap([1,1,1], liveColor,1001);
%     h = scatter(A(:,1),A(:,2),300*xOut(1:4:end)+0.0001,dCM(interp1(0:0.001:1,0:1000,xOut(1:4:end)./max(xOut(1:4:end)+0.00001),'nearest')+1,:) ,'filled');%xOut(2:4:end)
%     h.MarkerFaceAlpha = 0.5
    
    
    axis equal
    set(ax1,'xcolor','none','ycolor','none','XLim',[-10,10],'YLim',[-10,10]);
    %an = annotation('textbox',[0.03+0.33*(i-1),0.88,0.3,0.1],'String',{['TNF=' num2str(TNF(rang(i))) 'ng/ml']},'LineStyle','none');
end

%%
figure('color','w','Position',[100, 100, 1200, 400])
rang = [1 3 19];
for i=1:numel(rang)
    xOut = mean(xs{rang(i)});
    %ax1 = subplot(2,10,i)
    ax1 = axes('position',[0.31*(i-1),0.1,0.3,0.8]);
    %h = scatter(A(:,1),A(:,2),150,'k' );%xOut(2:4:end)
    %h.MarkerEdgeAlpha=0.1
    h1 = voronoi(B(:,1),B(:,2));
    h1(1).Color = 'none';
    h1(2).Color = [0 0 0 0.2];
    hold on
    vCM = makeColorMap([1,1,1],virusColor,1001);
    h = scatter(A(:,1),A(:,2),220*(xOut(2:4:end)+xOut(4:4:end)).^2+0.000001,vCM(interp1(0:0.001:1,0:1000,(xOut(2:4:end))./max(xOut(2:4:end)+0.00001),'nearest')+1,:) ,'filled');%xOut(2:4:end)
    h.MarkerFaceAlpha = 1
    hold on;
    dCM = makeColorMap([1,1,1], [0 0 0],1001);
    h = scatter(A(:,1),A(:,2),220*xOut(4:4:end).^2+0.000001,dCM(interp1(0:0.001:1,0:1000,xOut(4:4:end)./max(xOut(4:4:end)+0.00001),'nearest')+1,:) ,'filled');%xOut(2:4:end)
    h.MarkerFaceAlpha = 1
% 
%     dCM = makeColorMap([1,1,1], liveColor,1001);
%     h = scatter(A(:,1),A(:,2),300*xOut(1:4:end)+0.0001,dCM(interp1(0:0.001:1,0:1000,xOut(1:4:end)./max(xOut(1:4:end)+0.00001),'nearest')+1,:) ,'filled');%xOut(2:4:end)
%     h.MarkerFaceAlpha = 0.5
    
    
    axis equal
    set(ax1,'xcolor','none','ycolor','none','XLim',[-7,7],'YLim',[-7,7]);
    %an = annotation('textbox',[0.03+0.33*(i-1),0.88,0.3,0.1],'String',{['TNF=' num2str(TNF(rang(i))) 'ng/ml']},'LineStyle','none');
end

%%
figname = 'Fig3_Model_Averages'
set(gcf, 'PaperPositionMode','auto','color','w','InvertHardcopy','off','renderer', 'painters')
print(gcf,'-depsc','-r600',[savpath '/epss/' figname]);
print(gcf,'-dpng','-r600',[savpath '/pngs/' figname]);

%%
figure('color','w','Position',[100, 100,900, 600])
rang = [1 3 6];
for i=1:numel(rang)
    xOut = mean(xs{rang(i)});
    %ax1 = subplot(2,10,i)
    ax1 = axes('position',[0.03+0.31*(i-1),0.52,0.28,0.4]);
    %h = scatter(A(:,1),A(:,2),150,'k' );%xOut(2:4:end)
    %h.MarkerEdgeAlpha=0.1
    h1 = voronoi(B(:,1),B(:,2));
    h1(1).Color = 'none';
    h1(2).Color = [0 0 0 0.2];
    hold on
    vCM = makeColorMap([1,1,1],virusColor,1001);
    h = scatter(A(:,1),A(:,2),150*(xOut(2:4:end)+xOut(4:4:end)).^2+0.000001,vCM(interp1(0:0.001:1,0:1000,(xOut(2:4:end))./max(xOut(2:4:end)+0.00001),'nearest')+1,:) ,'filled');%xOut(2:4:end)
    h.MarkerFaceAlpha = 1
    hold on;
    dCM = makeColorMap([1,1,1], [0 0 0],1001);
    h = scatter(A(:,1),A(:,2),150*xOut(4:4:end).^2+0.000001,dCM(interp1(0:0.001:1,0:1000,xOut(4:4:end)./max(xOut(4:4:end)+0.00001),'nearest')+1,:) ,'filled');%xOut(2:4:end)
    h.MarkerFaceAlpha = 1
% 
     dCM = makeColorMap([1,1,1], liveColor,1001);
     h = scatter(A(:,1),A(:,2),150*xOut(1:4:end)+0.0001,dCM(interp1(0:0.001:1,0:1000,xOut(1:4:end)./max(xOut(1:4:end)+0.00001),'nearest')+1,:) ,'filled');%xOut(2:4:end)
     h.MarkerFaceAlpha = 0.5
    
    
    axis equal
    set(ax1,'xcolor','none','ycolor','none','XLim',[-7,7],'YLim',[-7,7]);
    %an = annotation('textbox',[0.03+0.33*(i-1),0.88,0.3,0.1],'String',{['TNF=' num2str(TNF(rang(i))) 'ng/ml']},'LineStyle','none');
end

rang = [9 13 19];
for i=1:numel(rang)
    xOut = mean(xs{rang(i)});
    %ax1 = subplot(2,10,i)
    ax1 = axes('position',[0.03+0.31*(i-1),0.03,0.28,0.4]);
    %h = scatter(A(:,1),A(:,2),150,'k' );%xOut(2:4:end)
    %h.MarkerEdgeAlpha=0.1
    h1 = voronoi(B(:,1),B(:,2));
    h1(1).Color = 'none';
    h1(2).Color = [0 0 0 0.2];
    hold on
    vCM = makeColorMap([1,1,1],virusColor,1001);
    h = scatter(A(:,1),A(:,2),150*(xOut(2:4:end)+xOut(4:4:end)).^2+0.000001,vCM(interp1(0:0.001:1,0:1000,(xOut(2:4:end))./max(xOut(2:4:end)+0.00001),'nearest')+1,:) ,'filled');%xOut(2:4:end)
    h.MarkerFaceAlpha = 1
    hold on;
    dCM = makeColorMap([1,1,1], [0 0 0],1001);
    h = scatter(A(:,1),A(:,2),150*xOut(4:4:end).^2+0.000001,dCM(interp1(0:0.001:1,0:1000,xOut(4:4:end)./max(xOut(4:4:end)+0.00001),'nearest')+1,:) ,'filled');%xOut(2:4:end)
    h.MarkerFaceAlpha = 1
% 
     dCM = makeColorMap([1,1,1], liveColor,1001);
     h = scatter(A(:,1),A(:,2),150*xOut(1:4:end)+0.0001,dCM(interp1(0:0.001:1,0:1000,xOut(1:4:end)./max(xOut(1:4:end)+0.00001),'nearest')+1,:) ,'filled');%xOut(2:4:end)
     h.MarkerFaceAlpha = 0.5
    
    axis equal
    set(ax1,'xcolor','none','ycolor','none','XLim',[-7,7],'YLim',[-7,7]);
    %an = annotation('textbox',[0.03+0.33*(i-1),0.88,0.3,0.1],'String',{['TNF=' num2str(TNF(rang(i))) 'ng/ml']},'LineStyle','none');
end
%%
figname = 'Fig3_Model_Averages_6panels_withLive'
set(gcf, 'PaperPositionMode','auto','color','w','InvertHardcopy','off','renderer', 'painters')
print(gcf,'-depsc','-r600',[savpath '/epss/' figname]);
print(gcf,'-dpng','-r600',[savpath '/pngs/' figname]);

%% Plot Initial Conditions
close all
clear h
VGR = 1/10;%viral growth rate
VI = 1/2; %viral infectivity

TNF = 0;
infDeathRate = HillFunction(HillFuncBeta,TNF)
basalDeathRate = Rectifier(RectifierBETA,infDeathRate)


xOut = x% = SSA_forSpatialSIDGrids(x,nnMatrix,A,VGR,VI,basalDeathRate,infDeathRate);

sum(reshape(xOut,4,[]),2)./sum(sum(reshape(xOut,4,[]),2))

%
figure('color','w','Position',[100,100, 450, 450],'renderer', 'painters')
ax1 = axes('position',[0.02,0.1,0.5,0.5]);
h(1) = scatter(A(logical(x(1:4:end)),1),A(logical(x(1:4:end)),2),pointSize,liveColor,'filled')
h(1).MarkerEdgeColor = liveColor;
hold on;
h(2) = scatter(A(logical(x(2:4:end)),1),A(logical(x(2:4:end)),2),pointSize,virusColor,'filled')
h(2).MarkerEdgeColor = virusColor;
h(3) = scatter(A(logical(x(3:4:end)),1),A(logical(x(3:4:end)),2),pointSize,deathColor,'filled')
h(3).MarkerEdgeColor = deathColor;
h(4) = scatter(A(logical(x(4:4:end)),1),A(logical(x(4:4:end)),2),pointSize,deathColor,'filled')
h(4).MarkerEdgeColor = virusColor;
for ii=1:4
    h(ii).LineWidth = lineSize;
end


axis equal
set(ax1,'xcolor','w','ycolor','w');
title({['Initial Conditions']})
%an = annotation('textbox',[0.7, 0.1, 0.14, 0.4],'String',{['TNF=' num2str(TNF) 'ng/ml']},'LineStyle','none')
shg
hl = legend('   Healthy','   Infected','   Dead Bystander','   Dead post infection');
hl.Location = 'northeastoutside';
hl.Position = [0.63    0.4    0.3    0.2046];
hl.Box = 'off';




%%
figname = 'Fig2_Model_Initial'
set(gcf, 'PaperPositionMode','auto','color','w','InvertHardcopy','off')
print(gcf,'-depsc','-r600',[savpath '/epss/' figname]);
print(gcf,'-dpng','-r600',[savpath '/pngs/' figname]);





%% Run model Examples
close all
VGR = 1/10;%viral growth rate
VI = 1/2; %viral infectivity
clear h

TNF = 0;
infDeathRate = HillFunction(HillFuncBeta,TNF)
basalDeathRate = Rectifier(RectifierBETA,infDeathRate)


xOut = SSA_forSpatialSIDGrids(x,nnMatrix,A,VGR,VI,basalDeathRate,infDeathRate);

sum(reshape(xOut,4,[]),2)./sum(sum(reshape(xOut,4,[]),2))

%
figure('color','w','Position',[100,100, 450, 450],'renderer','painters')
ax1 = axes('position',[0.02,0.1,0.5,0.5]);
h(1) = scatter(A(logical(xOut(1:4:end)),1),A(logical(xOut(1:4:end)),2),pointSize,liveColor,'filled')
h(1).MarkerEdgeColor = liveColor;
hold on;
h(2) = scatter(A(logical(xOut(2:4:end)),1),A(logical(xOut(2:4:end)),2),pointSize,virusColor,'filled')
h(2).MarkerEdgeColor = virusColor;
h(3) = scatter(A(logical(xOut(3:4:end)),1),A(logical(xOut(3:4:end)),2),pointSize,deathColor,'filled')
h(3).MarkerEdgeColor = deathColor;
h(4) = scatter(A(logical(xOut(4:4:end)),1),A(logical(xOut(4:4:end)),2),pointSize,deathColor,'filled')
h(4).MarkerEdgeColor = virusColor;
for ii=1:4
    h(ii).LineWidth = lineSize;
end

axis equal
set(ax1,'xcolor','w','ycolor','w');
title({['TNF\alpha = ' num2str(TNF) ' ng/ml']})
%an = annotation('textbox',[0.7, 0.1, 0.14, 0.4],'String',{['TNF=' num2str(TNF) 'ng/ml']},'LineStyle','none')
shg
hl = legend('   Healthy','   Infected','   Dead Bystander','   Dead post infection');
hl.Location = 'northeastoutside';
hl.Position = [0.63    0.4    0.3    0.2046];
hl.Box = 'off';

%%
figname = 'Fig2_Model_0TNF'
set(gcf, 'PaperPositionMode','auto','color','w','InvertHardcopy','off')
print(gcf,'-depsc','-r600',[savpath '/epss/' figname]);
print(gcf,'-dpng','-r600',[savpath '/pngs/' figname]);


%% Run model Examples
close all
VGR = 1/10;%viral growth rate
VI = 1/2; %viral infectivity
clear h

TNF = 5;
infDeathRate = HillFunction(HillFuncBeta,TNF)
basalDeathRate = Rectifier(RectifierBETA,infDeathRate)


xOut = SSA_forSpatialSIDGrids(x,nnMatrix,A,VGR,VI,basalDeathRate,infDeathRate);

sum(reshape(xOut,4,[]),2)./sum(sum(reshape(xOut,4,[]),2))

%
figure('color','w','Position',[100,100, 450, 450],'renderer','painters')
ax1 = axes('position',[0.02,0.1,0.5,0.5]);
h(1) = scatter(A(logical(xOut(1:4:end)),1),A(logical(xOut(1:4:end)),2),pointSize,liveColor,'filled')
h(1).MarkerEdgeColor = liveColor;
hold on;
h(2) = scatter(A(logical(xOut(2:4:end)),1),A(logical(xOut(2:4:end)),2),pointSize,virusColor,'filled')
h(2).MarkerEdgeColor = virusColor;
h(3) = scatter(A(logical(xOut(3:4:end)),1),A(logical(xOut(3:4:end)),2),pointSize,deathColor,'filled')
h(3).MarkerEdgeColor = deathColor;
h(4) = scatter(A(logical(xOut(4:4:end)),1),A(logical(xOut(4:4:end)),2),pointSize,deathColor,'filled')
h(4).MarkerEdgeColor = virusColor;
for ii=1:4
    h(ii).LineWidth = lineSize;
end

axis equal
set(ax1,'xcolor','w','ycolor','w');
title({['TNF\alpha = ' num2str(TNF) ' ng/ml']})
%an = annotation('textbox',[0.7, 0.1, 0.14, 0.4],'String',{['TNF=' num2str(TNF) 'ng/ml']},'LineStyle','none')
shg
hl = legend('   Healthy','   Infected','   Dead Bystander','   Dead post infection');
hl.Location = 'northeastoutside';
hl.Position = [0.63    0.4    0.3    0.2046];
hl.Box = 'off';

%%
figname = 'Fig2_Model_5TNF'
set(gcf, 'PaperPositionMode','auto','color','w','InvertHardcopy','off')
print(gcf,'-depsc','-r600',[savpath '/epss/' figname]);
print(gcf,'-dpng','-r600',[savpath '/pngs/' figname]);

%% Run model Examples
close all
VGR = 1/10;%viral growth rate
VI = 1/2; %viral infectivity
clear h

TNF = 100;
infDeathRate = HillFunction(HillFuncBeta,TNF)
basalDeathRate = Rectifier(RectifierBETA,infDeathRate)


xOut = SSA_forSpatialSIDGrids(x,nnMatrix,A,VGR,VI,basalDeathRate,infDeathRate);

sum(reshape(xOut,4,[]),2)./sum(sum(reshape(xOut,4,[]),2))

%
figure('color','w','Position',[100,100, 450, 450],'renderer','painters')
ax1 = axes('position',[0.02,0.1,0.5,0.5]);
h(1) = scatter(A(logical(xOut(1:4:end)),1),A(logical(xOut(1:4:end)),2),pointSize,liveColor,'filled')
h(1).MarkerEdgeColor = liveColor;
hold on;
h(2) = scatter(A(logical(xOut(2:4:end)),1),A(logical(xOut(2:4:end)),2),pointSize,virusColor,'filled')
h(2).MarkerEdgeColor = virusColor;
h(3) = scatter(A(logical(xOut(3:4:end)),1),A(logical(xOut(3:4:end)),2),pointSize,deathColor,'filled')
h(3).MarkerEdgeColor = deathColor;
h(4) = scatter(A(logical(xOut(4:4:end)),1),A(logical(xOut(4:4:end)),2),pointSize,deathColor,'filled')
h(4).MarkerEdgeColor = virusColor;
for ii=1:4
    h(ii).LineWidth = lineSize;
end

axis equal
set(ax1,'xcolor','w','ycolor','w');
title({['TNF\alpha = ' num2str(TNF) ' ng/ml']})
%an = annotation('textbox',[0.7, 0.1, 0.14, 0.4],'String',{['TNF=' num2str(TNF) 'ng/ml']},'LineStyle','none')
shg
hl = legend('   Healthy','   Infected','   Dead Bystander','   Dead post infection');
hl.Location = 'northeastoutside';
hl.Position = [0.63    0.4    0.3    0.2046];
hl.Box = 'off';

%%
figname = 'Fig2_Model_100TNF'
set(gcf, 'PaperPositionMode','auto','color','w','InvertHardcopy','off')
print(gcf,'-depsc','-r600',[savpath '/epss/' figname]);
print(gcf,'-dpng','-r600',[savpath '/pngs/' figname]);




%% Movies of model stats
close all
VGR = 1/10;%viral growth rate
VI = 1/2; %viral infectivity

TNF = 0;



movieName = 'modelExample0TNF'
xOut = SSA_forSpatialSIDGrids(x,nnMatrix,A,VGR,VI,basalDeathRate,infDeathRate,'plot', true, 'TNF', TNF, 'savpath', [savpath '/movies/' movieName]);


%% Movies of model stats
close all
VGR = 1/10;%viral growth rate
VI = 1/2; %viral infectivity

TNF = 5;



movieName = 'modelExample5TNF'
xOut = SSA_forSpatialSIDGrids(x,nnMatrix,A,VGR,VI,basalDeathRate,infDeathRate,'plot', true, 'TNF', TNF, 'savpath', [savpath '/movies/' movieName]);

%% Movies of model stats
close all
VGR = 1/10;%viral growth rate
VI = 1/2; %viral infectivity

TNF = 100;



movieName = 'modelExample100TNF'
xOut = SSA_forSpatialSIDGrids(x,nnMatrix,A,VGR,VI,basalDeathRate,infDeathRate,'plot', true, 'TNF', TNF, 'savpath', [savpath '/movies/' movieName]);




%% Get model stats
TNF = R.getData('modelTNF')
Stats = R.getData('modelSTATs')
StatsNoDeath = R.getData('modelSTATsNoDeath')


%% Plot Model statistics
close all
LinearRange = 0.4

figure('color','w','Position',[100,100, 450, 450])
ax = axes('Position', [0.15, 0.14, 0.7, 0.7])

a = cellfun(@(x) x(:,1), Stats,'uniformOutput',false)
yyaxis(ax,'left')
plot(asinh(TNF/LinearRange),median(cat(2,a{:})),'Color',liveColor);
set(gca,'YLim',[-1,100]/100,'YColor',liveColor)
ylabel('Percent Healthy (TN)')

b = cellfun(@(x) x(:,2), Stats,'uniformOutput',false)
yyaxis(ax,'right')
plot(asinh(TNF/LinearRange), median(cat(2,b{:})),'Color',virusColor);
set(gca,'YLim',[-1,100]/100,'YColor',virusColor)
ax.XLim=asinh([0 100]/LinearRange)
ax.YLim=[-1 100]/100
ylabel('Percent Infected (FN)')

xl = xlabel('[TNF] ng/ml')

title('Model')

hold on
%
yyaxis(ax,'left')
h = boxplot(cat(2,a{:}), 'Positions',asinh(TNF/LinearRange) , 'PlotStyle','compact', 'MedianStyle','line','Symbol','.','Widths',0.3, 'Labels',{'' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' },'Color', liveColor);
shg
set(gca,'YLim',[-1,100]/100)
ax.YTick = 0:0.2:1;
ax.YTickLabel = 100*[0:0.2:1];

yyaxis(ax,'right')
h = boxplot(cat(2,b{:}), 'Positions',asinh(TNF/LinearRange), 'PlotStyle','compact', 'MedianStyle','line','Symbol','.','Widths',0.3, 'Labels',  {'' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' },'Color', virusColor);
shg
set(gca,'YLim',[-1,100]/100)


a1 = sort(- (1:9)'*10.^(1:2));
a2 =(1:9)'*10.^(1:5);
a = [sort(a1(:)); 0; a2(:)]/1000;
emptyStr = {''; ''; '';''; ''; '';''; ''};
ax.XTickLabel = [emptyStr; ''; emptyStr; {''}; '0'; {''}; emptyStr; '  .1'; emptyStr; ' 1'; emptyStr; ' 10'; emptyStr; '    100'; emptyStr];
ax.XTick = asinh(a/LinearRange);

ax.YTick = 0:0.2:1;
ax.YTickLabel = 100*[0:0.2:1];


%%
figname = 'Fig3_ModelStats_LiveVsInfected'
set(gcf, 'PaperPositionMode','auto','color','w','InvertHardcopy','off')
print(gcf,'-depsc','-r600',[savpath '/epss/' figname]);
print(gcf,'-dpng','-r600',[savpath '/pngs/' figname]);


%% Plot Model statistics
close all
figure('color','w','Position',[100,100, 450, 450])
ax = axes('Position', [0.15, 0.14, 0.75, 0.75])
a = cellfun(@(x) x(:,1), Stats,'uniformOutput',false)
plot(asinh(TNF/LinearRange), median(cat(2,a{:})),'Color', liveColor);
set(gca,'YLim',[0,100])
ylabel('Percent of cells')
hold on
b = cellfun(@(x) x(:,2), Stats,'uniformOutput',false)
plot(asinh(TNF/LinearRange), median(cat(2,b{:})),'Color',virusColor);
set(gca,'YLim',[0,100])
hold on
c = cellfun(@(x) x(:,4), Stats,'uniformOutput',false)
plot(asinh(TNF/LinearRange), median(cat(2,c{:})),'Color',deadinfColor);
set(gca,'YLim',[0,100])


d = cellfun(@(x) x(:,3), Stats,'uniformOutput',false)
plot(asinh(TNF/LinearRange), median(cat(2,d{:})),'Color',deathColor);
set(gca,'YLim',[0,100])

hold on

h = boxplot(cat(2,a{:}), 'Positions',asinh(TNF/LinearRange), 'PlotStyle','compact', 'MedianStyle','line','Symbol','.','Widths',0.2, 'Labels', {'' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' },'Color', liveColor);
shg
xlabel('[TNF] ng/ml')
set(gca,'YLim',[-1,100]/100)
h = boxplot(cat(2,b{:}), 'Positions',asinh(TNF/LinearRange), 'PlotStyle','compact', 'MedianStyle','line','Symbol','.','Widths',0.2, 'Labels', {'' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' },'Color', virusColor);
shg
set(gca,'YLim',[-1,100]/100)


h = boxplot(cat(2,c{:}), 'Positions',asinh(TNF/LinearRange), 'PlotStyle','compact', 'MedianStyle','line','Symbol','.','Widths',0.2, 'Labels', {'' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' },'Color', deadinfColor);
shg
xlabel('[TNF] ng/ml')
set(gca,'YLim',[-1,100]/100)
h = boxplot(cat(2,d{:}), 'Positions',asinh(TNF/LinearRange), 'PlotStyle','compact', 'MedianStyle','line','Symbol','.','Widths',0.2, 'Labels',  {'' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' },'Color', deathColor);
shg
set(gca,'YLim',[-1,100]/100)

hleg = legend('   TN - Healthy','   FN - Infected','   TP - Dead infected','   FP - Dead uninfected');
hleg.Location = 'northwest';
hleg.Box = 'off';

ax.YTick = 0:0.2:1;
ax.YTickLabel = 100*[0:0.2:1];

a1 = sort(- (1:9)'*10.^(1:2));
a2 =(1:9)'*10.^(1:5);
a = [sort(a1(:)); 0; a2(:)]/1000;
emptyStr = {''; ''; '';''; ''; '';''; ''};
ax.XTickLabel = [emptyStr; ''; emptyStr; {''}; '0'; {''}; emptyStr; '  .1'; emptyStr; ' 1'; emptyStr; ' 10'; emptyStr; '    100'; emptyStr];
ax.XTick = asinh(a/LinearRange);

xl = xlabel('[TNF] ng/ml')
xl.Position(2) = xl.Position(2)



%%
figname = 'Fig3_ModelStats_AllTogether'
set(gcf, 'PaperPositionMode','auto','color','w','InvertHardcopy','off')
print(gcf,'-depsc','-r600',[savpath '/epss/' figname]);
print(gcf,'-dpng','-r600',[savpath '/pngs/' figname]);


%% Plot Model statistics
close all
[463 67 449 535]
figure('color','w','Position',[463 67 450 550])

ax = axes('Position', [0.14, 0.6, 0.33, 0.3])
a = cellfun(@(x) x(:,1), Stats,'uniformOutput',false)
plot(asinh(TNF/LinearRange), median(cat(2,a{:})),'Color', liveColor);
set(gca,'YLim',[0,1])
ylabel('TN - Healthy')
hold on
h = boxplot(cat(2,a{:}), 'Positions',asinh(TNF/LinearRange), 'PlotStyle','compact', 'MedianStyle','line','Symbol','.','Widths',0.2, 'Labels', num2str(TNF',3),'Color', liveColor);
shg
set(gca,'YLim',[-1,100]/100)
a1 = sort(- (1:9)'*10.^(1:2));
a2 =(1:9)'*10.^(1:5);
a = [sort(a1(:)); 0; a2(:)]/1000;
emptyStr = {''; ''; '';''; ''; '';''; ''};
ax.XTickLabel = [emptyStr; ''; emptyStr; {''}; '0'; {''}; emptyStr; '  .1'; emptyStr; ' 1'; emptyStr; ' 10'; emptyStr; '    100'; emptyStr];
ax.XTick = asinh(a/LinearRange);

xl = xlabel('[TNF] ng/ml');
xl.Position(2) = xl.Position(2)-5;
ax.YTick = [0 0.5 1]

ax = axes('Position', [0.64, 0.6, 0.33, 0.3])
b = cellfun(@(x) x(:,2), Stats,'uniformOutput',false)
plot(asinh(TNF/LinearRange), median(cat(2,b{:})),'Color', virusColor);
set(gca,'YLim',[0,1])
ylabel('FN - Infected')
hold on
h = boxplot(cat(2,b{:}), 'Positions',asinh(TNF/LinearRange), 'PlotStyle','compact', 'MedianStyle','line','Symbol','.','Widths',0.2, 'Labels', num2str(TNF',3),'Color', virusColor);
shg
set(gca,'YLim',[-1,100]/100)
a1 = sort(- (1:9)'*10.^(1:2));
a2 =(1:9)'*10.^(1:5);
a = [sort(a1(:)); 0; a2(:)]/1000;
emptyStr = {''; ''; '';''; ''; '';''; ''};
ax.XTickLabel = [emptyStr; ''; emptyStr; {''}; '0'; {''}; emptyStr; '  .1'; emptyStr; ' 1'; emptyStr; ' 10'; emptyStr; '    100'; emptyStr];
ax.XTick = asinh(a/LinearRange);

xl = xlabel('[TNF] ng/ml');
xl.Position(2) = xl.Position(2)-5;
ax.YTick = [0 0.5 1]

ax = axes('Position', [0.14, 0.15, 0.33, 0.3])
c = cellfun(@(x) x(:,4), Stats,'uniformOutput',false)
plot(asinh(TNF/LinearRange), median(cat(2,c{:})),'Color', deadinfColor);
set(gca,'YLim',[0,1])
hold on;
ylabel('TP - Dead infected')
h = boxplot(cat(2,c{:}), 'Positions',asinh(TNF/LinearRange), 'PlotStyle','compact', 'MedianStyle','line','Symbol','.','Widths',0.2, 'Labels', num2str(TNF',3),'Color', deadinfColor);
shg
set(gca,'YLim',[-1,100]/100)
a1 = sort(- (1:9)'*10.^(1:2));
a2 =(1:9)'*10.^(1:5);
a = [sort(a1(:)); 0; a2(:)]/1000;
emptyStr = {''; ''; '';''; ''; '';''; ''};
ax.XTickLabel = [emptyStr; ''; emptyStr; {''}; '0'; {''}; emptyStr; '  .1'; emptyStr; ' 1'; emptyStr; ' 10'; emptyStr; '    100'; emptyStr];
ax.XTick = asinh(a/LinearRange);

xl = xlabel('[TNF] ng/ml');
xl.Position(2) = xl.Position(2)-5;
ax.YTick = [0 0.5 1]

ax = axes('Position', [0.64, 0.15, 0.33, 0.3])
d = cellfun(@(x) x(:,3), Stats,'uniformOutput',false)
plot(asinh(TNF/LinearRange), median(cat(2,d{:})),'Color', deathColor);
set(gca,'YLim',[0,1])
ylabel('FP - Dead uninfected')
hold on;
h = boxplot(cat(2,d{:}), 'Positions',asinh(TNF/LinearRange), 'PlotStyle','compact', 'MedianStyle','line','Symbol','.','Widths',0.2, 'Labels', num2str(TNF',3),'Color', deathColor);
shg
set(gca,'YLim',[-1,100]/100)
a1 = sort(- (1:9)'*10.^(1:2));
a2 =(1:9)'*10.^(1:5);
a = [sort(a1(:)); 0; a2(:)]/1000;
emptyStr = {''; ''; '';''; ''; '';''; ''};
ax.XTickLabel = [emptyStr; ''; emptyStr; {''}; '0'; {''}; emptyStr; '  .1'; emptyStr; ' 1'; emptyStr; ' 10'; emptyStr; '    100'; emptyStr];
ax.XTick = asinh(a/LinearRange);

xl = xlabel('[TNF] ng/ml');
xl.Position(2) = xl.Position(2)-5;
ax.YTick = [0 0.5 1]



%%
figname = 'Fig3_ModelStats_AllSeparate'
set(gcf, 'PaperPositionMode','auto','color','w','InvertHardcopy','off')
print(gcf,'-depsc','-r600',[savpath '/epss/' figname]);
print(gcf,'-dpng','-r600',[savpath '/pngs/' figname]);


% 
% 
% %% Plot Model statistics compare to apoptosis inhibitor
% close all
% LinearRange = 0.4
% 
% figure('color','w','Position',[100,100, 450, 450])
% ax = axes('Position', [0.15, 0.14, 0.75, 0.75])
% 
% 
% 
% b = cellfun(@(x) x(:,2)+x(:,4), Stats,'uniformOutput',false)
% h0 = plot(asinh(TNF/LinearRange), median(cat(2,b{:})),'Color','k');
% set(gca,'YLim',[-1,100]/100)%,'YColor',virusColor)
% ax.XLim=asinh([0 100]/LinearRange)
% ylabel('Total Infected (TP+FN)')
% 
% xl = xlabel('[TNF] ng/ml')
% 
% title('Model')
% 
% hold on
% %h1 = boxplot(StatsNoDeath{1}(:,2), 'PlotStyle','compact', 'MedianStyle','line','Symbol','.','Widths',10,'Color',[0 0.5 0])
% 
% h = boxplot(cat(2,b{:}), 'Positions',asinh(TNF/LinearRange), 'PlotStyle','compact', 'MedianStyle','line','Symbol','.','Widths',0.3, 'Labels',  {'' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' },'Color', 'k');
% shg
% 
% 
% p = patch([0 asinh(150/LinearRange) asinh(150/LinearRange) 0], [prctile(StatsNoDeath{1}(:,2),1),prctile(StatsNoDeath{1}(:,2),1),prctile(StatsNoDeath{1}(:,2),99),prctile(StatsNoDeath{1}(:,2),99)],'r')
%  p.FaceAlpha = 0.1
%  p.EdgeColor = 'none'
%  
% p = patch([0 asinh(150/LinearRange) asinh(150/LinearRange) 0], [prctile(StatsNoDeath{1}(:,2),25),prctile(StatsNoDeath{1}(:,2),25),prctile(StatsNoDeath{1}(:,2),75),prctile(StatsNoDeath{1}(:,2),75)],'r')
%  p.FaceAlpha = 0.25
%  p.EdgeColor = 'none'
%  
%  h1 = plot([0 asinh(150/LinearRange)],[median(StatsNoDeath{1}(:,2)) median(StatsNoDeath{1}(:,2))],'Color','r')
% 
% a1 = sort(- (1:9)'*10.^(1:2));
% a2 =(1:9)'*10.^(1:5);
% a = [sort(a1(:)); 0; a2(:)]/1000;
% emptyStr = {''; ''; '';''; ''; '';''; ''};
% ax.XTickLabel = [emptyStr; ''; emptyStr; {''}; '0'; {''}; emptyStr; '  .1'; emptyStr; ' 1'; emptyStr; ' 10'; emptyStr; '    100'; emptyStr];
% ax.XTick = asinh(a/LinearRange);
% 
% ax.YTick = [0,0.1, 1];
% ax.YTickLabel = 100*[0,0.1, 1];
% 
% 
% ax.YScale = 'log'
% ax.YLim=[3 100]/100
% 
% 
% hleg = legend([h1, h0], 'TNF+Apoptosis inhibitor', 'TNF')
% hleg.Box='off'
% hleg.Position = [0.5541    0.6496    0.3311    0.1554]
% 
% %%
% figname = 'Fig3_ModelStats_AllinfVsNoApoptosis'
% set(gcf, 'PaperPositionMode','auto','color','w','InvertHardcopy','off','renderer','painters')
% print(gcf,'-depsc','-r600',[savpath '/epss/' figname]);
% print(gcf,'-dpng','-r600',[savpath '/pngs/' figname]);
% 
% 
% 
% 
% 
% 
% %% Get model stats
% TNF = R.getData('modelTNF')
% Stats = R.getData('modelSTATsInfEqualBasal')
% 
% 
% %% Plot Model statistics All equal basal
% close all
% LinearRange = 0.4
% 
% figure('color','w','Position',[100,100, 450, 450])
% ax = axes('Position', [0.15, 0.14, 0.75, 0.75])
% 
% a = cellfun(@(x) x(:,1), Stats,'uniformOutput',false)
% yyaxis(ax,'left')
% plot(asinh(TNF/LinearRange),median(cat(2,a{:})),'Color',liveColor);
% set(gca,'YLim',[-1,100]/100,'YColor',liveColor)
% ylabel('Percent Healthy (TN)')
% 
% b = cellfun(@(x) x(:,2), Stats,'uniformOutput',false)
% yyaxis(ax,'right')
% plot(asinh(TNF/LinearRange), median(cat(2,b{:})),'Color',virusColor);
% set(gca,'YLim',[-1,100]/100,'YColor',virusColor)
% ax.XLim=asinh([0 100]/LinearRange)
% ax.YLim=[-1 100]/100
% ylabel('Percent Infected (FN)')
% 
% xl = xlabel('[TNF] ng/ml')
% 
% title('Model')
% 
% hold on
% %
% yyaxis(ax,'left')
% h = boxplot(cat(2,a{:}), 'Positions',asinh(TNF/LinearRange) , 'PlotStyle','compact', 'MedianStyle','line','Symbol','.','Widths',0.3, 'Labels',{'' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' },'Color', liveColor);
% shg
% set(gca,'YLim',[-1,100]/100)
% ax.YTick = 0:0.2:1;
% ax.YTickLabel = 100*[0:0.2:1];
% 
% yyaxis(ax,'right')
% h = boxplot(cat(2,b{:}), 'Positions',asinh(TNF/LinearRange), 'PlotStyle','compact', 'MedianStyle','line','Symbol','.','Widths',0.3, 'Labels',  {'' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' },'Color', virusColor);
% shg
% set(gca,'YLim',[-1,100]/100)
% 
% 
% a1 = sort(- (1:9)'*10.^(1:2));
% a2 =(1:9)'*10.^(1:5);
% a = [sort(a1(:)); 0; a2(:)]/1000;
% emptyStr = {''; ''; '';''; ''; '';''; ''};
% ax.XTickLabel = [emptyStr; ''; emptyStr; {''}; '0'; {''}; emptyStr; '  .1'; emptyStr; ' 1'; emptyStr; ' 10'; emptyStr; '    100'; emptyStr];
% ax.XTick = asinh(a/LinearRange);
% 
% ax.YTick = 0:0.2:1;
% ax.YTickLabel = 100*[0:0.2:1];
% 
% 
% %%
% figname = 'Fig3_ModelStats_LiveVsInfected_AllEqualBasal'
% set(gcf, 'PaperPositionMode','auto','color','w','InvertHardcopy','off')
% print(gcf,'-depsc','-r600',[savpath '/epss/' figname]);
% print(gcf,'-dpng','-r600',[savpath '/pngs/' figname]);
% 
% 
% %% Plot Model statistics All equal basal
% close all
% figure('color','w','Position',[100,100, 450, 450])
% ax = axes('Position', [0.15, 0.14, 0.75, 0.75])
% a = cellfun(@(x) x(:,1), Stats,'uniformOutput',false)
% plot(asinh(TNF/LinearRange), median(cat(2,a{:})),'Color', liveColor);
% set(gca,'YLim',[0,100])
% ylabel('Percent of cells')
% hold on
% b = cellfun(@(x) x(:,2), Stats,'uniformOutput',false)
% plot(asinh(TNF/LinearRange), median(cat(2,b{:})),'Color',virusColor);
% set(gca,'YLim',[0,100])
% hold on
% c = cellfun(@(x) x(:,4), Stats,'uniformOutput',false)
% plot(asinh(TNF/LinearRange), median(cat(2,c{:})),'Color',deadinfColor);
% set(gca,'YLim',[0,100])
% 
% 
% d = cellfun(@(x) x(:,3), Stats,'uniformOutput',false)
% plot(asinh(TNF/LinearRange), median(cat(2,d{:})),'Color',deathColor);
% set(gca,'YLim',[0,100])
% 
% hold on
% 
% h = boxplot(cat(2,a{:}), 'Positions',asinh(TNF/LinearRange), 'PlotStyle','compact', 'MedianStyle','line','Symbol','.','Widths',0.2, 'Labels', {'' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' },'Color', liveColor);
% shg
% xlabel('[TNF] ng/ml')
% set(gca,'YLim',[-1,100]/100)
% h = boxplot(cat(2,b{:}), 'Positions',asinh(TNF/LinearRange), 'PlotStyle','compact', 'MedianStyle','line','Symbol','.','Widths',0.2, 'Labels', {'' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' },'Color', virusColor);
% shg
% set(gca,'YLim',[-1,100]/100)
% 
% 
% h = boxplot(cat(2,c{:}), 'Positions',asinh(TNF/LinearRange), 'PlotStyle','compact', 'MedianStyle','line','Symbol','.','Widths',0.2, 'Labels', {'' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' },'Color', deadinfColor);
% shg
% xlabel('[TNF] ng/ml')
% set(gca,'YLim',[-1,100]/100)
% h = boxplot(cat(2,d{:}), 'Positions',asinh(TNF/LinearRange), 'PlotStyle','compact', 'MedianStyle','line','Symbol','.','Widths',0.2, 'Labels',  {'' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' },'Color', deathColor);
% shg
% set(gca,'YLim',[-1,100]/100)
% 
% hleg = legend('   TN - Healthy','   FN - Infected','   TP - Dead infected','   FP - Dead uninfected');
% hleg.Location = 'northwest';
% hleg.Box = 'off';
% 
% ax.YTick = 0:0.2:1;
% ax.YTickLabel = 100*[0:0.2:1];
% 
% a1 = sort(- (1:9)'*10.^(1:2));
% a2 =(1:9)'*10.^(1:5);
% a = [sort(a1(:)); 0; a2(:)]/1000;
% emptyStr = {''; ''; '';''; ''; '';''; ''};
% ax.XTickLabel = [emptyStr; ''; emptyStr; {''}; '0'; {''}; emptyStr; '  .1'; emptyStr; ' 1'; emptyStr; ' 10'; emptyStr; '    100'; emptyStr];
% ax.XTick = asinh(a/LinearRange);
% 
% xl = xlabel('[TNF] ng/ml')
% xl.Position(2) = xl.Position(2)
% 
% 
% 
% %%
% figname = 'Fig3_ModelStats_AllTogether_AllEqualBasal'
% set(gcf, 'PaperPositionMode','auto','color','w','InvertHardcopy','off')
% print(gcf,'-depsc','-r600',[savpath '/epss/' figname]);
% print(gcf,'-dpng','-r600',[savpath '/pngs/' figname]);
% 
% 
% %% Plot Model statistics All equal basal
% close all
% [463 67 449 535]
% figure('color','w','Position',[463 67 450 550])
% 
% ax = axes('Position', [0.14, 0.6, 0.33, 0.3])
% a = cellfun(@(x) x(:,1), Stats,'uniformOutput',false)
% plot(asinh(TNF/LinearRange), median(cat(2,a{:})),'Color', liveColor);
% set(gca,'YLim',[0,1])
% ylabel('TN - Healthy')
% hold on
% h = boxplot(cat(2,a{:}), 'Positions',asinh(TNF/LinearRange), 'PlotStyle','compact', 'MedianStyle','line','Symbol','.','Widths',0.2, 'Labels', num2str(TNF',3),'Color', liveColor);
% shg
% set(gca,'YLim',[-1,100]/100)
% a1 = sort(- (1:9)'*10.^(1:2));
% a2 =(1:9)'*10.^(1:5);
% a = [sort(a1(:)); 0; a2(:)]/1000;
% emptyStr = {''; ''; '';''; ''; '';''; ''};
% ax.XTickLabel = [emptyStr; ''; emptyStr; {''}; '0'; {''}; emptyStr; '  .1'; emptyStr; ' 1'; emptyStr; ' 10'; emptyStr; '    100'; emptyStr];
% ax.XTick = asinh(a/LinearRange);
% 
% xl = xlabel('[TNF] ng/ml');
% xl.Position(2) = xl.Position(2)-5;
% ax.YTick = [0 1]
% 
% ax = axes('Position', [0.64, 0.6, 0.33, 0.3])
% b = cellfun(@(x) x(:,2), Stats,'uniformOutput',false)
% plot(asinh(TNF/LinearRange), median(cat(2,b{:})),'Color', virusColor);
% set(gca,'YLim',[0,1])
% ylabel('FN - Infected')
% hold on
% h = boxplot(cat(2,b{:}), 'Positions',asinh(TNF/LinearRange), 'PlotStyle','compact', 'MedianStyle','line','Symbol','.','Widths',0.2, 'Labels', num2str(TNF',3),'Color', virusColor);
% shg
% set(gca,'YLim',[-1,100]/100)
% a1 = sort(- (1:9)'*10.^(1:2));
% a2 =(1:9)'*10.^(1:5);
% a = [sort(a1(:)); 0; a2(:)]/1000;
% emptyStr = {''; ''; '';''; ''; '';''; ''};
% ax.XTickLabel = [emptyStr; ''; emptyStr; {''}; '0'; {''}; emptyStr; '  .1'; emptyStr; ' 1'; emptyStr; ' 10'; emptyStr; '    100'; emptyStr];
% ax.XTick = asinh(a/LinearRange);
% 
% xl = xlabel('[TNF] ng/ml');
% xl.Position(2) = xl.Position(2)-5;
% ax.YTick = [0 1]
% 
% ax = axes('Position', [0.14, 0.15, 0.33, 0.3])
% c = cellfun(@(x) x(:,4), Stats,'uniformOutput',false)
% plot(asinh(TNF/LinearRange), median(cat(2,c{:})),'Color', deadinfColor);
% set(gca,'YLim',[0,1])
% hold on;
% ylabel('TP - Dead infected')
% h = boxplot(cat(2,c{:}), 'Positions',asinh(TNF/LinearRange), 'PlotStyle','compact', 'MedianStyle','line','Symbol','.','Widths',0.2, 'Labels', num2str(TNF',3),'Color', deadinfColor);
% shg
% set(gca,'YLim',[-1,100]/100)
% a1 = sort(- (1:9)'*10.^(1:2));
% a2 =(1:9)'*10.^(1:5);
% a = [sort(a1(:)); 0; a2(:)]/1000;
% emptyStr = {''; ''; '';''; ''; '';''; ''};
% ax.XTickLabel = [emptyStr; ''; emptyStr; {''}; '0'; {''}; emptyStr; '  .1'; emptyStr; ' 1'; emptyStr; ' 10'; emptyStr; '    100'; emptyStr];
% ax.XTick = asinh(a/LinearRange);
% 
% xl = xlabel('[TNF] ng/ml');
% xl.Position(2) = xl.Position(2)-5;
% ax.YTick = [0 1]
% 
% ax = axes('Position', [0.64, 0.15, 0.33, 0.3])
% d = cellfun(@(x) x(:,3), Stats,'uniformOutput',false)
% plot(asinh(TNF/LinearRange), median(cat(2,d{:})),'Color', deathColor);
% set(gca,'YLim',[0,1])
% ylabel('FP - Dead uninfected')
% hold on;
% h = boxplot(cat(2,d{:}), 'Positions',asinh(TNF/LinearRange), 'PlotStyle','compact', 'MedianStyle','line','Symbol','.','Widths',0.2, 'Labels', num2str(TNF',3),'Color', deathColor);
% shg
% set(gca,'YLim',[-1,100]/100)
% a1 = sort(- (1:9)'*10.^(1:2));
% a2 =(1:9)'*10.^(1:5);
% a = [sort(a1(:)); 0; a2(:)]/1000;
% emptyStr = {''; ''; '';''; ''; '';''; ''};
% ax.XTickLabel = [emptyStr; ''; emptyStr; {''}; '0'; {''}; emptyStr; '  .1'; emptyStr; ' 1'; emptyStr; ' 10'; emptyStr; '    100'; emptyStr];
% ax.XTick = asinh(a/LinearRange);
% 
% xl = xlabel('[TNF] ng/ml');
% xl.Position(2) = xl.Position(2)-5;
% ax.YTick = [0 1]
% 
% 
% 
% %%
% figname = 'Fig3_ModelStats_AllSeparate_AllEqualBasal'
% set(gcf, 'PaperPositionMode','auto','color','w','InvertHardcopy','off')
% print(gcf,'-depsc','-r600',[savpath '/epss/' figname]);
% print(gcf,'-dpng','-r600',[savpath '/pngs/' figname]);
% 
% 
% 
% %% Get model stats
% TNF = R.getData('modelTNF')
% Stats = R.getData('modelSTATsBasalis0')
% 
% 
% %% Plot Model statistics All equal basal
% close all
% LinearRange = 0.4
% 
% figure('color','w','Position',[100,100, 450, 450])
% ax = axes('Position', [0.15, 0.14, 0.7, 0.7])
% 
% a = cellfun(@(x) x(:,1), Stats,'uniformOutput',false)
% yyaxis(ax,'left')
% plot(asinh(TNF/LinearRange),median(cat(2,a{:})),'Color',liveColor);
% set(gca,'YLim',[-1,100]/100,'YColor',liveColor)
% ylabel('Percent Healthy (TN)')
% 
% b = cellfun(@(x) x(:,2), Stats,'uniformOutput',false)
% yyaxis(ax,'right')
% plot(asinh(TNF/LinearRange), median(cat(2,b{:})),'Color',virusColor);
% set(gca,'YLim',[-1,100]/100,'YColor',virusColor)
% ax.XLim=asinh([0 100]/LinearRange)
% ax.YLim=[-1 100]/100
% ylabel('Percent Infected (FN)')
% 
% xl = xlabel('[TNF] ng/ml')
% 
% title('Model')
% 
% hold on
% %
% yyaxis(ax,'left')
% h = boxplot(cat(2,a{:}), 'Positions',asinh(TNF/LinearRange) , 'PlotStyle','compact', 'MedianStyle','line','Symbol','.','Widths',0.3, 'Labels',{'' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' },'Color', liveColor);
% shg
% set(gca,'YLim',[-1,100]/100)
% ax.YTick = 0:0.2:1;
% ax.YTickLabel = 100*[0:0.2:1];
% 
% yyaxis(ax,'right')
% h = boxplot(cat(2,b{:}), 'Positions',asinh(TNF/LinearRange), 'PlotStyle','compact', 'MedianStyle','line','Symbol','.','Widths',0.3, 'Labels',  {'' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' },'Color', virusColor);
% shg
% set(gca,'YLim',[-1,100]/100)
% 
% 
% a1 = sort(- (1:9)'*10.^(1:2));
% a2 =(1:9)'*10.^(1:5);
% a = [sort(a1(:)); 0; a2(:)]/1000;
% emptyStr = {''; ''; '';''; ''; '';''; ''};
% ax.XTickLabel = [emptyStr; ''; emptyStr; {''}; '0'; {''}; emptyStr; '  .1'; emptyStr; ' 1'; emptyStr; ' 10'; emptyStr; '    100'; emptyStr];
% ax.XTick = asinh(a/LinearRange);
% 
% ax.YTick = 0:0.2:1;
% ax.YTickLabel = 100*[0:0.2:1];
% 
% 
% %%
% figname = 'Fig3_ModelStats_LiveVsInfected_BasalEqualZero'
% set(gcf, 'PaperPositionMode','auto','color','w','InvertHardcopy','off')
% print(gcf,'-depsc','-r600',[savpath '/epss/' figname]);
% print(gcf,'-dpng','-r600',[savpath '/pngs/' figname]);
% 
% 
% %% Plot Model statistics All equal basal
% close all
% figure('color','w','Position',[100,100, 450, 450])
% ax = axes('Position', [0.15, 0.14, 0.75, 0.75])
% a = cellfun(@(x) x(:,1), Stats,'uniformOutput',false)
% plot(asinh(TNF/LinearRange), median(cat(2,a{:})),'Color', liveColor);
% set(gca,'YLim',[0,100])
% ylabel('Percent of cells')
% hold on
% b = cellfun(@(x) x(:,2), Stats,'uniformOutput',false)
% plot(asinh(TNF/LinearRange), median(cat(2,b{:})),'Color',virusColor);
% set(gca,'YLim',[0,100])
% hold on
% c = cellfun(@(x) x(:,4), Stats,'uniformOutput',false)
% plot(asinh(TNF/LinearRange), median(cat(2,c{:})),'Color',deadinfColor);
% set(gca,'YLim',[0,100])
% 
% 
% d = cellfun(@(x) x(:,3), Stats,'uniformOutput',false)
% plot(asinh(TNF/LinearRange), median(cat(2,d{:})),'Color',deathColor);
% set(gca,'YLim',[0,100])
% 
% hold on
% 
% h = boxplot(cat(2,a{:}), 'Positions',asinh(TNF/LinearRange), 'PlotStyle','compact', 'MedianStyle','line','Symbol','.','Widths',0.2, 'Labels', {'' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' },'Color', liveColor);
% shg
% xlabel('[TNF] ng/ml')
% set(gca,'YLim',[-1,100]/100)
% h = boxplot(cat(2,b{:}), 'Positions',asinh(TNF/LinearRange), 'PlotStyle','compact', 'MedianStyle','line','Symbol','.','Widths',0.2, 'Labels', {'' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' },'Color', virusColor);
% shg
% set(gca,'YLim',[-1,100]/100)
% 
% 
% h = boxplot(cat(2,c{:}), 'Positions',asinh(TNF/LinearRange), 'PlotStyle','compact', 'MedianStyle','line','Symbol','.','Widths',0.2, 'Labels', {'' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' },'Color', deadinfColor);
% shg
% xlabel('[TNF] ng/ml')
% set(gca,'YLim',[-1,100]/100)
% h = boxplot(cat(2,d{:}), 'Positions',asinh(TNF/LinearRange), 'PlotStyle','compact', 'MedianStyle','line','Symbol','.','Widths',0.2, 'Labels',  {'' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' },'Color', deathColor);
% shg
% set(gca,'YLim',[-1,100]/100)
% 
% hleg = legend('   TN - Healthy','   FN - Infected','   TP - Dead infected','   FP - Dead uninfected');
% hleg.Location = 'northwest';
% hleg.Box = 'off';
% 
% ax.YTick = 0:0.2:1;
% ax.YTickLabel = 100*[0:0.2:1];
% 
% a1 = sort(- (1:9)'*10.^(1:2));
% a2 =(1:9)'*10.^(1:5);
% a = [sort(a1(:)); 0; a2(:)]/1000;
% emptyStr = {''; ''; '';''; ''; '';''; ''};
% ax.XTickLabel = [emptyStr; ''; emptyStr; {''}; '0'; {''}; emptyStr; '  .1'; emptyStr; ' 1'; emptyStr; ' 10'; emptyStr; '    100'; emptyStr];
% ax.XTick = asinh(a/LinearRange);
% 
% xl = xlabel('[TNF] ng/ml')
% xl.Position(2) = xl.Position(2)
% 
% 
% 
% %%
% figname = 'Fig3_ModelStats_AllTogether_BasalEqualZero'
% set(gcf, 'PaperPositionMode','auto','color','w','InvertHardcopy','off')
% print(gcf,'-depsc','-r600',[savpath '/epss/' figname]);
% print(gcf,'-dpng','-r600',[savpath '/pngs/' figname]);


% %% Plot Model statistics All equal basal
% close all
% [463 67 449 535]
% figure('color','w','Position',[463 67 450 550])
% 
% ax = axes('Position', [0.14, 0.6, 0.33, 0.3])
% a = cellfun(@(x) x(:,1), Stats,'uniformOutput',false)
% plot(asinh(TNF/LinearRange), median(cat(2,a{:})),'Color', liveColor);
% set(gca,'YLim',[0,1])
% ylabel('TN - Healthy')
% hold on
% h = boxplot(cat(2,a{:}), 'Positions',asinh(TNF/LinearRange), 'PlotStyle','compact', 'MedianStyle','line','Symbol','.','Widths',0.2, 'Labels', num2str(TNF',3),'Color', liveColor);
% shg
% set(gca,'YLim',[-1,100]/100)
% a1 = sort(- (1:9)'*10.^(1:2));
% a2 =(1:9)'*10.^(1:5);
% a = [sort(a1(:)); 0; a2(:)]/1000;
% emptyStr = {''; ''; '';''; ''; '';''; ''};
% ax.XTickLabel = [emptyStr; ''; emptyStr; {''}; '0'; {''}; emptyStr; '  .1'; emptyStr; ' 1'; emptyStr; ' 10'; emptyStr; '    100'; emptyStr];
% ax.XTick = asinh(a/LinearRange);
% 
% xl = xlabel('[TNF] ng/ml');
% xl.Position(2) = xl.Position(2)-5;
% ax.YTick = [0 1]
% 
% ax = axes('Position', [0.64, 0.6, 0.33, 0.3])
% b = cellfun(@(x) x(:,2), Stats,'uniformOutput',false)
% plot(asinh(TNF/LinearRange), median(cat(2,b{:})),'Color', virusColor);
% set(gca,'YLim',[0,1])
% ylabel('FN - Infected')
% hold on
% h = boxplot(cat(2,b{:}), 'Positions',asinh(TNF/LinearRange), 'PlotStyle','compact', 'MedianStyle','line','Symbol','.','Widths',0.2, 'Labels', num2str(TNF',3),'Color', virusColor);
% shg
% set(gca,'YLim',[-1,100]/100)
% a1 = sort(- (1:9)'*10.^(1:2));
% a2 =(1:9)'*10.^(1:5);
% a = [sort(a1(:)); 0; a2(:)]/1000;
% emptyStr = {''; ''; '';''; ''; '';''; ''};
% ax.XTickLabel = [emptyStr; ''; emptyStr; {''}; '0'; {''}; emptyStr; '  .1'; emptyStr; ' 1'; emptyStr; ' 10'; emptyStr; '    100'; emptyStr];
% ax.XTick = asinh(a/LinearRange);
% 
% xl = xlabel('[TNF] ng/ml');
% xl.Position(2) = xl.Position(2)-5;
% ax.YTick = [0 1]
% 
% ax = axes('Position', [0.14, 0.15, 0.33, 0.3])
% c = cellfun(@(x) x(:,4), Stats,'uniformOutput',false)
% plot(asinh(TNF/LinearRange), median(cat(2,c{:})),'Color', deadinfColor);
% set(gca,'YLim',[0,1])
% hold on;
% ylabel('TP - Dead infected')
% h = boxplot(cat(2,c{:}), 'Positions',asinh(TNF/LinearRange), 'PlotStyle','compact', 'MedianStyle','line','Symbol','.','Widths',0.2, 'Labels', num2str(TNF',3),'Color', deadinfColor);
% shg
% set(gca,'YLim',[-1,100]/100)
% a1 = sort(- (1:9)'*10.^(1:2));
% a2 =(1:9)'*10.^(1:5);
% a = [sort(a1(:)); 0; a2(:)]/1000;
% emptyStr = {''; ''; '';''; ''; '';''; ''};
% ax.XTickLabel = [emptyStr; ''; emptyStr; {''}; '0'; {''}; emptyStr; '  .1'; emptyStr; ' 1'; emptyStr; ' 10'; emptyStr; '    100'; emptyStr];
% ax.XTick = asinh(a/LinearRange);
% 
% xl = xlabel('[TNF] ng/ml');
% xl.Position(2) = xl.Position(2)-5;
% ax.YTick = [0 1]
% 
% ax = axes('Position', [0.64, 0.15, 0.33, 0.3])
% d = cellfun(@(x) x(:,3), Stats,'uniformOutput',false)
% plot(asinh(TNF/LinearRange), median(cat(2,d{:})),'Color', deathColor);
% set(gca,'YLim',[0,1])
% ylabel('FP - Dead uninfected')
% hold on;
% h = boxplot(cat(2,d{:}), 'Positions',asinh(TNF/LinearRange), 'PlotStyle','compact', 'MedianStyle','line','Symbol','.','Widths',0.2, 'Labels', num2str(TNF',3),'Color', deathColor);
% shg
% set(gca,'YLim',[-1,100]/100)
% a1 = sort(- (1:9)'*10.^(1:2));
% a2 =(1:9)'*10.^(1:5);
% a = [sort(a1(:)); 0; a2(:)]/1000;
% emptyStr = {''; ''; '';''; ''; '';''; ''};
% ax.XTickLabel = [emptyStr; ''; emptyStr; {''}; '0'; {''}; emptyStr; '  .1'; emptyStr; ' 1'; emptyStr; ' 10'; emptyStr; '    100'; emptyStr];
% ax.XTick = asinh(a/LinearRange);
% 
% xl = xlabel('[TNF] ng/ml');
% xl.Position(2) = xl.Position(2)-5;
% ax.YTick = [0 1]
% 
% 
% 
% %%
% figname = 'Fig3_ModelStats_AllSeparate_BasalEqualZero'
% set(gcf, 'PaperPositionMode','auto','color','w','InvertHardcopy','off')
% print(gcf,'-depsc','-r600',[savpath '/epss/' figname]);
% print(gcf,'-dpng','-r600',[savpath '/pngs/' figname]);
% 



%% spread in cell diameters as a function of dT
figure('color','w','Position',[100,100, 450, 450])
ax2 = axes('Position', [0.15 0.15 0.7 0.7])

b = cellfun(@(x) x(:,2)+x(:,4), Stats,'uniformOutput',false)
mean(cat(2,b{:}));
A = triangleGrid([-10 -10 10 10], [0,0], 1);
nCells = size(A,1);

ViralSpreadInCD = sqrt(nCells*mean(cat(2,b{:})))/pi;
err_ViralSpreadInCD = std(cat(2,b{:})).*sqrt(nCells./mean(cat(2,b{:})))/pi;

infDeathRate = interp1(cTNF,filtfilt([1,1,1],3,infDeathRatesExp),TNF)

h1 = ploterr(1./infDeathRate,ViralSpreadInCD,[],err_ViralSpreadInCD,'-o','logx');
BarTicklength(h1(2),0)
for i=1:length(h1) 
    h1(i).Color = virusColor;
    h1(i).MarkerSize=7
    h1(i).MarkerFaceColor = virusColor;
    
end
h1(2).Color = 'none'
xlabel('dT (h)')
ylabel('Viral spread (cd)')
set(gca,'XLim',[4.5,70],'XTick', [5 10 50])

%%
figname = 'Fig3_dTvsSpread_inCellDiameter_Model_NoError'
set(gcf, 'PaperPositionMode','auto','color','w','InvertHardcopy','off')
print(gcf,'-depsc','-r600',[savpath '/epss/' figname]);
print(gcf,'-dpng','-r600',[savpath '/pngs/' figname]);












%% Load low-moi with replicates

BaseStr = regexprep([char(ispc.*'Z:\Images2019\') char(isunix.*'/bigstore/Images2019/')],char(0),'');
Usr = 'Jen';
Project = 'NFkBDynamics';
Dataset = 'TNF_Titr_HSV_MOI1_withReplic_Rep3_2019Sep03';
acquisition = 1;

% Get MD of raw data
acqname = ['acq_' num2str(acquisition)];
fpath = [BaseStr Usr filesep Project filesep Dataset filesep acqname];
MD=Metadata(fpath);
Wells = unique(MD.getSpecificMetadata('Position'));
frames = unique(cell2mat(MD.getSpecificMetadata('frame')));
dt = 1/3;%h

cTNF = ([100*(1/2).^[0:8], 0]);

% Load Results

R = MultiPositionSingleCellVirusResults(fpath)
Wells = R.PosNames;
frames = R.Frames;
t = [R.Frames-1]'.*dt;

%% Plot Sweet Spot
figure('color','w','Position',[100,100, 450, 450])
ax2 = axes('Position', [0.15 0.15 0.7 0.7])

%axes('Position', [0.12, 0.14, 0.37, 0.75])

t = R.Frames'/3;
TNFracatT = [];
FNFracatT = [];
FPFracatT = [];
TPFracatT = [];
range = [31:60]
tzeva = flipud(parula(numel(range)));
for i=1:numel(range);
    pos = R.PosNames{range(i)};
    %deathMat = (R.putInMat(R.getTracks(R.PosNames{i}),'DeadSum'));
    infMat = R.putInMat(R.getTracks(pos),'Infected');
    deathMat = R.putInMat(R.getTracks(pos),'Dead');
    timept = 58; %48h +10h for initial infection to start (model starts with an infecting cell)
    TNFrac = sum(single(infMat==0).*single(deathMat==0))./(sum(deathMat(:,40)==1)+sum(deathMat(:,40)==0));
    FNFrac = sum(single(infMat==1).*single(deathMat==0))./(sum(deathMat(:,40)==1)+sum(deathMat(:,40)==0));
    FPFrac = sum(single(infMat==0).*single(deathMat==1))./(sum(deathMat(:,40)==1)+sum(deathMat(:,40)==0));
    TPFrac = sum(single(infMat==1).*single(deathMat==1))./(sum(deathMat(:,40)==1)+sum(deathMat(:,40)==0));
    
    TNFracatT = [TNFracatT, TNFrac(t==timept)];
    FNFracatT = [FNFracatT, FNFrac(t==timept)];
    FPFracatT = [FPFracatT, FPFrac(t==timept)];
    TPFracatT = [TPFracatT, TPFrac(t==timept)];
end

TNFracatT = reshape(TNFracatT,[],3)';
FNFracatT = reshape(FNFracatT,[],3)';
FPFracatT = reshape(FPFracatT,[],3)';
TPFracatT = reshape(TPFracatT,[],3)';

LinearRange = 0.4

yyaxis(ax2,'left')
h1 = ploterr(asinh(cTNF/LinearRange),mean(TNFracatT),[],std(TNFracatT),'-');
for i=1:length(h1) h1(i).Color = liveColor; end
hold on
set(gca,'YLim',[-1,100]/100,'YTick',[0 1],'YColor',liveColor)
%plot(asinh(cTNF/LinearRange),mean(FracHealthyat30),'-ob')
ylabel('Fraction Healthy (TN)')

xl = xlabel('[TNF\alpha] ng/ml')
title('Experiment')

ax2.YTick = 0:0.2:1;
ax2.YTickLabel = 100*[0:0.2:1];


%plot(asinh(cTNF/LinearRange),mean(FracInfat30),'-or')
yyaxis(ax2,'right')
h2 = ploterr(asinh(cTNF/LinearRange),mean(FNFracatT),[],std(FNFracatT),'-');
for i=1:length(h2) h2(i).Color = virusColor; end
set(gca,'YLim',[-1,100]/100,'YTick',[0 1],'YColor',virusColor)
ylabel('Fraction Infected (FN)')


%plot(asinh(cTNF/LinearRange),mean(FracInfat30),'-or')
% h3 = ploterr(asinh(cTNF/LinearRange),mean(FPFracatT),[],std(FPFracatT),'-');
% for i=1:length(h3) h3(i).Color = deathColor; end
% 
% %plot(asinh(cTNF/LinearRange),mean(FracInfat30),'-or')
% h4 = ploterr(asinh(cTNF/LinearRange),mean(TPFracatT),[],std(TPFracatT),'-');
% for i=1:length(h4) h4(i).Color = deadinfColor; end


a1 = sort(- (1:9)'*10.^(1:2));
a2 =(1:9)'*10.^(1:5);
a = [sort(a1(:)); 0; a2(:)]/1000;
emptyStr = {''; ''; '';''; ''; '';''; ''};
ax2.XTickLabel = [emptyStr; ''; emptyStr; {''}; '0'; {''}; emptyStr; '  .1'; emptyStr; ' 1'; emptyStr; ' 10'; emptyStr; '    100'; emptyStr];
ax2.XTick = asinh(a/LinearRange);



%ax2.XTick = asinh(fliplr([100*(1/10).^[0:9], 0])/LinearRange);
%ax2.XTickLabel= num2str(fliplr([100*(1/10).^[0:9], 0])','%3.1f');
ax2.XLim=asinh([-0.1 100]/LinearRange)
ax2.YLim=[-1 100]/100
ax2.YTick = 0:0.2:1;
ax2.YTickLabel = 100*[0:0.2:1];

ax2.XLim=[-0.3000 6.5146]



%%
figname = 'Fig3_LowMOI_SweetSPot_LivevsInfected'
set(gcf, 'PaperPositionMode','auto','color','w','InvertHardcopy','off')
print(gcf,'-depsc','-r600',[savpath '/epss/' figname]);
print(gcf,'-dpng','-r600',[savpath '/pngs/' figname]);

%%% Plot Sweet Spot
% %%
% figure('color','w','Position',[100,100, 450, 450])
% ax2 = axes('Position', [0.15 0.15 0.7 0.7])
% 
% LinearRange = 0.4
% 
% tp = 48*3-1
%     cla
% hold on
% TP = R.getData('TP');
% TP = cat(1,TP{:});
% 
% 
% FP = R.getData('FP');
% FP = cat(1,FP{:});
% 
% 
% TN = R.getData('TN');
% TN = cat(1,TN{:});
% 
% 
% FN = R.getData('FN');
% FN = cat(1,FN{:});
% 
% totalCells = TP+FP+FN+TN;
% tc = reshape(totalCells(:,tp),6,[])'
% 
% TP = reshape(TP(:,tp),6,[])';
% h1 = ploterr(asinh(cTNF/LinearRange),mean(TP(:,4:6)')./mean(mean(tc(:,4:6)')),[],std(TP(:,4:6)')./mean(mean(tc(:,4:6)')))
% for i=1:length(h1) h1(i).Color = deadinfColor; end
% BarTicklength(h1(2),0);
% FP = reshape(FP(:,tp),6,[])';
% h2 = ploterr(asinh(cTNF/LinearRange),mean(FP(:,4:6)')./mean(mean(tc(:,4:6)')),[],std(FP(:,4:6)')./mean(mean(tc(:,4:6)')))
% for i=1:length(h2) h2(i).Color = deathColor; end
% BarTicklength(h2(2),0)
% 
% TN = reshape(TN(:,tp),6,[])';
% h3 = ploterr(asinh(cTNF/LinearRange),mean(TN(:,4:6)')./mean(mean(tc(:,4:6)')),[],std(TN(:,4:6)')./mean(mean(tc(:,4:6)')))
% for i=1:length(h3) h3(i).Color = liveColor; end
% BarTicklength(h3(2),0)
% 
% FN = reshape(FN(:,tp),6,[])';
% h4 = ploterr(asinh(cTNF/LinearRange),mean(FN(:,4:6)')./mean(mean(tc(:,4:6)')),[],std(FN(:,4:6)')./mean(mean(tc(:,4:6)')))
% for i=1:length(h4) h4(i).Color = virusColor; end
% BarTicklength(h4(2),0)
% 
% set(gca,'YLim',[-1,100]/100,'YTick',[0 1])
% 
% a1 = sort(- (1:9)'*10.^(1:2));
% a2 =(1:9)'*10.^(1:5);
% a = [sort(a1(:)); 0; a2(:)]/1000;
% emptyStr = {''; ''; '';''; ''; '';''; ''};
% ax2.XTickLabel = [emptyStr; ''; emptyStr; {''}; '0'; {''}; emptyStr; '  .1'; emptyStr; ' 1'; emptyStr; ' 10'; emptyStr; '    100'; emptyStr];
% ax2.XTick = asinh(a/LinearRange);
% 
% ax2.XLim=asinh([-0.1 100]/LinearRange)
% ax2.YLim=[-1 100]/100
% ax2.YTick = 0:0.2:1;
% ax2.YTickLabel = 100*[0:0.2:1];
% ax2.Box = 'on';
% 
% xlabel('[TNF] ng/ml')
% ylabel('Percent of cells')
% 
% hleg = legend([h3(2) h4(2) h1(2) h2(2)],'   TN - Healthy','   FN - Infected','   TP - Dead infected','   FP - Dead uninfected');
% hleg.Box='off';
% hleg.Position(1) = hleg.Position(1);
% shg
%pause;


% %%
% figure('color','w','Position',[100,100, 450, 450])
% ax2 = axes('Position', [0.15 0.15 0.7 0.7])
% 
% LinearRange = 0.4
% 
% tp = 48*3+1
%     cla
% hold on
% TP = R.getData('TP');
% TP = cat(1,TP{:});
% 
% 
% FP = R.getData('FP');
% FP = cat(1,FP{:});
% 
% 
% TN = R.getData('TN');
% TN = cat(1,TN{:});
% 
% 
% FN = R.getData('FN');
% FN = cat(1,FN{:});
% 
% totalCells = TP+FP+FN+TN;
% tc = reshape(totalCells(:,tp),6,[])'
% 
% TP = reshape(TP(:,tp),6,[])';
% h1 = ploterr(asinh(cTNF/LinearRange),mean(TP(:,4:6)'./tc(:,4:6)'),[],std(TP(:,4:6)'./tc(:,4:6)'))
% for i=1:length(h1) h1(i).Color = deadinfColor; end
% BarTicklength(h1(2),0);
% FP = reshape(FP(:,tp),6,[])';
% h2 = ploterr(asinh(cTNF/LinearRange),mean(FP(:,4:6)'./tc(:,4:6)'),[],std(FP(:,4:6)'./tc(:,4:6)'))
% for i=1:length(h2) h2(i).Color = deathColor; end
% BarTicklength(h2(2),0)
% 
% TN = reshape(TN(:,tp),6,[])';
% h3 = ploterr(asinh(cTNF/LinearRange),mean(TN(:,4:6)'./tc(:,4:6)'),[],std(TN(:,4:6)'./tc(:,4:6)'))
% for i=1:length(h3) h3(i).Color = liveColor; end
% BarTicklength(h3(2),0)
% 
% FN = reshape(FN(:,tp),6,[])';
% h4 = ploterr(asinh(cTNF/LinearRange),mean(FN(:,4:6)'./tc(:,4:6)'),[],std(FN(:,4:6)'./tc(:,4:6)'))
% for i=1:length(h4) h4(i).Color = virusColor; end
% BarTicklength(h4(2),0)
% 
% set(gca,'YLim',[-1,100]/100,'YTick',[0 1])
% 
% a1 = sort(- (1:9)'*10.^(1:2));
% a2 =(1:9)'*10.^(1:5);
% a = [sort(a1(:)); 0; a2(:)]/1000;
% emptyStr = {''; ''; '';''; ''; '';''; ''};
% ax2.XTickLabel = [emptyStr; ''; emptyStr; {''}; '0'; {''}; emptyStr; '  .1'; emptyStr; ' 1'; emptyStr; ' 10'; emptyStr; '    100'; emptyStr];
% ax2.XTick = asinh(a/LinearRange);
% 
% ax2.XLim=asinh([-0.1 100]/LinearRange)
% ax2.YLim=[-1 100]/100
% ax2.YTick = 0:0.2:1;
% ax2.YTickLabel = 100*[0:0.2:1];
% ax2.Box = 'on';
% 
% xlabel('[TNF] ng/ml')
% ylabel('Percent of cells')
% 
% hleg = legend([h3(2) h4(2) h1(2) h2(2)],'   TN - Healthy','   FN - Infected','   TP - Dead infected','   FP - Dead uninfected');
% hleg.Box='off';
% hleg.Position(1) = hleg.Position(1);
% shg
% %pause;


%%
%% Plot Sweet Spot All
figure('color','w','Position',[100,100, 450, 450])
ax2 = axes('Position', [0.15 0.15 0.7 0.7])

%axes('Position', [0.12, 0.14, 0.37, 0.75])

t = R.Frames'/3;
TNFracatT = [];
FNFracatT = [];
FPFracatT = [];
TPFracatT = [];
range = [31:60]
tzeva = flipud(parula(numel(range)));
tp0=40;
for i=1:numel(range);
    pos = R.PosNames{range(i)};
    %deathMat = (R.putInMat(R.getTracks(R.PosNames{i}),'DeadSum'));
    infMat = R.putInMat(R.getTracks(pos),'Infected');
    deathMat = R.putInMat(R.getTracks(pos),'Dead');
    timept = 58;%48h +10h for initial infection to start, model starts with infectious cell
    TNFrac = sum(single(infMat==0).*single(deathMat==0))./(sum(deathMat(:,tp0)==1)+sum(deathMat(:,tp0)==0));
    FNFrac = sum(single(infMat==1).*single(deathMat==0))./(sum(deathMat(:,tp0)==1)+sum(deathMat(:,tp0)==0));
    FPFrac = sum(single(infMat==0).*single(deathMat==1))./(sum(deathMat(:,tp0)==1)+sum(deathMat(:,tp0)==0));
    TPFrac = sum(single(infMat==1).*single(deathMat==1))./(sum(deathMat(:,tp0)==1)+sum(deathMat(:,tp0)==0));
    
    TNFracatT1 = TNFrac(t==timept);
    FNFracatT1 = FNFrac(t==timept);
    TPFracatT1 = TPFrac(t==timept);
    FPFracatT1 = FPFrac(t==timept)+max(1-TNFracatT1-FNFracatT1-TPFracatT1-FPFrac(t==timept),0);%disappeared cells are presumed dead
    
    
    TNFracatT = [TNFracatT, TNFracatT1];
    FNFracatT = [FNFracatT, FNFracatT1];
    FPFracatT = [FPFracatT, FPFracatT1];
    TPFracatT = [TPFracatT, TPFracatT1];
end

TNFracatT = reshape(TNFracatT,[],3)';
FNFracatT = reshape(FNFracatT,[],3)';
FPFracatT = reshape(FPFracatT,[],3)';
TPFracatT = reshape(TPFracatT,[],3)';

LinearRange = 0.4

h1 = ploterr(asinh(cTNF/LinearRange),mean(TNFracatT),[],std(TNFracatT),'-');
for i=1:length(h1) h1(i).Color = liveColor; h1(i).MarkerFaceColor = liveColor; end
hold on
%plot(asinh(cTNF/LinearRange),mean(FracHealthyat30),'-ob')
BarTicklength(h1(2),0)
h1(1).Marker='o';
h1(1).MarkerSize=7;

xlabel('[TNF] ng/ml')
ylabel('Fraction of initial population')



%plot(asinh(cTNF/LinearRange),mean(FracInfat30),'-or')
h2 = ploterr(asinh(cTNF/LinearRange),mean(FNFracatT),[],std(FNFracatT),'-');
for i=1:length(h2) h2(i).Color = virusColor;h2(i).MarkerFaceColor = virusColor; end

BarTicklength(h2(2),0)
h2(1).Marker='o';
h2(1).MarkerSize=7;

%plot(asinh(cTNF/LinearRange),mean(FracInfat30),'-or')
 h3 = ploterr(asinh(cTNF/LinearRange),mean(FPFracatT),[],std(FPFracatT),'-');
 for i=1:length(h3) h3(i).Color = deathColor;h3(i).MarkerFaceColor = deathColor; end
 BarTicklength(h3(2),0)
h3(1).Marker='o';
h3(1).MarkerSize=7;
% 
% %plot(asinh(cTNF/LinearRange),mean(FracInfat30),'-or')
 h4 = ploterr(asinh(cTNF/LinearRange),mean(TPFracatT),[],std(TPFracatT),'-');
 for i=1:length(h4) h4(i).Color = deadinfColor; h4(i).MarkerFaceColor = deadinfColor; end

BarTicklength(h4(2),0)
h4(1).Marker='o';
h4(1).MarkerSize=7;
a1 = sort(- (1:9)'*10.^(1:2));
a2 =(1:9)'*10.^(1:5);
a = [sort(a1(:)); 0; a2(:)]/1000;
emptyStr = {''; ''; '';''; ''; '';''; ''};
ax2.XTickLabel = [emptyStr; ''; emptyStr; {''}; '0'; {''}; emptyStr; '  .1'; emptyStr; ' 1'; emptyStr; ' 10'; emptyStr; '    100'; emptyStr];
ax2.XTick = asinh(a/LinearRange);


set(gca,'YLim',[-1,100]/100,'YTick',[0 1])

%ax2.XTick = asinh(fliplr([100*(1/10).^[0:9], 0])/LinearRange);
%ax2.XTickLabel= num2str(fliplr([100*(1/10).^[0:9], 0])','%3.1f');

ax2.YTick = 0:0.2:1;
ax2.YTickLabel = [0:0.2:1];
ax2.XLim=asinh([-0.1 100]/LinearRange)
ax2.YLim=[-1 100]/100


ax2.XLim=[-0.3000 6.5146]


%%
figname = 'Fig3_LowMOI_SweetSPot_AllTogether'
set(gcf, 'PaperPositionMode','auto','color','w','InvertHardcopy','off','Renderer', 'painters')
print(gcf,'-depsc','-r600',[savpath '/epss/' figname]);
print(gcf,'-dpng','-r600',[savpath '/pngs/' figname]);




%% Plot Model v Exp, Rsq


alldata = [mean(TNFracatT);mean(FNFracatT); mean(TPFracatT);mean(FPFracatT)]';
cTNF;

a = cellfun(@(x) x(:,1), Stats,'uniformOutput',false);
b = cellfun(@(x) x(:,2), Stats,'uniformOutput',false);
c = cellfun(@(x) x(:,4), Stats,'uniformOutput',false);
d = cellfun(@(x) x(:,3), Stats,'uniformOutput',false);
TNF;
allmodel = [mean(cat(2,a{:}));mean(cat(2,b{:}));mean(cat(2,c{:}));mean(cat(2,d{:}))]';

%
y=alldata;
yfit = interp1(TNF, allmodel, cTNF);


SStot = sum(sum((bsxfun(@minus, y,median(y))).^2));                            % Total Sum-Of-Squares
SSres = sum((y(:)-yfit(:)).^2);                         % Residual Sum-Of-Squares
Rsq = 1-SSres/SStot                             % R^2


figure('color','w','Position',[100,100, 900, 450])
ax2 = axes('Position', [0.18 0.18 0.7 0.7])

h = plot(asinh(TNF/0.4),allmodel,'--', asinh(cTNF/0.4),alldata,'-o','MarkerSize',8)
h(1).Color = liveColor;
h(5).Color = liveColor;
h(5).MarkerFaceColor = liveColor;
h(2).Color = virusColor;
h(6).Color = virusColor;
h(6).MarkerFaceColor = virusColor;
h(3).Color = deadinfColor;
h(7).Color = deadinfColor;
h(7).MarkerFaceColor = deadinfColor;
h(4).Color = deathColor;
h(8).Color = deathColor;
h(8).MarkerFaceColor = deathColor;


a1 = sort(- (1:9)'*10.^(1:2));
a2 =(1:9)'*10.^(1:5);
a = [sort(a1(:)); 0; a2(:)]/1000;
emptyStr = {''; ''; '';''; ''; '';''; ''};
ax2.XTickLabel = [emptyStr; ''; emptyStr; {''}; '0'; {''}; emptyStr; '  .1'; emptyStr; ' 1'; emptyStr; ' 10'; emptyStr; '    100'; emptyStr];
ax2.XTick = asinh(a/LinearRange);

xlabel('[TNF] ng/ml')
ylabel('Fraction of initial population')


ax2.YTick = 0:0.2:1;
ax2.YTickLabel = [0:0.2:1];
ax2.XLim=asinh([0 100]/LinearRange)
ax2.YLim=[-1 100]/100


a = annotation('textbox', [.3 .5 .3 .3],'String',['R^2=' num2str(Rsq,'%.3f')])
a.EdgeColor='None'


    ColumnNames = {'Model','Data'};
    RowNames = {'   TN - Healthy','   FN - Infected','   TP - Dead infected','   FP - Dead uninfected'};
    gridlegend(RowNames,ColumnNames,'Location','eastOutside','box','off');

%%
figname = 'FigS4_LowMOI_SweetSPot_AllTogether_CompareToModelWithRsq'
set(gcf, 'PaperPositionMode','auto','color','w','InvertHardcopy','off','Renderer', 'painters')
print(gcf,'-depsc','-r600',[savpath '/epss/' figname]);
print(gcf,'-dpng','-r600',[savpath '/pngs/' figname]);

%% Get model stats
TNF = R.getData('modelTNF')
Stats = R.getData('modelSTATs')

allParams = R.getData('ExpFitParam');
allParams = cat(1,allParams{:})

allParams = allParams(:,2);
allParams = reshape(allParams,[],6);

infDeathRatesExp = mean(allParams(:,4:6),2);
uninfDeathRatesExp = mean(allParams(:,1:3),2);






% %% Experimental sweet spot: Load results from low MOI expt
% %%do post analysis
%
% BaseStr = regexprep([char(ispc.*'Z:\Images2019\') char(isunix.*'/bigstore/Images2019/')],char(0),'');
% Usr = 'Jen';
% Project = 'NFkBDynamics';
% Dataset = 'TNF_DRUGS_HSVMOI1_2019Jun20';
% acquisition = 1;
%
% % Get MD of raw data
% acqname = ['acq_' num2str(acquisition)];
% fpath = [BaseStr Usr filesep Project filesep Dataset filesep acqname];
% MD=Metadata(fpath);
% Wells = unique(MD.getSpecificMetadata('Position'));
% frames = unique(cell2mat(MD.getSpecificMetadata('frame')));
% dt = 1/3;%h
%
% cTNF = ([100*(1/2).^[0:8], 0]);
% % Load Results
%
% R = MultiPositionSingleCellVirusResults(fpath)
% Wells = R.PosNames;
% frames = R.Frames;
% t = [R.Frames-1]'.*dt;


% %% Plot Sweet spot
% figure('color','w','Position',[100,100, 450, 450])
% %axes('Position', [0.12, 0.14, 0.37, 0.75])
%
% t = R.Frames'/3;
% FracHealthyat30 = []
% FracInfat30 = []
% FracDeadInf30 = []
% FracDeadBys30 = []
%
% range = [31:40]
% tzeva = flipud(parula(numel(range)));
% for i=1:numel(range);
%     pos = R.PosNames{range(i)};
% %deathMat = (R.putInMat(R.getTracks(R.PosNames{i}),'DeadSum'));
% infMat = R.putInMat(R.getTracks(pos),'Infected');
% deathMat = R.putInMat(R.getTracks(pos),'Dead');
% timept = 34;
% survivorFrac = sum(single(infMat==0).*single(deathMat==0))./(sum(deathMat==1)+sum(deathMat==0));
% infFrac = sum(single(infMat==1).*single(deathMat==0))./(sum(deathMat==1)+sum(deathMat==0));
% DeadInfFrac = sum(single(infMat==1).*single(deathMat==1))./(sum(deathMat==1)+sum(deathMat==0));
% DeadByFrac = sum(single(infMat==0).*single(deathMat==1))./(sum(deathMat==1)+sum(deathMat==0));
%
% FracHealthyat30 = [FracHealthyat30, survivorFrac(t==timept)];
% FracInfat30 = [FracInfat30, infFrac(t==timept)];
% FracDeadInf30 = [FracDeadInf30, DeadInfFrac(t==timept)];
% FracDeadBys30 = [FracDeadBys30, DeadByFrac(t==timept)];
% %h=semilogy(t, (survivorFrac),'o');h
% %h.Color = tzeva(i,:);
% %hold all
%
% %BETA = R.getData('ExpFitParam',pos);
% % y = Exponent(BETA,t);
% % h1 = semilogy(t, y)
% % h1.Color = tzeva(i,:);
% % h1.LineWidth=1.5
%
% end
% % xlabel('time(h)')
% %set(gca,'XLim',[0,max(t)],'YLim',[0.05,1],'YTick',[0.1 0.25 0.5 0.75 1],'YTickLabel',100*[],'TickLength', [0.03 0.01])
% %title('+HSV-I')
%
% % cb = colorbar('Position', [0.9 0.1400 0.0324 0.7500]);
% % cb.Ticks = asinh(fliplr([100*(1/2).^[0:9], 0]))./asinh(100);
% % cb.TickLabels= num2str(fliplr([100*(1/2).^[0:9], 0])','%3.1f');
% % cb.Label.String = '[TNF\alpha] (ng/ml)'
% % cb.Label.FontSize=12
% % cb.Label.Position(1) = 0.2
% % cb.Label.Color='w'
% nTNF = numel(cTNF);
% ax = axes('Position', [0.15, 0.2, 0.7, 0.7])
% LinearRange = 0.4
% %h1 = scatter(nTNF:-1:1,FracHealthyat30,[],tzeva,'filled');
% %hold on
% yyaxis(ax,'left')
%
% h1 = plot(nTNF:-1:1,100*FracHealthyat30,'-ob')
% h1.MarkerFaceColor = 'b';
% set(gca,'YLim',[-1,100],'Ycolor','b')
% ylabel('%Healthy')
%
% yyaxis(ax,'right')
%
% %h1 = scatter(nTNF:-1:1,FracInfat30,[],tzeva,'filled');
% hold on
% h2 = plot(nTNF:-1:1,100*FracInfat30,'-or')
% h2.MarkerFaceColor = 'r'
% set(gca,'YLim',[-1,100],'Ycolor','r')
% ylabel('%Infected')
%
%
% ax.XTick = 1:nTNF;
% ax.XTickLabel= num2str(fliplr(cTNF)','%3.1f');
% ax.XTickLabelRotation = 90
% ax.XLim=[1 nTNF]
% ax.YLim=[-1 100]
%
% xl = xlabel('[TNF\alpha] ng/ml')
% %xl.Position(2) = xl.Position(2)-30
% title('Experiment')
%
% %%
% figure(gcf)
% set(gcf, 'PaperPositionMode','auto','color','w','InvertHardcopy','off')
% print(gcf,'-depsc','-r600',[savpath 'SweetSpotExperimental'])
%

% %% Plot Sweet spot
% figure('color','w','Position',[100,100, 450, 450])
% %axes('Position', [0.12, 0.14, 0.37, 0.75])
%
% LinearRange = 0.4;
% a1 = sort(- (1:9)'*10.^(1:2));
% a2 =(1:9)'*10.^(1:5);
% a = [sort(a1(:)); 0; a2(:)]/1000;
% emptyStr = {''; ''; '';''; ''; '';''; ''};
% XTickLabel = [emptyStr; ''; emptyStr; {''}; '0'; {''}; emptyStr; ' 0.1'; emptyStr; ' 1'; emptyStr; ' 10'; emptyStr; '    100'; emptyStr];
% XTicks = asinh(a/LinearRange);
%
% t = R.Frames'/3;
% FracHealthyat30 = []
% FracInfat30 = []
% FracDeadInf30 = []
% FracDeadBys30 = []
%
% range = [31:40]
% tzeva = flipud(parula(numel(range)));
% for i=1:numel(range);
%     pos = R.PosNames{range(i)};
% %deathMat = (R.putInMat(R.getTracks(R.PosNames{i}),'DeadSum'));
% infMat = R.putInMat(R.getTracks(pos),'Infected');
% deathMat = R.putInMat(R.getTracks(pos),'Dead');
% timept = 34;
% survivorFrac = sum(single(infMat==0).*single(deathMat==0))./(sum(deathMat==1)+sum(deathMat==0));
% infFrac = sum(single(infMat==1).*single(deathMat==0))./(sum(deathMat==1)+sum(deathMat==0));
% DeadInfFrac = sum(single(infMat==1).*single(deathMat==1))./(sum(deathMat==1)+sum(deathMat==0));
% DeadByFrac = sum(single(infMat==0).*single(deathMat==1))./(sum(deathMat==1)+sum(deathMat==0));
%
% FracHealthyat30 = [FracHealthyat30, survivorFrac(t==timept)];
% FracInfat30 = [FracInfat30, infFrac(t==timept)];
% FracDeadInf30 = [FracDeadInf30, DeadInfFrac(t==timept)];
% FracDeadBys30 = [FracDeadBys30, DeadByFrac(t==timept)];
% %h=semilogy(t, (survivorFrac),'o');h
% %h.Color = tzeva(i,:);
% %hold all
%
% %BETA = R.getData('ExpFitParam',pos);
% % y = Exponent(BETA,t);
% % h1 = semilogy(t, y)
% % h1.Color = tzeva(i,:);
% % h1.LineWidth=1.5
%
% end
% % xlabel('time(h)')
% %set(gca,'XLim',[0,max(t)],'YLim',[0.05,1],'YTick',[0.1 0.25 0.5 0.75 1],'YTickLabel',100*[],'TickLength', [0.03 0.01])
% %title('+HSV-I')
%
% % cb = colorbar('Position', [0.9 0.1400 0.0324 0.7500]);
% % cb.Ticks = asinh(fliplr([100*(1/2).^[0:9], 0]))./asinh(100);
% % cb.TickLabels= num2str(fliplr([100*(1/2).^[0:9], 0])','%3.1f');
% % cb.Label.String = '[TNF\alpha] (ng/ml)'
% % cb.Label.FontSize=12
% % cb.Label.Position(1) = 0.2
% % cb.Label.Color='w'
% nTNF = numel(cTNF);
% ax = axes('Position', [0.15, 0.2, 0.7, 0.7])
% %h1 = scatter(nTNF:-1:1,FracHealthyat30,[],tzeva,'filled');
% %hold on
%
% h1 = plot(asinh(cTNF/LinearRange),100*FracHealthyat30,'-ob')
% h1.MarkerFaceColor = 'b';
% set(gca,'YLim',[-1,100])
% ylabel('Fraction')
% hold on
%
% %h1 = scatter(nTNF:-1:1,FracInfat30,[],tzeva,'filled');
% % hold on
% % h2 = plot(asinh(cTNF/LinearRange),100*FracInfat30,'-or')
% % h2.MarkerFaceColor = 'r'
% %
% %
% % h1 = plot(asinh(cTNF/LinearRange),100*FracDeadBys30,'-ok')
% % h1.MarkerFaceColor = 'k';
% % set(gca,'YLim',[-1,100])
% % xl = xlabel('[TNF] ng/ml');
% %
% %
% %
% % set(gca, 'fontsize', 17, 'xtick', XTicks, 'xticklabel', XTickLabel, 'xlim', [-0.5 8],'YTick',[0:20:100],'YTickLabel',[0:20:100]/100)
% %
% % set(gca,'XLim',asinh([-.1,150]/LinearRange))
% %
% %
% %
% % xl = xlabel('[TNF\alpha] ng/ml')
% % %xl.Position(2) = xl.Position(2)-30
% % title('Experiment')
% %
% % hleg = legend('Healthy', 'Infected', 'Erroneous dead');
% % hleg.Box='off';
% % hleg.Position(1) = hleg.Position(1)-0.3;
% % %%
% % figure(gcf)
% % set(gcf, 'PaperPositionMode','auto','color','w','InvertHardcopy','off')
% % print(gcf,'-depsc','-r600',[savpath 'SweetSpotExperimentalInfectedHealthyFP'])
%
% %% Plot statistics
% LinearRange = 0.1;
% a1 = sort(- (1:9)'*10.^(1:2));
% a2 =(1:9)'*10.^(1:5);
% a = [sort(a1(:)); 0; a2(:)]/1000;
% emptyStr = {''; ''; '';''; ''; '';''; ''};
% XTickLabel = [emptyStr; ''; emptyStr; {''}; '0'; {''}; emptyStr; ' 0.1'; emptyStr; ' 1'; emptyStr; ' 10'; emptyStr; '    100'; emptyStr];
% XTicks = asinh(a/LinearRange);
%
%
% close all
% ax = subplot(2,2,1)
%
% h1 = plot(asinh(cTNF/LinearRange),100*FracHealthyat30,'-ob')
% h1.MarkerFaceColor = 'b';
% set(gca,'YLim',[-1,100])
% ylabel('%Healthy')
% xl = xlabel('[TNF] ng/ml');
% set(gca, 'fontsize', 17, 'xtick', XTicks, 'xticklabel', XTickLabel, 'xlim', [-0.5 8])
%
% ax = subplot(2,2,2)
% h1 = plot(asinh(cTNF/LinearRange),100*FracInfat30,'-or')
% h1.MarkerFaceColor = 'r';
% set(gca,'YLim',[-1,100])
% ylabel('%Infected')
% xl = xlabel('[TNF] ng/ml');
% set(gca, 'fontsize', 17, 'xtick', XTicks, 'xticklabel', XTickLabel, 'xlim', [-0.5 8])
%
% ax = subplot(2,2,3)
% h1 = plot(asinh(cTNF/LinearRange),100*FracDeadInf30,'-og')
% h1.MarkerFaceColor = 'g';
% set(gca,'YLim',[-1,100])
% ylabel('%Dead Infected')
% xl = xlabel('[TNF] ng/ml');
% set(gca, 'fontsize', 17, 'xtick', XTicks, 'xticklabel', XTickLabel, 'xlim', [-0.5 8])
%
% ax = subplot(2,2,4)
% h1 = plot(asinh(cTNF/LinearRange),100*FracDeadBys30,'-ok')
% h1.MarkerFaceColor = 'k';
% set(gca,'YLim',[-1,100])
% ylabel('%Dead bystander')
% xl = xlabel('[TNF] ng/ml');
% set(gca, 'fontsize', 17, 'xtick', XTicks, 'xticklabel', XTickLabel, 'xlim', [-0.5 8])
% %%
% figure(gcf)
% set(gcf, 'PaperPositionMode','auto','color','w','InvertHardcopy','off')
% print(gcf,'-depsc','-r600',[savpath 'SweetSpotExperimentalAllSeparate'])
%
%
%



%% Heatmaps of infection/death high MOI

%% Load results high MOI
BaseStr = regexprep([char(ispc.*'Z:\Images2019\') char(isunix.*'/bigstore/Images2019/')],char(0),'');
Usr = 'Jen';
Project = 'NFkBDynamics';
Dataset = 'TNF_SpeedvAcc_HSVMOI10_Rpt_2019Jun14';
acquisition = 1;
% Get MD of raw data
acqname = ['acq_' num2str(acquisition)];
fpath = [BaseStr Usr filesep Project filesep Dataset filesep acqname];
MD=Metadata(fpath);
Wells = unique(MD.getSpecificMetadata('Position'));
frames = unique(cell2mat(MD.getSpecificMetadata('frame')));
dt = 1/3;%h

cTNF = ([100*(1/2).^[0:8], 0]);
% Load Results

R = MultiPositionSingleCellVirusResults(fpath)
Wells = R.PosNames;
frames = R.Frames;
t = [R.Frames-1]'.*dt;



%%
figure
i=9
pos = R.PosNames{i};
Tracks = R.getTracks(pos)

inds = find(logical(arrayfun(@(x) x.CellDies.*x.CellsGetInfected, Tracks)));
deathMat = R.heatmapTracks(pos, 'DeathTrack');
virMat = R.heatmapTracks(pos, 'VirusTrack');
TDies = arrayfun(@(x) x.T(x.indWhenCellDies), Tracks,'UniformOutput',false);
Tinf = arrayfun(@(x) x.T(x.indWhenCellGetsInfected), Tracks,'UniformOutput',false);


RGB = cat(3,imadjust(deathMat(inds,:),[0.02,0.05]), imadjust(virMat(inds,:), [0.02, 0.05]), imadjust(deathMat(inds,:),[0.02,0.05]));
dT=[];
for j=1:numel(inds)
    RGB(j, Tinf{inds(j)},:)=1;
    RGB(j, TDies{inds(j)},:)=[1 0 0];
    dT = [dT TDies{inds(j)}-Tinf{inds(j)}];
end
histogram(dT,40)


figure

imagesc(RGB); shg






%%
figure
i=1
pos = R.PosNames{i};
Tracks = R.getTracks(pos)

%inds = find(logical(arrayfun(@(x) x.CellDies.*x.CellsGetInfected, Tracks)));
inds = 1:numel(Tracks)
deathMat = R.heatmapTracks(pos, 'DeathTrack');
virMat = R.heatmapTracks(pos, 'VirusTrack');
TDies = arrayfun(@(x) x.T(x.indWhenCellDies), Tracks,'UniformOutput',false);
Tinf = arrayfun(@(x) x.T(x.indWhenCellGetsInfected), Tracks,'UniformOutput',false);


RGB = cat(3,imadjust(virMat(inds,:),[0.02,0.08]), imadjust(virMat(inds,:), [0.02, 0.08]), imadjust(virMat(inds,:),[0.02,0.08]));
dT=[];
for j=1:numel(inds)
    if ~isempty(Tinf{inds(j)})
        RGB(j, Tinf{inds(j)},:)=[0 1 0];
    end
    if ~isempty(TDies{inds(j)})
        RGB(j, TDies{inds(j)},:)=[1 0 0];
    end
    %dT = [dT TDies{inds(j)}-Tinf{inds(j)}];
end
%mean(dT)

histogram(cat(2,Tinf{:})); hold on;
histogram(cat(2,TDies{:}))

figure

imshow(RGB); shg


%%
figure
i=1

subplot(3,1,1)
pos = R.PosNames{i};
Tracks = R.getTracks(pos)

%inds = find(logical(arrayfun(@(x) x.CellDies.*x.CellsGetInfected, Tracks)));
inds = 1:numel(Tracks)
deathMat = R.heatmapTracks(pos, 'DeathTrack');
virMat = R.heatmapTracks(pos, 'VirusTrack');
TDies = arrayfun(@(x) x.T(x.indWhenCellDies), Tracks,'UniformOutput',false);
Tinf = arrayfun(@(x) x.T(x.indWhenCellGetsInfected), Tracks,'UniformOutput',false);


histogram(cat(2,Tinf{:}),20,'normalization','probability'); hold on;
histogram(cat(2,TDies{:}),20,'normalization','probability')

subplot(3,1,2)

i=16
pos = R.PosNames{i};
Tracks = R.getTracks(pos)

%inds = find(logical(arrayfun(@(x) x.CellDies.*x.CellsGetInfected, Tracks)));
inds = 1:numel(Tracks)
deathMat = R.heatmapTracks(pos, 'DeathTrack');
virMat = R.heatmapTracks(pos, 'VirusTrack');
TDies = arrayfun(@(x) x.T(x.indWhenCellDies), Tracks,'UniformOutput',false);
Tinf = arrayfun(@(x) x.T(x.indWhenCellGetsInfected), Tracks,'UniformOutput',false);


histogram(cat(2,Tinf{:}),20,'normalization','probability'); hold on;
histogram(cat(2,TDies{:}),20,'normalization','probability')


subplot(3,1,3)

i=10
pos = R.PosNames{i};
Tracks = R.getTracks(pos)

%inds = find(logical(arrayfun(@(x) x.CellDies.*x.CellsGetInfected, Tracks)));
inds = 1:numel(Tracks)
deathMat = R.heatmapTracks(pos, 'DeathTrack');
virMat = R.heatmapTracks(pos, 'VirusTrack');
TDies = arrayfun(@(x) x.T(x.indWhenCellDies), Tracks,'UniformOutput',false);
Tinf = arrayfun(@(x) x.T(x.indWhenCellGetsInfected), Tracks,'UniformOutput',false);


histogram(cat(2,Tinf{:}),20,'normalization','probability'); hold on;
histogram(cat(2,TDies{:}),20,'normalization','probability')


%%
for i=1:20
    pos = R.PosNames{i};
    Tracks = R.getTracks(pos);
    
    inds = find(logical(arrayfun(@(x) x.CellDies.*x.CellsGetInfected, Tracks)));
    TDies = arrayfun(@(x) x.T(x.indWhenCellDies), Tracks,'UniformOutput',false);
    
    dT=[];
    for j=1:numel(inds)
        dT = [dT TDies{inds(j)}];
    end
    
    mean(dT)
end









%%%%%%%%%%
%%%%%%%%%%
%%%%%%%%%%
%% Fig4: Spim stuff
%% Load MD
BaseStr = '/RazorScopeData/RazorScopeSets/';
Usr = 'Alon';
Project = 'CorneaHSV';
Dataset = 'Infection48hHwithOUT_TNF_Region_2019Nov19';
acquisition = 3;
% Get MD of raw data
acqname = ['acq_' num2str(acquisition)];
fpath = [BaseStr Usr filesep Project filesep Dataset filesep acqname ];
ImageJ;
%% Get metadata and load results object
SpimMD = SpimMetadata(fpath)
R = MulticolorSPIMResults(fpath)

tpToPlot = 49;
slices = 185:220;

D = R.SPIMTimepoints{tpToPlot}.BlockLbls{1}.stkshow('SpimMD',SpimMD,'channel','Red','maxproj',true,'slices', slices,'refTP', tpToPlot,'title','Red');
ij.IJ.setMinAndMax(2^16*prctile(D(D~=0),5), 2^16*prctile(D(D~=0),99.99))
ij.IJ.run('Invert LUT');

D = R.SPIMTimepoints{tpToPlot}.BlockLbls{1}.stkshow('SpimMD',SpimMD,'channel','DeepBlue','maxproj',true,'slices', slices,'refTP', tpToPlot,'title','DeepBlue');
ij.IJ.setMinAndMax(2^16*prctile(D(D~=0),5), 2^16*prctile(D(D~=0),99.9))
ij.IJ.run('purples')

ij.IJ.run('Merge Channels...', 'c4=Red c2=DeepBlue create');
ij.IJ.run('Properties...', 'unit=um pixel_width=0.7606 pixel_height=0.7606 voxel_depth=3');
ij.IJ.run('Scale Bar...', 'width=100 height=8 font=28 color=Black background=None location=[Lower Right] hide overlay');

ij.IJ.run('Flatten');
%ij.IJ.saveAs('Tiff', '/bigstore/GeneralStorage/Alon/Figures/DecisionPaper2019/Figures010420/SPIMStuff/Fig4_24hNoTNFVirusOnly.tif');


%% Load MD
BaseStr = '/RazorScopeData/RazorScopeSets/';
Usr = 'Alon';
Project = 'CorneaHSV';
Dataset = 'Infection48hHwith_TNF_Region_2019Dec11';
acquisition = 2;
% Get MD of raw data
acqname = ['acq_' num2str(acquisition)];
fpath = [BaseStr Usr filesep Project filesep Dataset filesep acqname ];
% Get metadata and load results object
%MD=Metadata(fpath);
SpimMD = SpimMetadata(fpath)
R = MulticolorSPIMResults(fpath)

%
tpToPlot = 49;
slices = 195:225;

D = R.SPIMTimepoints{tpToPlot}.BlockLbls{4}.stkshow('SpimMD',SpimMD,'channel','Red','maxproj',true,'slices', slices,'refTP', tpToPlot,'title','Red');
ij.IJ.setMinAndMax(2^16*prctile(D(D~=0),5), 2^16*prctile(D(D~=0),99.99))
ij.IJ.run('Invert LUT');

D = R.SPIMTimepoints{tpToPlot}.BlockLbls{4}.stkshow('SpimMD',SpimMD,'channel','DeepBlue','maxproj',true,'slices', slices,'refTP', tpToPlot,'title','DeepBlue');
ij.IJ.setMinAndMax(2^16*prctile(D(D~=0),5), 2^16*prctile(D(D~=0),99.9))
ij.IJ.run('purples')

ij.IJ.run('Merge Channels...', 'c4=Red c2=DeepBlue create');
ij.IJ.run('Properties...', 'unit=um pixel_width=0.7606 pixel_height=0.7606 voxel_depth=3');

ij.IJ.run('Scale Bar...', 'width=100 height=8 font=28 color=Black background=None location=[Lower Right] hide overlay');
ij.IJ.run('Flatten');
%ij.IJ.saveAs('Tiff', '/bigstore/GeneralStorage/Alon/Figures/DecisionPaper2019/Figures010420/SPIMStuff/Fig4_24hHighTNFVirusOnly.tif');



%%
fpath = '/RazorScopeData/RazorScopeImages/Alon/CorneaHSV/SingleTP48h_2019Dec14'
MD = Metadata(fpath)
%%

rng = [2 3 6 8]
hists = cell(4,1);
edges = cell(4,1);
totalVirus = []
for i=1:numel(rng)    
Red = stkread(MD,'acq', ['acq_' num2str(rng(i))],'Channel','Red');
Green = stkread(MD,'acq', ['acq_' num2str(rng(i))],'Channel','Green');

[hists{i}, edges{i}] = histcounts(log(Green(Red>.0025)),200);
totalVirus = [totalVirus sum(sum(sum(Red(Red>.0025))))];
end

%totalVirus = cellfun(@sum, hists);
totalVirus = reshape(totalVirus,2,[]);
fracDead = cellfun(@(x) sum(x(140:end))./sum(x), hists);
fracDead = reshape(fracDead,2,[]);
%% Panel total virus
close all
figure('color','w','Position',[100,100, 300, 300])
ax = axes('Position', [0.24, 0.18, 0.7, 0.7])

%[hBar , hErr] = barwitherr(std(totalVirus), mean(totalVirus))

hBar = bar([1,2],mean(totalVirus))
hold on
hErr = errorbar([1,2],mean(totalVirus),std(totalVirus),'.')

hBar.FaceColor = virusColor
hErr.LineWidth = 2;
hErr.Color = 'k';

%stupid matlab
mult = 0;                               % 5 times as long
b = hErr.Bar;                           % hidden property/handle
drawnow                                 % populate b's properties
vd = b.VertexData;
N = numel(hErr.XData);                           % number of error bars
capLength = vd(1,2*N+2,1) - vd(1,1,1);  % assumes equal length on all
newLength = capLength * mult;
leftInds = N*2+1:2:N*6;
rightInds = N*2+2:2:N*6;
vd(1,leftInds,1) = [hErr.XData-newLength, hErr.XData-newLength];
vd(1,rightInds,1) = [hErr.XData+newLength, hErr.XData+newLength];
b.VertexData = vd;

ylabel('Total virus (a.u.)')
xlabel('[TNF] ng/ml')

ax.XTickLabel = {'0', '50'}
ax.YLim = [0 1.5*10^4]
ax.YTick = 0:5000:20000;
ax.YTickLabel = 0:5:20;
ax.XLim = [0.5 2.5]

shg

%%
figname = 'Fig4_SingleTP_TotalVirus'
set(gcf, 'PaperPositionMode','auto','color','w','InvertHardcopy','off')
print(gcf,'-depsc','-r600',[savpath '/epss/' figname]);
print(gcf,'-dpng','-r600',[savpath '/pngs/' figname]);

%% Panel total virus
close all
figure('color','w','Position',[100,100, 300, 300])
ax = axes('Position', [0.24, 0.18, 0.7, 0.7])

%[hBar , hErr] = barwitherr(std(totalVirus), mean(totalVirus))

hBar = bar([1,2],mean(totalVirus))
hold on
%hErr = errorbar([1,2],mean(totalVirus),std(totalVirus),'.')

hBar.FaceColor = virusColor
%hErr.LineWidth = 2;
%hErr.Color = 'k';

h2=scatter([1,1,2,2],totalVirus(:),50,'r');
h2.MarkerFaceColor='r'


%stupid matlab
mult = 0;                               % 5 times as long
%b = hErr.Bar;                           % hidden property/handle
% drawnow                                 % populate b's properties
% vd = b.VertexData;
% N = numel(hErr.XData);                           % number of error bars
% capLength = vd(1,2*N+2,1) - vd(1,1,1);  % assumes equal length on all
% newLength = capLength * mult;
% leftInds = N*2+1:2:N*6;
% rightInds = N*2+2:2:N*6;
% vd(1,leftInds,1) = [hErr.XData-newLength, hErr.XData-newLength];
% vd(1,rightInds,1) = [hErr.XData+newLength, hErr.XData+newLength];
% b.VertexData = vd;

ylabel('Total virus (a.u.)')
xlabel('[TNF] ng/ml')

ax.XTickLabel = {'0', '50'}
ax.YLim = [0 1.5*10^4]
ax.YTick = 0:5000:20000;
ax.YTickLabel = 0:5:20;
ax.XLim = [0.5 2.5]

shg

%%
figname = 'Fig4_SingleTP_TotalVirus_wIndPoints'
set(gcf, 'PaperPositionMode','auto','color','w','InvertHardcopy','off')
print(gcf,'-depsc','-r600',[savpath '/epss/' figname]);
print(gcf,'-dpng','-r600',[savpath '/pngs/' figname]);
%% Panel percent dead virus
figure('color','w','Position',[100,100, 300, 300])
ax = axes('Position', [0.24, 0.18, 0.7, 0.7])

%[hBar , hErr] = barwitherr(std(totalVirus), mean(totalVirus))

hBar = bar([1,2],mean(fracDead))
hold on
hErr = errorbar([1,2],mean(fracDead),std(fracDead),'.')

hBar.FaceColor = deathColor
hErr.LineWidth = 2;
hErr.Color = 'k';

%stupid matlab
mult = 0;                               % 5 times as long
b = hErr.Bar;                           % hidden property/handle
drawnow                                 % populate b's properties
vd = b.VertexData;
N = numel(hErr.XData);                           % number of error bars
capLength = vd(1,2*N+2,1) - vd(1,1,1);  % assumes equal length on all
newLength = capLength * mult;
leftInds = N*2+1:2:N*6;
rightInds = N*2+2:2:N*6;
vd(1,leftInds,1) = [hErr.XData-newLength, hErr.XData-newLength];
vd(1,rightInds,1) = [hErr.XData+newLength, hErr.XData+newLength];
b.VertexData = vd;
ylabel('Percent dead of infected')
xlabel('[TNF] ng/ml')

ax.XTickLabel = {'0', '50'}
ax.YLim = [0 1]
ax.YTick = 0:0.5:1;
ax.XLim = [0.5 2.5]

shg
%%
figname = 'Fig4_SingleTP_PercentDeadOfInfected'
set(gcf, 'PaperPositionMode','auto','color','w','InvertHardcopy','off')
print(gcf,'-depsc','-r600',[savpath '/epss/' figname]);
print(gcf,'-dpng','-r600',[savpath '/pngs/' figname]);

%% Panel percent dead virus
figure('color','w','Position',[100,100, 300, 300])
ax = axes('Position', [0.24, 0.18, 0.7, 0.7])

%[hBar , hErr] = barwitherr(std(totalVirus), mean(totalVirus))

hBar = bar([1,2],mean(fracDead))
hold on
%hErr = errorbar([1,2],mean(fracDead),std(fracDead),'.')

hBar.FaceColor = [0.6 0.6 0.6]
%hErr.LineWidth = 2;
%hErr.Color = 'k';
h2=scatter([1,1,2,2],fracDead(:),50,'r');
h2.MarkerFaceColor='r'

% %stupid matlab
% mult = 0;                               % 5 times as long
% b = hErr.Bar;                           % hidden property/handle
% drawnow                                 % populate b's properties
% vd = b.VertexData;
% N = numel(hErr.XData);                           % number of error bars
% capLength = vd(1,2*N+2,1) - vd(1,1,1);  % assumes equal length on all
% newLength = capLength * mult;
% leftInds = N*2+1:2:N*6;
% rightInds = N*2+2:2:N*6;
% vd(1,leftInds,1) = [hErr.XData-newLength, hErr.XData-newLength];
% vd(1,rightInds,1) = [hErr.XData+newLength, hErr.XData+newLength];
% b.VertexData = vd;
 ylabel('Percent dead of infected')
 xlabel('[TNF] ng/ml')

ax.XTickLabel = {'0', '50'}
ax.YLim = [0 1]
ax.YTick = 0:0.5:1;
ax.YTickLabel=100*[0:0.5:1];
ax.XLim = [0.5 2.5]

shg

%%
figname = 'Fig4_SingleTP_PercentDeadOfInfected_withPoints'
set(gcf, 'PaperPositionMode','auto','color','w','InvertHardcopy','off')
print(gcf,'-depsc','-r600',[savpath '/epss/' figname]);
print(gcf,'-dpng','-r600',[savpath '/pngs/' figname]);

%% Panel frac dead
[hBar , hErr] = barwitherr(std(fracDead), mean(fracDead))


%%
a = []
for i=1:97; a = [a sum(WellCells{i}.BlockLbls{4}.Int90Prctile{3}(find(WellCells{i}.BlockLbls{4}.Int90Prctile{2}>0.08))>0.2)]; end
plot(Smoothing(a,'neigh',10)); shg
%
hold on
a = [];
for i=1:97; a = [a sum(WellCells{i}.BlockLbls{4}.Int90Prctile{2}>0.08)]; end
plot(Smoothing(a,'neigh',10)); shg




%% Load MD
BaseStr = '/RazorScopeData/RazorScopeSets/';
Usr = 'Alon';
Project = 'CorneaHSV';
Dataset = 'Infection48hHwithOUT_TNF_Region_2019Nov19';
acquisition = 3;
% Get MD of raw data
acqname = ['acq_' num2str(acquisition)];
fpath = [BaseStr Usr filesep Project filesep Dataset filesep acqname ];
SpimMD = SpimMetadata(fpath)
R = MulticolorSPIMResults(fpath)
%% Get metadata and load results object


refTP = 49
slices = 185:220;

Nuclei = cell(numel(R.Frames),1);
Dead = cell(numel(R.Frames),1);
Virus = cell(numel(R.Frames),1);
%%
parfor tpToPlot = 1:numel(R.Frames);
SpimMD = SpimMetadata(fpath);
%D = R.SPIMTimepoints{tpToPlot}.BlockLbls{1}.stkread('SpimMD',SpimMD,'channel','Red','maxproj',true,'slices', slices,'refTP', refTP,'title','Red');
D = SpimMD.stkread(tpToPlot,0,'Red','maxproj',true,'slices', slices,'refTP', refTP);

%ij.IJ.setMinAndMax(2^16*prctile(D(D~=0),5), 2^16*prctile(D(D~=0),99.99))
%ij.IJ.run('Invert LUT');
Nuclei{tpToPlot} = D;

%D = R.SPIMTimepoints{tpToPlot}.BlockLbls{1}.stkread('SpimMD',SpimMD,'channel','DeepBlue','maxproj',true,'slices', slices,'refTP', refTP,'title','DeepBlue');
D = SpimMD.stkread(tpToPlot,0,'DeepBlue','maxproj',true,'slices', slices,'refTP', refTP);
%ij.IJ.setMinAndMax(2^16*prctile(D(D~=0),5), 2^16*prctile(D(D~=0),99.9))
%ij.IJ.run('purples')
Virus{tpToPlot} = D;

%D = R.SPIMTimepoints{tpToPlot}.BlockLbls{1}.stkread('SpimMD',SpimMD,'channel','Green','maxproj',true,'slices', slices,'refTP', refTP,'title','Green');
D = SpimMD.stkread(tpToPlot,0,'Green','maxproj',true,'slices', slices,'refTP', refTP);
%ij.IJ.setMinAndMax(2^16*prctile(D(D~=0),5), 2^16*prctile(D(D~=0),99.9))
%ij.IJ.run('purples')
Dead{tpToPlot} = D;

end

%%
ij.IJ.run('Merge Channels...', 'c4=Red c2=DeepBlue create');
ij.IJ.run('Properties...', 'unit=um pixel_width=0.7606 pixel_height=0.7606 voxel_depth=3');
ij.IJ.run('Scale Bar...', 'width=100 height=8 font=28 color=Black background=None location=[Lower Right] hide overlay');

ij.IJ.run('Flatten');
ij.IJ.saveAs('Tiff', '/bigstore/GeneralStorage/Alon/Figures/DecisionPaper2019/Figures010420/SPIMStuff/Fig4_24hNoTNFVirusOnly.tif');



%% Figure 4, quantify corneas bystander deatn w/wo TNF

MD = Metadata('/bigstore/Images2019/Jen/CorneaCCM/TNF_Death_T1_24h_2019Dec13/acq_1/ActualSAR');
Data = MD.stkread('Channel','Yellow'); %load all Yellow
        
 TotalDeath = fliplr(reshape(sum(reshape(squeeze(sum(sum(Data))),[],4)),[],2));
 deathColor = [255,165,0]/255; %orange

 close all
figure('color','w','Position',[100,100, 300, 300])
ax = axes('Position', [0.24, 0.18, 0.7, 0.7])

%[hBar , hErr] = barwitherr(std(totalVirus), mean(totalVirus))

hBar = bar([1,2],mean(TotalDeath))
hold on
hErr = errorbar([1,2],mean(TotalDeath),std(TotalDeath),'.')

hBar.FaceColor = deathColor
hErr.LineWidth = 2;
hErr.Color = 'k';

%stupid matlab
mult = 0;                               % 5 times as long
b = hErr.Bar;                           % hidden property/handle
drawnow                                 % populate b's properties
vd = b.VertexData;
N = numel(hErr.XData);                           % number of error bars
capLength = vd(1,2*N+2,1) - vd(1,1,1);  % assumes equal length on all
newLength = capLength * mult;
leftInds = N*2+1:2:N*6;
rightInds = N*2+2:2:N*6;
vd(1,leftInds,1) = [hErr.XData-newLength, hErr.XData-newLength];
vd(1,rightInds,1) = [hErr.XData+newLength, hErr.XData+newLength];
b.VertexData = vd;

ylabel('Total death (a.u.)')
xlabel('[TNF] ng/ml')

ax.XTickLabel = {'0', '50'}
ax.YLim = [0 1.7*10^5]
ax.YTick = 0:50000:200000;
ax.YTickLabel = 0:5:20;
ax.XLim = [0.5 2.5]

shg

%%
figname = 'Fig4_CorneaDeathW_WO_TNF'
set(gcf, 'PaperPositionMode','auto','color','w','InvertHardcopy','off')
print(gcf,'-depsc','-r600',[savpath '/epss/' figname]);
print(gcf,'-dpng','-r600',[savpath '/pngs/' figname]);



%% Figure 4, quantify corneas bystander deatn w/wo TNF

MD = Metadata('/bigstore/Images2019/Jen/CorneaCCM/TNF_Death_T1_24h_2019Dec13/acq_1/ActualSAR');
Data = MD.stkread('Channel','Yellow'); %load all Yellow
        
totdeath=squeeze(sum(sum(Data)));
 TotalDeath = fliplr(reshape(sum(reshape(totdeath,[],4)),[],2));
 deathColor = [255,165,0]/255; %orange

 close all
figure('color','w','Position',[100,100, 300, 300])
ax = axes('Position', [0.24, 0.18, 0.7, 0.7])

%[hBar , hErr] = barwitherr(std(totalVirus), mean(totalVirus))

hBar = bar([1,2],mean(TotalDeath))
hold on
%hErr = errorbar([1,2],mean(TotalDeath),std(TotalDeath),'.')

hBar.FaceColor = deathColor
%hErr.LineWidth = 2;
%hErr.Color = 'k';

h2=scatter([1,1,2,2],TotalDeath(:),50,'r');
h2.MarkerFaceColor='r'
%stupid matlab
% mult = 0;                               % 5 times as long
% b = hErr.Bar;                           % hidden property/handle
% drawnow                                 % populate b's properties
% vd = b.VertexData;
% N = numel(hErr.XData);                           % number of error bars
% capLength = vd(1,2*N+2,1) - vd(1,1,1);  % assumes equal length on all
% newLength = capLength * mult;
% leftInds = N*2+1:2:N*6;
% rightInds = N*2+2:2:N*6;
% vd(1,leftInds,1) = [hErr.XData-newLength, hErr.XData-newLength];
% vd(1,rightInds,1) = [hErr.XData+newLength, hErr.XData+newLength];
% b.VertexData = vd;

ylabel('Total death (a.u.)')
xlabel('[TNF] ng/ml')

ax.XTickLabel = {'0', '50'}
ax.YLim = [0 1.7*10^5]
ax.YTick = 0:50000:200000;
ax.YTickLabel = 0:5:20;
ax.XLim = [0.5 2.5]

shg

%%
figname = 'Fig4_CorneaDeathW_WO_TNF_withPoints'
set(gcf, 'PaperPositionMode','auto','color','w','InvertHardcopy','off')
print(gcf,'-depsc','-r600',[savpath '/epss/' figname]);
print(gcf,'-dpng','-r600',[savpath '/pngs/' figname]);



%% Load MD 
BaseStr = '/RazorScopeData/RazorScopeSets/';
Usr = 'Alon';
Project = 'CorneaHSV';
Dataset = 'Infection48hHwithOUT_TNF_Region_2019Nov19';
acquisition = 3;
% Get MD of raw data
acqname = ['acq_' num2str(acquisition)];
fpath = [BaseStr Usr filesep Project filesep Dataset filesep acqname ];
% Get metadata

MD=Metadata(fpath);
% Make imageTP
dirname = [MD.pth '/Projections/'];

tp = 73;

ij.IJ.run('Close All')
pause(0.5)
flist = getFlistbyPattern('Red',dirname);
flist = flist(tp);
D = []
for i=1:numel(flist)
    D = cat(3,D,imread([dirname flist{i}]));
end
D = single(D)./2^16;
stkshow(D,'title', 'Red_1')
ij.IJ.run('Duplicate...', 'title=Red duplicate');
ij.IJ.setMinAndMax(2^16*prctile(D(D~=0),10), 2^16*prctile(D(D~=0),99.9))
ij.IJ.run('Invert LUT')
ij.IJ.selectWindow('Red_1');
ij.IJ.run('Close');
pause(0.1)

flist = getFlistbyPattern('DeepBlue',dirname);
flist = flist(tp);
D = []
for i=1:numel(flist)
    D = cat(3,D,imread([dirname flist{i}]));
end
D = single(D)./2^16;
stkshow(D,'title', 'DeepBlue_1')
ij.IJ.run('Duplicate...', 'title=DeepBlue duplicate');
ij.IJ.setMinAndMax(4800, 2^16*prctile(D(D~=0),99.9)); %2^16*prctile(D(D~=0),90)v
ij.IJ.run('purples')
ij.IJ.run('Apply LUT')
ij.IJ.selectWindow('DeepBlue_1');
ij.IJ.run('Close');
pause(0.1)

flist = getFlistbyPattern('Green',dirname);
flist = flist(tp);
D = []
for i=1:numel(flist)
    D = cat(3,D,imread([dirname flist{i}]));
end
D = single(D)./2^16;
stkshow(D,'title', 'Green_1')
ij.IJ.run('Duplicate...', 'title=Green duplicate');
ij.IJ.setMinAndMax(12800, 2^16*prctile(D(D~=0),99.9))
ij.IJ.run('oranges')
ij.IJ.run('Apply LUT')

ij.IJ.selectWindow('Green_1');
ij.IJ.run('Close');
pause(0.1)

ij.IJ.run('Merge Channels...', 'c2=DeepBlue c3=Green create');
ij.IJ.selectWindow('Composite');

ij.IJ.run('Stack to RGB', 'slices');

ij.IJ.selectWindow('Composite');
ij.IJ.run('Close');
pause(0.1)

ij.IJ.selectWindow('Composite (RGB)');
% overlay
for i=1:size(D,3)
    ij.IJ.selectWindow('Red');
    ij.IJ.setSlice(i);
    ij.IJ.selectWindow('Composite (RGB)');
    ij.IJ.setSlice(i);
    ij.IJ.selectWindow('Red');
    ij.IJ.run('Add Image...', 'image=[Composite (RGB)] x=0 y=0 opacity=100 zero');
end
% set scale and flatten

ij.IJ.selectWindow('Red');
ij.IJ.run('Properties...', 'unit=um pixel_width=0.7606 pixel_height=0.7606 voxel_depth=3');
ij.IJ.run('Scale Bar...', 'width=100 height=4 font=28 color=black background=None location=[Lower Right] hide overlay');
ij.IJ.setForegroundColor(0, 0, 0);
ij.IJ.setBackgroundColor(0, 0, 0);
ij.IJ.run('Label...', ['format=0 starting=36 interval=0.5 x=10 y=50 font=30 text=[h] range=1-' num2str(size(D,3))]);

ij.IJ.run('Flatten','stack');
%%
ij.IJ.saveAs('Tiff', [savpath '/SPIMStuff/Fig4Spim_Trans_NoTNF.tif']);



%% Load MD 
BaseStr = '/RazorScopeData/RazorScopeSets/';
Usr = 'Alon';
Project = 'CorneaHSV';
Dataset = 'Infection48hHwithOUT_TNF_Region_2019Nov19';
acquisition = 3;
% Get MD of raw data
acqname = ['acq_' num2str(acquisition)];
fpath = [BaseStr Usr filesep Project filesep Dataset filesep acqname ];
% Get metadata

MD=Metadata(fpath);
% Make imageTP
dirname = [MD.pth '/Projections/'];

tp = 73;

ij.IJ.run('Close All')
pause(0.5)
flist = getFlistbyPattern('Red',dirname);
flist = flist(tp);
D = []
for i=1:numel(flist)
    D = cat(3,D,imread([dirname flist{i}]));
end
D = single(D)./2^16;
stkshow(D,'title', 'Red_1')
ij.IJ.run('Duplicate...', 'title=Red duplicate');
ij.IJ.setMinAndMax(2^16*prctile(D(D~=0),10), 2^16*prctile(D(D~=0),99.9))
ij.IJ.run('Invert LUT')
ij.IJ.selectWindow('Red_1');
ij.IJ.run('Close');
pause(0.1)

flist = getFlistbyPattern('DeepBlue',dirname);
flist = flist(tp);
D = []
for i=1:numel(flist)
    D = cat(3,D,imread([dirname flist{i}]));
end
D = single(D)./2^16;
stkshow(D,'title', 'DeepBlue_1')
ij.IJ.run('Duplicate...', 'title=DeepBlue duplicate');
ij.IJ.setMinAndMax(4800, 2^16*prctile(D(D~=0),99.9)); %2^16*prctile(D(D~=0),90)v
ij.IJ.run('purples')
ij.IJ.run('Apply LUT')
ij.IJ.selectWindow('DeepBlue_1');
ij.IJ.run('Close');
pause(0.1)

% flist = getFlistbyPattern('Green',dirname);
% flist = flist(tp);
% D = []
% for i=1:numel(flist)
%     D = cat(3,D,imread([dirname flist{i}]));
% end
% D = single(D)./2^16;
% stkshow(D,'title', 'Green_1')
% ij.IJ.run('Duplicate...', 'title=Green duplicate');
% ij.IJ.setMinAndMax(12800, 2^16*prctile(D(D~=0),99.9))
% ij.IJ.run('oranges')
% ij.IJ.run('Apply LUT')
% 
% ij.IJ.selectWindow('Green_1');
% ij.IJ.run('Close');
% pause(0.1)
% 
% ij.IJ.run('Merge Channels...', 'c2=DeepBlue c3=Green create');
ij.IJ.selectWindow('DeepBlue');

%ij.IJ.run('Stack to RGB', 'slices');

%ij.IJ.selectWindow('DeepBlue');
%ij.IJ.run('Close');
%pause(0.1)

%ij.IJ.selectWindow('DeepBlue (RGB)');
% overlay
for i=1:size(D,3)
    ij.IJ.selectWindow('Red');
    ij.IJ.setSlice(i);
    ij.IJ.selectWindow('DeepBlue');
    ij.IJ.setSlice(i);
    ij.IJ.selectWindow('Red');
    ij.IJ.run('Add Image...', 'image=[DeepBlue] x=0 y=0 opacity=100 zero');
end
% set scale and flatten

ij.IJ.selectWindow('Red');
ij.IJ.run('Properties...', 'unit=um pixel_width=0.7606 pixel_height=0.7606 voxel_depth=3');
ij.IJ.run('Scale Bar...', 'width=100 height=4 font=28 color=black background=None location=[Lower Right] hide overlay');
ij.IJ.setForegroundColor(0, 0, 0);
ij.IJ.setBackgroundColor(0, 0, 0);
ij.IJ.run('Label...', ['format=0 starting=36 interval=0.5 x=10 y=50 font=30 text=[h] range=1-' num2str(size(D,3))]);

ij.IJ.run('Flatten','stack');
%%
ij.IJ.saveAs('Tiff', [savpath '/SPIMStuff/Fig4Spim_Trans_NoTNF_NoDeath.tif']);



%% Load MD
BaseStr = '/RazorScopeData/RazorScopeSets/';
Usr = 'Alon';
Project = 'CorneaHSV';
Dataset = 'Infection48hHwith_TNF_Region_2019Dec11';
acquisition = 2;

% Get MD of raw data
acqname = ['acq_' num2str(acquisition)];
fpath = [BaseStr Usr filesep Project filesep Dataset filesep acqname ];
% Get metadata

MD=Metadata(fpath);
% Make imageTP
dirname = [MD.pth '/Projections/'];

tp = 73;

ij.IJ.run('Close All')
pause(0.5)
flist = getFlistbyPattern('Red',dirname);
flist = flist(tp);
D = []
for i=1:numel(flist)
    D = cat(3,D,imread([dirname flist{i}]));
end
D = single(D)./2^16;
stkshow(D,'title', 'Red_1')
ij.IJ.run('Duplicate...', 'title=Red duplicate');
ij.IJ.setMinAndMax(2^16*prctile(D(D~=0),10), 2^16*prctile(D(D~=0),99.9))
ij.IJ.run('Invert LUT')
ij.IJ.selectWindow('Red_1');
ij.IJ.run('Close');
pause(0.1)

flist = getFlistbyPattern('DeepBlue',dirname);
flist = flist(tp);
D = []
for i=1:numel(flist)
    D = cat(3,D,imread([dirname flist{i}]));
end
D = single(D)./2^16;
stkshow(D,'title', 'DeepBlue_1')
ij.IJ.run('Duplicate...', 'title=DeepBlue duplicate');
ij.IJ.setMinAndMax(4800, 2^16*prctile(D(D~=0),99.9)); %2^16*prctile(D(D~=0),90)v
ij.IJ.run('purples')
ij.IJ.run('Apply LUT')
ij.IJ.selectWindow('DeepBlue_1');
ij.IJ.run('Close');
pause(0.1)

flist = getFlistbyPattern('Green',dirname);
flist = flist(tp);
D = []
for i=1:numel(flist)
    D = cat(3,D,imread([dirname flist{i}]));
end
D = single(D)./2^16;
stkshow(D,'title', 'Green_1')
ij.IJ.run('Duplicate...', 'title=Green duplicate');
ij.IJ.setMinAndMax(12800, 2^16*prctile(D(D~=0),99.9))
ij.IJ.run('oranges')
ij.IJ.run('Apply LUT')

ij.IJ.selectWindow('Green_1');
ij.IJ.run('Close');
pause(0.1)

ij.IJ.run('Merge Channels...', 'c2=DeepBlue c3=Green create');
ij.IJ.selectWindow('Composite');

ij.IJ.run('Stack to RGB', 'slices');

ij.IJ.selectWindow('Composite');
ij.IJ.run('Close');
pause(0.1)

ij.IJ.selectWindow('Composite (RGB)');
% overlay
for i=1:size(D,3)
    ij.IJ.selectWindow('Red');
    ij.IJ.setSlice(i);
    ij.IJ.selectWindow('Composite (RGB)');
    ij.IJ.setSlice(i);
    ij.IJ.selectWindow('Red');
    ij.IJ.run('Add Image...', 'image=[Composite (RGB)] x=0 y=0 opacity=100 zero');
end
% set scale and flatten

ij.IJ.selectWindow('Red');
ij.IJ.run('Properties...', 'unit=um pixel_width=0.7606 pixel_height=0.7606 voxel_depth=3');
ij.IJ.run('Scale Bar...', 'width=100 height=4 font=28 color=black background=None location=[Lower Right] hide overlay');
ij.IJ.setForegroundColor(0, 0, 0);
ij.IJ.setBackgroundColor(0, 0, 0);
ij.IJ.run('Label...', ['format=0 starting=36 interval=0.5 x=10 y=50 font=30 text=[h] range=1-' num2str(size(D,3))]);

ij.IJ.run('Flatten','stack');
%%
ij.IJ.saveAs('Tiff', [savpath '/SPIMStuff/Fig4Spim_Trans_WithTNF.tif']);



%% Load MD 
BaseStr = '/RazorScopeData/RazorScopeSets/';
Usr = 'Alon';
Project = 'CorneaHSV';
Dataset = 'Infection48hHwith_TNF_Region_2019Dec11';
acquisition = 2;
% Get MD of raw data
acqname = ['acq_' num2str(acquisition)];
fpath = [BaseStr Usr filesep Project filesep Dataset filesep acqname ];
% Get metadata

MD=Metadata(fpath);
% Make imageTP
dirname = [MD.pth '/Projections/'];

tp = 73;

ij.IJ.run('Close All')
pause(0.5)
flist = getFlistbyPattern('Red',dirname);
flist = flist(tp);
D = []
for i=1:numel(flist)
    D = cat(3,D,imread([dirname flist{i}]));
end
D = single(D)./2^16;
stkshow(D,'title', 'Red_1')
ij.IJ.run('Duplicate...', 'title=Red duplicate');
ij.IJ.setMinAndMax(2^16*prctile(D(D~=0),10), 2^16*prctile(D(D~=0),99.9))
ij.IJ.run('Invert LUT')
ij.IJ.selectWindow('Red_1');
ij.IJ.run('Close');
pause(0.1)

flist = getFlistbyPattern('DeepBlue',dirname);
flist = flist(tp);
D = []
for i=1:numel(flist)
    D = cat(3,D,imread([dirname flist{i}]));
end
D = single(D)./2^16;
stkshow(D,'title', 'DeepBlue_1')
ij.IJ.run('Duplicate...', 'title=DeepBlue duplicate');
ij.IJ.setMinAndMax(4800, 2^16*prctile(D(D~=0),99.9)); %2^16*prctile(D(D~=0),90)v
ij.IJ.run('purples')
ij.IJ.run('Apply LUT')
ij.IJ.selectWindow('DeepBlue_1');
ij.IJ.run('Close');
pause(0.1)

% flist = getFlistbyPattern('Green',dirname);
% flist = flist(tp);
% D = []
% for i=1:numel(flist)
%     D = cat(3,D,imread([dirname flist{i}]));
% end
% D = single(D)./2^16;
% stkshow(D,'title', 'Green_1')
% ij.IJ.run('Duplicate...', 'title=Green duplicate');
% ij.IJ.setMinAndMax(12800, 2^16*prctile(D(D~=0),99.9))
% ij.IJ.run('oranges')
% ij.IJ.run('Apply LUT')
% 
% ij.IJ.selectWindow('Green_1');
% ij.IJ.run('Close');
% pause(0.1)
% 
% ij.IJ.run('Merge Channels...', 'c2=DeepBlue c3=Green create');
ij.IJ.selectWindow('DeepBlue');

%ij.IJ.run('Stack to RGB', 'slices');

%ij.IJ.selectWindow('DeepBlue');
%ij.IJ.run('Close');
%pause(0.1)

%ij.IJ.selectWindow('DeepBlue (RGB)');
% overlay
for i=1:size(D,3)
    ij.IJ.selectWindow('Red');
    ij.IJ.setSlice(i);
    ij.IJ.selectWindow('DeepBlue');
    ij.IJ.setSlice(i);
    ij.IJ.selectWindow('Red');
    ij.IJ.run('Add Image...', 'image=[DeepBlue] x=0 y=0 opacity=100 zero');
end
% set scale and flatten

ij.IJ.selectWindow('Red');
ij.IJ.run('Properties...', 'unit=um pixel_width=0.7606 pixel_height=0.7606 voxel_depth=3');
ij.IJ.run('Scale Bar...', 'width=100 height=4 font=28 color=black background=None location=[Lower Right] hide overlay');
ij.IJ.setForegroundColor(0, 0, 0);
ij.IJ.setBackgroundColor(0, 0, 0);
ij.IJ.run('Label...', ['format=0 starting=36 interval=0.5 x=10 y=50 font=30 text=[h] range=1-' num2str(size(D,3))]);

ij.IJ.run('Flatten','stack');
%%
ij.IJ.saveAs('Tiff', [savpath '/SPIMStuff/Fig4Spim_Trans_WithTNF_NoDeath.tif']);




%% Colorbars
orangeCM = makeColorMap([1 1 1], [255,165,0]/255, 100)
figure('color','w','Position',[100,100, 300, 300])

ax = axes('Position', [0.24, 0.18, 0.7, 0.7])
ax.Visible='off'

cb = colorbar('Position', [0.5 0.1400 0.08 0.7500]);
cb.Ticks = '';
cb.TickLabels= '';
cb.Label.String = 'Death'
cb.Label.FontSize=17
cb.Label.Position(1) = 1
cb.Label.Color='k'
colormap(orangeCM)

%%
figname = 'oranges_colorbar'
set(gcf, 'PaperPositionMode','auto','color','w','InvertHardcopy','off')
print(gcf,'-depsc','-r600',[savpath '/epss/' figname]);
print(gcf,'-dpng','-r600',[savpath '/pngs/' figname]);

%% Colorbars
orangeCM = fire(145)
figure('color','w','Position',[100,100, 300, 300])

ax = axes('Position', [0.24, 0.18, 0.7, 0.7])
ax.Visible='off'

cb = colorbar('Position', [0.5 0.1400 0.08 0.7500]);
cb.Ticks = [0:48:145]./145;
cb.TickLabels= [0:48:145]/2;
cb.Label.String = 'Time (h)'
cb.Label.FontSize=17
cb.Label.Position(1) = 1
cb.Label.Color='w'
colormap(orangeCM)

%%
figname = 'fire_time_colorbar'
set(gcf, 'PaperPositionMode','auto','color','w','InvertHardcopy','off')
print(gcf,'-depsc','-r600',[savpath '/epss/' figname]);
print(gcf,'-dpng','-r600',[savpath '/pngs/' figname]);

%%
close all
figure('color','w','Position',[100,100, 550, 450])
ax = axes('Position', [0,0,450/550,1])

a = zeros(512);
d = strel('disk',15,0);
d = d.Neighborhood;

a(236:235+size(d,1), 251:250+size(d,2))=d;
a(242:241+size(d,1), 226:225+size(d,2))=d;
a = imclose(a,d);
As = zeros(512,512,15);
As(:,:,1) = a;
for i=2:10;
    d = strel('disk',5+5*i,0);
    As(:,:,i) = imdilate(As(:,:,i-1),d);
end

im = imagesc(sum(As,3)); shg
cmap = fire(30);

colormap([flipud(cmap(2:end-3,:))]);
axis equal
freezeColors()

orangeCM = fire(145)
ax.Color = 'none'
ax.XColor='none'
ax.YColor = 'none'

ax = axes('Position', [0.24, 0.18, 0.7, 0.7])
ax.Visible='off'

cb = colorbar('Position', [0.85 0.1400 0.08 0.7500]);
cb.Ticks = [0:48:145]./145;
cb.TickLabels= [0:48:145]/2;
cb.Color = 'k'
cb.Label.String = 'Time of first infection(h)'
cb.Label.FontSize=17
cb.Label.Position(1) = 0.25
cb.Label.Position(2) = 0.45
cb.Label.Color='w'
colormap(orangeCM)

%%
figname = 'fire_SpreadCartoon_v2'
set(gcf, 'PaperPositionMode','auto','color','w','InvertHardcopy','off')
print(gcf,'-depsc','-r600',[savpath '/epss/' figname]);
print(gcf,'-dpng','-r600',[savpath '/pngs/' figname]);








%% Fret stuff
fpath='/bigstore/Images2020/Jen/NFkBDynamics/Virus_Casp8Biosensor_2020Mar03/acq_1'
R = MultiPositionSingleCellVirusResults(fpath);
%%
dt=7

%% Treated dead
% %TNF:
ptrn1 = ['0[2-5]']
ix1 = regexp(R.PosNames, ptrn1);
ix=~cellfun('isempty',ix1);
posListDo = R.PosNames(ix);
Tracks = R.getTracks(posListDo)



inxToUse = cat(2,Tracks.CellDies);
Tracks = Tracks(find(inxToUse'));

YellowThresh = exp(meanThresh(log(cat(2,Tracks.YellowTrack))));
inxToUse = arrayfun(@(x) (nanmean(x.YellowTrack)./(0.5*YellowThresh))>1,Tracks);
Tracks = Tracks(find(inxToUse'));
%inxToUse = arrayfun(@(x) sum(((x.RatioTrack)-min(x.RatioTrack))./(nanmean(x.RatioTrack))>1),Tracks)>2;
%Tracks = Tracks(find(inxToUse'));

%% plot hist

figure('color','w','Position',[100,100, 300, 300])
ax = axes('Position', [0.2, 0.18, 0.7, 0.7])
a = (arrayfun(@(x) x.indWhenCellDies-x.indWhenC8Triggers,Tracks));
h = histogram(a,[0:2:30]);
h.FaceColor = [0 0 0]
medA = median(a(a>=0))/3;

xlabel('Time from Cas8 to death(h)','interpreter', 'tex','fontsize', 14)
ylabel('# cells','interpreter', 'tex','fontsize', 14)
set(gca, 'xtick', [0:6:40], 'XTickLabel',[0:6:40]/3)
text(5, 90, ['Median = ' num2str(medA,2) 'h'] )

%% Save
figname = 'Fig_3SFretdTStatistics'
set(gcf, 'PaperPositionMode','auto','color','w','InvertHardcopy','off','renderer','painters')
print(gcf,'-depsc','-r600',[savpath '/epss/' figname]);
print(gcf,'-dpng','-r600',[savpath '/pngs/' figname]);



%% Plot single
close all
figure(1)
set(gcf,'Position',[151 288 600 300], 'color', 'w')
for trackNum =14%:numel(Tracks);
    %cla
    
    
    YY = [InterpNANs(Tracks(trackNum).RatioTrack(1:end))]./100;%./Tracks(trackNum).NucTrack(1:end))];
    XX = [ones(1,numel(Tracks(trackNum).T))', Tracks(trackNum).T'/3];
    
    
    mSlope = [];
    linFitParams = [];
    rng = floor(dt/2)+1:length(YY)-floor(dt/2);
    for i=1:numel(rng)
        t0=rng(i);
        y= YY(t0-floor(dt/2):t0+floor(dt/2));
        x= XX(t0-floor(dt/2):t0+floor(dt/2),:);
        theta = (x'*x)\x'*y;%linear fit
        linFitParams = [linFitParams, theta] ;
        mSlope = [mSlope theta(2)];
    end
    Tracks(trackNum).slope = mSlope;
    [~,ind] = max(conv(Tracks(trackNum).slope,fliplr([-1 -1 -1 1 1 1]),'same'));
    Tracks(trackNum).indWhenC8Triggers = ind;
    ax1 = subplot(1,2,2)
    plot(XX(rng,2),mSlope,'-k')
    hold all
    h = scatter(XX(rng,2), mSlope,[],exp([Tracks(trackNum).DeathTrack(rng)]'),'filled');
    set(h,'MarkerFaceAlpha', 0.5);
    colormap(makeColorMap([0,0,1],[1,0,0]));
    
    %shg
    tzeva = viridis(numel(mSlope));
    %hold off;
    
    %hold off
    ax2 = subplot(1,2,1)
    
    for i=1:numel(rng)
        t0=rng(i);
        y= YY(t0-floor(dt/2):t0+floor(dt/2));
        x= XX(t0-floor(dt/2):t0+floor(dt/2),2);
        p=linFitParams(:,i);
        
        h = plot(x,polyval(flipud(p),x),'-','markerEdgeColor', 'none','markerFaceColor', tzeva(i,:)','Color', tzeva(i,:)');
        h(1).Color(4)=0.5;
        hold all;
        
    end
    
    h = scatter(XX(1:rng(1)-1,2), YY(1:rng(1)-1),[],tzeva(1,:),'filled');
    set(h,'MarkerFaceAlpha', 0.5)
    hold on;
    h = scatter(XX(rng(end)+1:end,2), YY(rng(end)+1:end),[],tzeva(end,:),'filled');
    set(h,'MarkerFaceAlpha', 0.5)
    hold on;
    h = scatter(XX(rng,2), YY(rng),[],tzeva,'filled');
    set(h,'MarkerFaceAlpha', 0.5)
    
    xlabel('time(h)','interpreter', 'latex','fontsize', 14)
    ylabel('Ratio $R_{Cas8}=\frac{Cyan}{Yellow}$','interpreter', 'latex','fontsize', 14)
    %set(gca,'ylim', [0.1 3])
    set(ax1,'xlim',get(gca,'xlim'));
    
    %hold off
    %pause;%(0.3)
    
end


subplot(1,2,1)
xlabel('time(h)','interpreter', 'latex')
ylabel('Ratio $R_{Cas8}=\frac{Cyan}{Yellow}$','interpreter', 'latex')
set(gca,'ylim', [0.1 150]./100,'xlim', [0, 20])
set(ax1,'xlim',get(gca,'xlim'));

subplot(1,2,2)

c = colorbar();
c.Position = [0.9 0.18 0.0254 0.65];
c.Ticks = [1.02 2.4917];
c.TickLabels={'Alive' 'Dead'};
c.FontSize=14

yl = ylabel(c,'DAPI','interpreter', 'latex','fontsize', 10);
yl.Position = yl.Position-[0.3 0 0];
xlabel('time(h)','interpreter', 'latex')
ylabel('Slope $\frac{dR_{Cas8}}{dt}$','interpreter', 'latex')
%set(gca,'ylim', [-0.1, 1],'xlim', [0, 15], 'clim', exp([0,0.08]))
    ax1.Position  = [0.56    0.18    0.3    0.6]
    ax2.Position  = [0.13    0.18    0.3    0.6]
%% Save
figname = 'Fig_3SFretSingleTrack'
set(gcf, 'PaperPositionMode','auto','color','w','InvertHardcopy','off','renderer','painters')
print(gcf,'-depsc','-r600',[savpath '/epss/' figname]);
print(gcf,'-dpng','-r600',[savpath '/pngs/' figname]);

