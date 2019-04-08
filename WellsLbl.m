classdef WellsLbl < handle %class of single cell processing of a whole cornea at a single timepoint
    properties
        PosName
        acq
        Frame
        pth
        ImageDims
        channels
        
        Centroids
        Intensities
        Int90Prctile       
        nzAreas
        Areas
        num
        
        %CCNuclei
        %CCVoronoi
        %Properties related to tracking
        %Link12Mat
        Link21Mat
    end
    
    properties (Transient = true)
        verbose = true;
    end
    
    properties (Dependent = true)
        filenames
    end
    
    
    
    methods
        
        %         function epiRank = epiRank(W)
        %             [h,x] = hist(W.epiScore,250);
        %             dataCDF = cumsum(h)/sum(h);
        %             epiRank = interp1(x,dataCDF,W.epiScore);
        %             epiRank(isnan(W.epiScore))=NaN;
        %         end
        

        
        
        
        function scatter(W,varargin)
            channelToShow = ParseInputs('channel', W.channels{1}, varargin);
            ratio = ParseInputs('ratio', false, varargin);
            rangeToPlot = ParseInputs('range', 1:W.num, varargin);
            if ~ratio
            indChNuc = find(strcmp(W.channels,channelToShow));
                %tzeva = viridis(length(W.Centroids));
                Ints = W.Int90Prctile{indChNuc};
            %Ctoplot = (Ints-0.9*min(Ints))./(max(Ints)-0.9*min(Ints));
            Ctoplot = Ints;
            h = scatter(W.Centroids(rangeToPlot,1),W.Centroids(rangeToPlot,2),[],Ctoplot(rangeToPlot),'filled');
            set(gca,'xlim',[-200 3000],'ylim',[-200 2500],'clim',[0 1],'color','k','ydir', 'reverse')
            colormap('viridis')
            shg
            shg
            else
            channel1 = ParseInputs('channel1', 'Cyan', varargin);
            channel2 = ParseInputs('channel2', 'Yellow', varargin);
            indCh1 = find(strcmp(W.channels,channel1));
            indCh2 = find(strcmp(W.channels,channel2));

                Ints = W.Intensities{indCh1}./W.Intensities{indCh2};
            %Ctoplot = (Ints-prctile(Ints,5))./(prctile(Ints,100)-prctile(Ints,5));
            h = scatter(W.Centroids(rangeToPlot,1),W.Centroids(rangeToPlot,2),[],Ints(rangeToPlot),'filled');
            set(gca,'xlim',[-200 3000],'ylim',[-200 2500],'clim',[0 3],'color','k','ydir', 'reverse')
            colormap('viridis')
            shg
            end

        end
        
        
        
        
        
        function stkshow(W)
            MD=Metadata(W.pth);
            
            Data =  stkread(MD,'Channel','DeepBlue', 'flatfieldcorrection', false, 'frame', W.Frame, 'Position', W.PosName,'register',false);
            RChannel=Data;
            GChannel=Data;
            BChannel=Data;
            [h,x] = hist(log(datasample(Data(:),min(1000000,numel(Data(:))))),1000);
            maxC = exp(x(find(cumsum(h)./sum(h)>0.995,1,'first')));
            cmap = parula(W.CC.NumObjects)*maxC*1.5;%scale colormap so that data and centroids are all visible
            
            for i=1:W.CC.NumObjects
                indexToChange = W.CC.PixelIdxList{i};
                RChannel(indexToChange)=cmap(i,1); %Replace actual pixel values with relevant RGB value
                GChannel(indexToChange)=cmap(i,2);
                BChannel(indexToChange)=cmap(i,3);
                i
            end
            RGB = cat(3,RChannel,GChannel,BChannel);%combine into a single RGB image and show. stkshow is sloooooooooooowing me down.
            stkshow(RGB);
            MIJ.selectWindow('RGB');
            MIJ.run('Stack to Hyperstack...', ['order=xyzct channels=3 slices=' num2str(size(Data,3)) ' frames=1 display=Composite']);
            %stkshow(RChannel)
            %stkshow(GChannel)
            %stkshow(BChannel)
            %MIJ.run('Merge Channels...', 'c1=RChannel c2=GChannel c3=BChannel create');
        end
        
        function Density = Density(W,Nbins,J1)
%            MD=Metadata(W.pth);
%            PixelSize = MD.unique('PixelSize');
            xx = W.Centroids(:,1);
            yy = W.Centroids(:,2);
            if nargin==2
                J1=1:size(xx,1);
            end
            J = W.Jepi;
            [DensityMatrix, Bins] = hist3([xx(J) yy(J)], [Nbins Nbins]);
            %BinSize = diff(Bins{1}(1:2))*diff(Bins{2}(1:2))*(PixelSize^2);%microns^2
            %DensityMatrix = DensityMatrix/BinSize; %cells per micron^2 in xy
            DensityMatrix = DensityMatrix./mean(DensityMatrix(:));
            DensityMatrix = imgaussfilt(DensityMatrix,1);
            Density = interp2(Bins{2}, Bins{1}, DensityMatrix, yy(J1), xx(J1), 'spline');
        end
        
        
        function scattershow(W,varargin)
            MD=Metadata(W.pth);
            channelToShow = ParseInputs('channel', 'DeepBlue', varargin);
            Data =  stkread(MD,'Channel',channelToShow, 'flatfieldcorrection', false, 'frame', W.Frame, 'Position', W.PosName,'register',false);
            RChannel=zeros(size(Data));
            GChannel=zeros(size(Data));
            BChannel=zeros(size(Data));
            [h,x] = hist(log(datasample(Data(:),min(1000000,numel(Data(:))))),1000);
            maxC = exp(x(find(cumsum(h)./sum(h)>0.995,1,'first')));

                range = 1:W.num;
            
            cmap = viridis(numel(range))*maxC*1.5;%scale colormap so that data and centroids are all visible

            for i=1:numel(range)
                indexToChange = W.CC.PixelIdxList{range(i)};
                RChannel(indexToChange)=cmap(i,1); %Replace actual pixel values with relevant RGB value
                GChannel(indexToChange)=cmap(i,2);
                BChannel(indexToChange)=cmap(i,3);
                i
            end
            RGB = cat(3,RChannel,GChannel,BChannel);%combine into a single RGB image and show. stkshow is sloooooooooooowing me down.
            stkshow(RGB);
            MIJ.selectWindow('RGB');
            MIJ.run('Stack to Hyperstack...', ['order=xyzct channels=3 slices=' num2str(size(Data,3)) ' frames=1 display=Composite']);
        end
        
        function r = ratioChannels(W,Ch1, Ch2, range)
            indCh1 = find(strcmp(W.channels,Ch1));
            indCh2 = find(strcmp(W.channels,Ch2));
            if nargin==3
                range = 1:W.num;
            elseif any(range>W.num)
                range = 1:W.num;
            end
            r = W.Intensities{indCh1}(range)./W.Intensities{indCh2}(range);
        end
        
    end
end
