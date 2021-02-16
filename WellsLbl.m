classdef WellsLbl < handle %class of single cell processing of a well at a single timepoint
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
        
        
        
        
        
        function h = scatter(W,varargin)
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
                
                climits = [prctile(Ints,1), prctile(Ints,99)];
                set(gca,'xlim',[-200 3000],'ylim',[-200 2500],'clim',climits,'color','w','ydir', 'reverse')
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
                set(gca,'xlim',[-200 3000],'ylim',[-200 2500],'clim',[0 3],'color','w','ydir', 'reverse')
                colormap('viridis')
                shg
            end
            
        end
        
        
        
        
        
        function stkshow(W,varargin)
            MD=Metadata(W.pth);
            channelToShow = ParseInputs('channel', 'DeepBlue', varargin);
            Data =  stkread(MD,'Channel',channelToShow, 'flatfieldcorrection', false, 'frame', W.Frame, 'Position', W.PosName,'register',false);
            stkshow(Data);
        end
        
        function D = Density(W)
            
            [~,D] = knnsearch(W.Centroids,W.Centroids,'K', 11);
            D = 1./mean(D(:,2:end),2);
            %h = scatter(W.Centroids(:,1),W.Centroids(:,2),[],D,'filled');
            %climits = [prctile(D,1), prctile(D,99)];
            %climits = [0.0042    0.0113];
            %set(gca,'xlim',[-200 3000],'ylim',[-200 2500],'clim',climits,'color','w','ydir', 'reverse')
            %             xx = W.Centroids(:,1);
            %             yy = W.Centroids(:,2);
            %             if nargin==2
            %                 J1=1:size(xx,1);
            %             end
            %             [DensityMatrix, Bins] = hist3([xx yy], [Nbins Nbins]);
            %             %BinSize = diff(Bins{1}(1:2))*diff(Bins{2}(1:2))*(PixelSize^2);%microns^2
            %             %DensityMatrix = DensityMatrix/BinSize; %cells per micron^2 in xy
            %             DensityMatrix = DensityMatrix./mean(DensityMatrix(:));
            %             DensityMatrix = imgaussfilt(DensityMatrix,1);
            %             Density = interp2(Bins{2}, Bins{1}, DensityMatrix, yy(J1), xx(J1), 'spline');
        end
        
        
        function scattershow(W,varargin)
            try
                MIJ.run('Close All');
            catch
            end
            MD=Metadata(W.pth,[],1);
            channelToShow = ParseInputs('channel', 'DeepBlue', varargin);
            showData = ParseInputs('showData', true, varargin);
            
            if ~isempty(regexpi(channelToShow,'_peri'))
                channelToShowimg = channelToShow(1:regexpi(channelToShow,'_peri')-1);
            else
                channelToShowimg = channelToShow;
            end
            Data =  stkread(MD,'Channel',channelToShowimg, 'flatfieldcorrection', false, 'frame', W.Frame, 'Position', W.PosName,'register',false);
            RChannel=zeros(size(Data));
            GChannel=zeros(size(Data));
            BChannel=zeros(size(Data));
            [h,x] = hist(log(datasample(Data(:),min(1000000,numel(Data(:))))),1000);
            maxC = exp(x(find(cumsum(h)./sum(h)>0.995,1,'first')));
            
            range = 1:W.num;
            
            indChNuc = find(strcmp(W.channels,channelToShow));
            %tzeva = viridis(length(W.Centroids));
            Ints = W.Int90Prctile{indChNuc};
            %Ctoplot = (Ints-0.9*min(Ints))./(max(Ints)-0.9*min(Ints));

            climits = [prctile(Ints,1), prctile(Ints,99)];
            Ints = (Ints-climits(1))/(climits(2)-climits(1));
            
            cmap = viridis(1024)*maxC*1.5;%scale colormap so that data and centroids are all visible
            
            CtoShow = interp1(linspace(0,1,1024), 1:1024,Ints,'nearest','extrap');
            Tforms = MD.getSpecificMetadata('driftTform','Position',W.PosName, 'frame', W.Frame);
            for i=1:numel(range)
                Cent = W.Centroids(range(i),:);
                if ~isnan(Cent)
                    
                    if ~isempty(Tforms)
                        dY = Tforms{1}(7);
                        dX = Tforms{1}(8);
                        n = size(Cent,1);
                        Cent = Cent-repmat([dY, dX],n,1);
                        
                    else
                        warning('No drift correction found')
                    end
                    Cent = round(Cent);
                    indexToChange = sub2ind(size(Data),min(max(Cent(:,2),1),size(Data,1)),min(max(Cent(:,1),1),size(Data,2)));
                    
                    RChannel(indexToChange)=cmap(CtoShow(i),1); %Replace actual pixel values with relevant RGB value
                    GChannel(indexToChange)=cmap(CtoShow(i),2);
                    BChannel(indexToChange)=cmap(CtoShow(i),3);
                    i;
                end
            end
            se = strel('disk',7);
            RGB = cat(4,imdilate(RChannel,se),imdilate(GChannel,se),imdilate(BChannel,se));%combine into a single RGB image and show. stkshow is sloooooooooooowing me down.
            stkshow(RGB);
            MIJ.selectWindow('RGB');
            %MIJ.run('Stack to Hyperstack...', ['order=xyzct channels=3 slices=' num2str(size(Data,3)) ' frames=1 display=Composite']);
            MIJ.run('Stack to RGB');
            if showData
                stkshow(Data)
                MIJ.run('Add Image...', 'image=[RGB (RGB)] x=0 y=0 opacity=100 zero');
            end
            ij.IJ.selectWindow('RGB');
            ij.IJ.run('Close');
            ij.IJ.selectWindow('RGB (RGB)');
            ij.IJ.run('Close');

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
        
        
        
        function r = Ints(W,Ch1, range)
            
            indCh1 = find(strcmp(W.channels,Ch1));
            if nargin==2
                range = 1:W.num;
            elseif any(range>W.num)
                range = 1:W.num;
            end
            range(range==0)=[];
            
            r = W.Intensities{indCh1}(range);
        end
        
        
        function r = Ints90(W,Ch1, range)
            indCh1 = find(strcmp(W.channels,Ch1));
            if nargin==2
                range = 1:W.num;
            elseif any(range>W.num)
                range = 1:W.num;
            end
            range(range==0)=[];

            r = W.Int90Prctile{indCh1}(range);
        end
        
        
        
    end
end
