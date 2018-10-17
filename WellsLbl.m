classdef WellsLbl < handle %class of single cell processing of a whole cornea at a single timepoint
    properties
        PosName
        Frame
        pth
        ImageDims
        
        
        Centroids
        Intensities
        Int90Prctile
        channels
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
            channelToShow = ParseInputs('channel', 'virus', varargin);
            if strcmp(channelToShow,'virus')
                %tzeva = viridis(length(W.Centroids));
                Ctoplot = (W.Virus90Prctile-0.08)./(0.13-0.08);
                h = scatter(W.Centroids(:,1),W.Centroids(:,2),[],Ctoplot,'filled');
                set(gca,'xlim',[-200 3000],'ylim',[-200 2500],'clim',[0 1],'color','k')
                colormap('viridis')
            elseif strcmp(channelToShow,'nuclei')
                %tzeva = plasma(length(W.Centroids));
                Ctoplot = (W.Nuclei90Prctile-0.15)./(0.8-0.15);

                h = scatter(W.Centroids(:,1),W.Centroids(:,2),[],Ctoplot,'linewidth',2);
                set(gca,'xlim',[-200 3000],'ylim',[-200 2500],'clim',[0,1],'color','k')
                colormap(makeColorMap([0 0 0], [1 0 0]))
                
                %             else
                %                 range = varargin{1};
                %
                %                 if any(range>W.num)
                %                     range = range(1):W.num;
                %                     disp('range out of bounds, drawing all points up to the total # of cells.')
                %                     scatter(W.Centroids(range,1),W.Centroids(range,2),[],tzeva(range,:));
                %                 else
                %                     scatter(W.Centroids(range,1),W.Centroids(range,2),[],tzeva(range,:));
                %                 end;
                
            end
            %scatter(Centroids(:,1),Centroids(:,2),[],parula(length(Centroids)));
            %hold on
            shg
            
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
            
            Data =  stkread(MD,'Channel','DeepBlue', 'flatfieldcorrection', false, 'frame', W.Frame, 'Position', W.PosName,'register',false);
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
        
    end
end
