function welllbl = WellsConstructor_v4(fpath, Well,frame,NucChannel,varargin)
% FF is a structure array - FF.channel, FF.img

    %init WellsLbl object
    welllbl = WellsLbl;
    welllbl.PosName = Well;
    welllbl.Frame = frame;
    welllbl.pth = fpath;
    
    periRingFlag = ParseInputs('PeriRing', false, varargin);
    perisize = ParseInputs('perisize', 5, varargin);

    
    %get the data for this well in all channels
    i=frame;
    
        MD=Metadata(fpath,[],1);

        if isempty(MD.Values)
        MD=Metadata(fpath);
        end
        

        Channels = MD.unique('Channel');
        indChNuc = find(strcmp(Channels,NucChannel));

    Data = struct;
    for i=1:numel(Channels)
        img = stkread(MD,'Channel',Channels(i), 'flatfieldcorrection', false,'blindflatfield',false, 'frame', frame, 'Position', Well,'register',false);
        Data(i).channel = Channels(i);
        %FFToUse =  FF(find(strcmp(arrayfun(@(x) x.channel, FF,'uniformoutput',false),Channels(i)))).img;
        %Data(i).img = img - FFToUse +max(FFToUse(:));
        %img = img - FFToUse;
        %img(img<0)=0;
        %Data(i).img = img;
        
        if i==indChNuc
            FFimg = squeeze(awt2Dlite(img,8));
            img = sum(FFimg(:,:,2:end-1),3)+max(max(FFimg(:,:,end)));
            Data(i).img = img;
        else
            Data(i).img = backgroundSubtraction(img);
        end
        
    end
    
    
    
    %size of image
    welllbl.ImageDims = size(img);
    
    %get the nuclear channel and segment
    
    NucData = Data(indChNuc).img;
    %SizeEst =     EstSizeImg(NucData);
    [L, voronoiCells] = SegmentCellsImg_v2(NucData,varargin{:});
        
    %init ss measurements


    
    
    Ints = cell(numel(Channels),1);
    Ints90P = cell(numel(Channels),1);
    
    A = regionprops(L,'Area');
    
    CCVoronoi = bwconncomp(~voronoiCells.RegionBounds); 
      
    S = regionprops(CCVoronoi,NucData.*(L>0),'MeanIntensity','WeightedCentroid','Area','PixelValues');
    % filter CCNuclei keep only cells with volume in some range
    Areas = cat(1, S.Area);
    Centroids = cat(1,S.WeightedCentroid);
    nzAreas = arrayfun(@(x) nnz(x.PixelValues),S);
    
    NuclearIntensities = cat(1, S.MeanIntensity).*Areas;
    Nuclei90Prctile = arrayfun(@(x) prctile(x.PixelValues(x.PixelValues~=0),99.9),S);
    
    %Add to Int cell array
    Ints{indChNuc}=NuclearIntensities;
    Ints90P{indChNuc}= Nuclei90Prctile;
    
    
    
    indOtherChannels = find(~strcmp(Channels,NucChannel));
    for i=indOtherChannels'
        DataOther = Data(i).img;
        if periRingFlag
            S = regionprops(CCVoronoi,DataOther.*(L>0),'MeanIntensity','PixelValues');
        else
            S = regionprops(CCVoronoi,DataOther,'MeanIntensity','PixelValues');
        end
        Intensities = cat(1, S.MeanIntensity).*Areas;
        %prctile(S(454).PixelValues,95);
        Int90Prctile = arrayfun(@(x) prctile(x.PixelValues,99.9),S);
        Ints{i}=Intensities;
        Ints90P{i}= Int90Prctile;
    end

    if periRingFlag
        se = strel('disk',perisize,0);
        PeriL = imdilate(L,se)-L;
        IntsPeri = cell(numel(Channels),1);
        Ints90PPeri = cell(numel(Channels),1);
        for i=indOtherChannels'
            DataOther = Data(i).img;
            S = regionprops(CCVoronoi,DataOther.*(PeriL>0),'MeanIntensity','PixelValues');
            Intensities = cat(1, S.MeanIntensity).*Areas;
            %prctile(S(454).PixelValues,95);
            Int90Prctile = arrayfun(@(x) prctile(x.PixelValues,99.9),S);
            IntsPeri{i}=Intensities;
            Ints90PPeri{i}= Int90Prctile;
        end
    end
    
    
    %Apply drift correction to centroids!
    Tforms = MD.getSpecificMetadata('driftTform','Position',Well, 'frame', frame);
    if ~isempty(Tforms)
        dY = Tforms{1}(7);
        dX = Tforms{1}(8);
        n = size(Centroids,1);
        Centroids = Centroids+repmat([dY, dX],n,1);
    else
        warning('No drift correction found')
    end
    welllbl.Centroids = Centroids;
    welllbl.Areas = Areas;
    welllbl.nzAreas = nzAreas;
    if periRingFlag
        welllbl.Intensities = [Ints; IntsPeri(indOtherChannels)];
        welllbl.Int90Prctile = [Ints90P; Ints90PPeri(indOtherChannels)];
        welllbl.channels = [Channels; cellfun(@(x) [x '_peri'],Channels(indOtherChannels),'UniformOutput',false)];
    else
        welllbl.Intensities = Ints;
        welllbl.Int90Prctile = Ints90P;
        welllbl.channels = Channels;
    end
    
    welllbl.num = numel(NuclearIntensities);

    frame
    
    
    
    
end
