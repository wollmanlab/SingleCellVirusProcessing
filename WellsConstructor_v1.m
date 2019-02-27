function welllbl = WellsConstructor_v1(fpath, Well,frame,NucChannel)
% FF is a structure array - FF.channel, FF.img

    %init WellsLbl object
    welllbl = WellsLbl;
    welllbl.PosName = Well;
    welllbl.Frame = frame;
    welllbl.pth = fpath;
    
    %get the data for this well in all channels
    i=frame;
    
    MD=Metadata(fpath,[],1);
    
    if isempty(MD.Values)
        MD=Metadata(fpath);
    end
    
    Channels = MD.unique('Channel');
    Data = struct;
    for i=1:numel(Channels)
    img = stkread(MD,'Channel',Channels(i), 'flatfieldcorrection', false,'blindflatfield',false, 'frame', frame, 'Position', Well,'register',false);
    Data(i).channel = Channels(i);
    %img = img - FFToUse;
    %img(img<0)=0;
    %Data(i).img = img;
    indChNuc = find(strcmp(Channels,NucChannel));
   
    if i==indChNuc
        FFimg = squeeze(awt2Dlite(img,7));
        img = sum(FFimg(:,:,1:end-1),3);
        Data(i).img = img + mean(mean(FFimg(:,:,end)));
    else
        Data(i).img = backgroundSubtraction(img);
    end
    
    end
    
    
    
    %size of image
    welllbl.ImageDims = size(img);
    
    %get the nuclear channel and segment
    
    NucData = Data(indChNuc).img;
    %SizeEst =     EstSizeImg(NucData);
    [L, voronoiCells] = SegmentCellsThreshImg(NucData, 'EstSize',8,'hthresh',0.1,'gmthresh',0.2);
        
    %init ss measurements


    
    
    Ints = cell(numel(Channels),1);
    Ints90P = cell(numel(Channels),1);
    
    A = regionprops(L,'Area');
    
    CCVoronoi = bwconncomp(~voronoiCells.RegionBounds); 
      
    S = regionprops(CCVoronoi,NucData.*(L>0),'MeanIntensity','Centroid','Area','PixelValues');
    % filter CCNuclei keep only cells with volume in some range
    Areas = cat(1, S.Area);
    Centroids = cat(1,S.Centroid);
    nzAreas = arrayfun(@(x) nnz(x.PixelValues),S);
    
    NuclearIntensities = cat(1, S.MeanIntensity);
    Nuclei90Prctile = arrayfun(@(x) prctile(x.PixelValues(x.PixelValues~=0),90),S);
    
    %Add to Int cell array
    Ints{indChNuc}=NuclearIntensities;
    Ints90P{indChNuc}= Nuclei90Prctile;
    
    
    
    indOtherChannels = find(~strcmp(Channels,NucChannel));
    for i=indOtherChannels'
        DataOther = Data(i).img;
        S = regionprops(CCVoronoi,DataOther,'MeanIntensity','PixelValues');
        Intensities = cat(1, S.MeanIntensity);
        %prctile(S(454).PixelValues,95);
        Int90Prctile = arrayfun(@(x) prctile(x.PixelValues,90),S);
        Ints{i}=Intensities;
        Ints90P{i}= Int90Prctile;
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
    welllbl.Intensities = Ints;
    welllbl.Int90Prctile = Ints90P;
    
    welllbl.channels = Channels;
    
    welllbl.num = numel(NuclearIntensities);


    frame
    
    
    
    
end
