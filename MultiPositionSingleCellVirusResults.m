classdef MultiPositionSingleCellVirusResults < MultiPositionResults
    properties
        Frames = {};
        WellLbls
        Tracks
        figSavePath
    end
    
    properties (Dependent = true)
        TimeVecs
    end
    
    methods (Static)
        function R = loadobj(S)
            R = MultiPositionSingleCellVirusResults;
            R = reload(R,S);
        end
    end
    
    
    methods
        
        function S = saveobj(R)
            S = toStruct(R);
        end
        
        function R = reload(R,S)
            R = reload@MultiPositionResults(R,S);
            %             if isfield(S,'PIVlbl') % is the fields exists load them, if not,
            %                 % they will be filled with default values
            %                 % effectivly upcasting an object.
            %                 %R.WoundLbl = S.WoundLbl;
            %                 R.PIVlbl = S.PIVlbl;
            %             end
            if isfield(S,'WellLbls') % is the fields exists load them, if not,
                % they will be filled with default values
                % effectivly upcasting an object.
                %R.WoundLbl = S.WoundLbl;
                R.WellLbls = S.WellLbls;
            end
            if isfield(S,'Frames') % is the fields exists load them, if not,
                % they will be filled with default values
                % effectivly upcasting an object.
                R.Frames = S.Frames;
            end
            if isfield(S,'Tracks') % is the fields exists load them, if not,
                % they will be filled with default values
                % effectivly upcasting an object.
                R.Tracks = S.Tracks;
            end
        end
        
        function S = toStruct(R)
            % call the superclass method to start the transition to a
            % struct
            S = toStruct@MultiPositionResults(R);
            
            %S.WoundLbl = R.WoundLbl;
            %S.PIVlbl = R.PIVlbl;
            S.WellLbls = R.WellLbls;
            S.Tracks = R.Tracks;
            S.Frames = R.Frames;
            %             S.TimeVecs = R.TimeVecs;
        end
        
        function R = MultiPositionSingleCellVirusResults(pth,reset) %constructor
            if nargin==0
                pth='';
                reset=false;
            end
            if nargin==1
                reset=false;
            end
            R@MultiPositionResults(pth,reset);
        end
        
        %function TimeVecs = get.TimeVecs(R)
        %    TimeVecs = getTimeVecs(R);
        %end
        
        function TimeVecs = get.TimeVecs(R)
            for i=1:R.Np
                pos = R.PosNames{i};
                ix = ismember(R.PosNames,pos);
                if ~isempty(R.WellLbls{ix})
                    TimeVecs{ix}=cell2mat(R.WellLbls{ix}.Tvec);
                else
                    TimeVecs{ix} = [];
                end
            end;
        end
        
        function T = getTimeVecByPosition(R,pos)
            assertPositionExist(R,pos);
            ix = ismember(R.PosNames,pos);
            T =  R.TimeVecs{ix};
        end
        
        
        
        
        
        
        
        function setTracks(R,trcks,pos,varargin)
            if numel(R.Tracks)<R.Np
                R.Tracks{R.Np}=[];
            end
            
            arg.redo = false;
            arg = parseVarargin(varargin,arg);
            
            % do some checks
            %assert(isa(trcks,'struct'),'second variable trcks must be of Class struct');
            assertPositionExist(R,pos);
            if iscell(pos)
                if numel(pos)>1
                    error('One position at a time, please.')
                end
            else
                pos={pos};
            end
            psname = R.PosNames;
            ix = ismember(psname,pos);
            
            R.Tracks{ix}=trcks;
        end
        
        function PL = getTracks(R,pos)
            assertPositionExist(R,pos);
            ix = ismember(R.PosNames,pos);
            PL = cat(1,R.Tracks{ix});
        end
        
        
        
        
        
        function setWellsLbl(R,CLbl,pos,varargin)
            if numel(R.WellLbls)<R.Np
                R.WellLbls{R.Np}=[];
            end
            
            arg.redo = false;
            arg = parseVarargin(varargin,arg);
            
            % do some checks
            assert(isa(CLbl{1},'WellsLbl'),'second variable PLbl must be of Class WellsLbl');
            assertPositionExist(R,pos);
            if iscell(pos)
                if numel(pos)>1
                    error('One position at a time, please.')
                end
            else
                pos={pos};
            end
            psname = R.PosNames;
            ix = ismember(psname,pos);
            
            R.WellLbls{ix}=CLbl;
        end
        
        function PL = getWellsLbl(R,pos)
            assertPositionExist(R,pos);
            ix = ismember(R.PosNames,pos);
            PL = R.WellLbls{ix};
        end
        
        
        function W = mergeWellLbls(R,Wvec)
            W = WellsLbl;
            %first, stuff that are constant for an experiment
            if numel(unique(arrayfun(@(x) x.pth,Wvec,'uniformOutput' ,0)))~=1
                error('Only merge stuff from the same experiment please!')
            end
            W.pth = Wvec(1).pth;
            W.ImageDims = Wvec(1).ImageDims;
            W.channels =  Wvec(1).channels;
            
            %next, properties of positions
            Pos = cell(numel(Wvec),1);
            pth = cell(numel(Wvec),1);
            Frame = [];
            for i=1:numel(Wvec)
                Pos{i} = Wvec(i).PosName;
                Frame = [Frame Wvec(i).Frame];
            end
            W.PosName=Pos;
            W.Frame=Frame;
            
            %Finally, data
     
            a = arrayfun(@(x) x.Centroids,Wvec,'uniformOutput' ,0);
            W.Centroids = cat(1,a{:});
            a = arrayfun(@(x) x.Intensities,Wvec,'uniformOutput' ,0);
            a = cellfun(@(x) cat(1,x), cellfun(@(x) cat(2,x{:}),a,'uniformOutput',0),'UniformOutput',0);
            W.Intensities = mat2cell(cat(1,a{:}),size(cat(1,a{:}),1),ones(numel(W.channels),1));
            a = arrayfun(@(x) x.Int90Prctile,Wvec,'uniformOutput' ,0);
            a = cellfun(@(x) cat(1,x), cellfun(@(x) cat(2,x{:}),a,'uniformOutput',0),'UniformOutput',0);
            W.Intensities = mat2cell(cat(1,a{:}),size(cat(1,a{:}),1),ones(numel(W.channels),1));
           
            a = arrayfun(@(x) x.nzAreas,Wvec,'uniformOutput' ,0);
            W.nzAreas = cat(1,a{:});
            a = arrayfun(@(x) x.Areas,Wvec,'uniformOutput' ,0);
            W.Areas = cat(1,a{:});
            W.num = numel(W.Areas);
            
        end
        
        
        
        function T = mergeTracks(R, Tvec)
            T = Tvec{1};
            for i=2:numel(Tvec)
                T = cat(1,T,Tvec{i});
            end
        end
        
        
        
        
        function Link(R,pos,trckChnl, varargin)
            %% Now, we'll see how well we can track between adjacent frames (lap)
            
            %init assignment matrices
            %Link12MatCell = {};
            Link21MatCell = {};
            
            searchRadius = ParseInputs('searchRadius', 60, varargin);
            intMultiplier = ParseInputs('intMultiplier', 10^2, varargin);
            % maxAmpRatio = 1.5;
            WellCells = R.getWellsLbl(pos);
            inds = find(cellfun(@(x) ~isempty(x), WellCells));
            indtrckChnl = find(strcmp(trckChnl,WellCells{1}.channels));
            
            for i=inds(1:end-1)';
                %build cost function
                if (isempty(WellCells{i}.Centroids)) || (isempty(WellCells{i+1}.Centroids))
                    Link21MatCell{i} = [];
                else
                    Dists  = createDistanceMatrix(WellCells{i}.Centroids,WellCells{i+1}.Centroids);
                    costMat = Dists;
                    
                    for jj = indtrckChnl
                    d1 = (WellCells{i}.Int90Prctile{jj}-nanmean(WellCells{i}.Int90Prctile{jj}))./nanstd(WellCells{i}.Int90Prctile{jj});
                    a1 = repmat(d1,1,numel(WellCells{i+1}.Int90Prctile{jj}));
                    d2 = (WellCells{i+1}.Int90Prctile{jj}-nanmean(WellCells{i+1}.Int90Prctile{jj}))./nanstd(WellCells{i+1}.Int90Prctile{jj});
                    a2 = repmat(d2',numel(WellCells{i}.Int90Prctile{jj}),1);
                    VirDists = abs(a2-a1);
                    
                    costMat = costMat+intMultiplier*VirDists;%;
                    end
                    %divide the larger of the two amplitudes by the smaller value
                    %               ampRatio = a1./a2;
                    %               J = ampRatio < 1;
                    %               ampRatio(J) = 1./ampRatio(J);
                    
                    
                    
                    costMat(Dists>searchRadius) = 0;
                    %             costMat(ampRatio>maxAmpRatio) = 0;
                    costMat(isnan(costMat)) = 0;
                    
                    costMat = sparse(double(costMat));
                    %
                    [~, Links21] = lap(costMat,[],[],1, 200);
                    %Make matrix of connections found
                    Link21Mat =repmat(Links21(1:length(WellCells{i+1}.Centroids(:,1))),1,length(WellCells{i}.Centroids(:,1)));
                    
                    LinkMat2 = meshgrid(1:length(WellCells{i}.Centroids(:,1)), 1:length(WellCells{i+1}.Centroids(:,1)));
                    
                    Link21Mat = LinkMat2==Link21Mat;
                    Link21MatCell{i}=sparse(Link21Mat);
                end
            end
            
            for i=inds(1:end-1)';
                WellCells{i}.Link21Mat = Link21MatCell{i};
            end
            R.setWellsLbl(WellCells,pos)
        end
        
        function closeGaps(R,pos,trckChnl,NucChnl,varargin)
            
            
            
            maxTimeJump = ParseInputs('maxTimeJump', 4, varargin);
            maxStep = ParseInputs('maxStep', 20, varargin);
            mergeSplit = ParseInputs('mergeSplit', 0, varargin);
            maxAmpDiff = ParseInputs('maxAmpDiff', 0.05, varargin);
            whereToStart = ParseInputs('wheretostart', 1, varargin);
            whereToEnd = ParseInputs('wheretoend', numel(R.Frames), varargin);

            
            WellCells = R.getWellsLbl(pos);
            WellCells = WellCells(whereToStart:whereToEnd);

            indtrckChnl = find(strcmp(trckChnl,WellCells{1}.channels));
            indNucChnl = find(strcmp(NucChnl,WellCells{1}.channels));
            
            
            inds = find(cellfun(@(x) ~isempty(x), WellCells));
            %% Book keeping: Make all tracks fragmants
            numFeatures = cellfun(@(x) x.num, WellCells(inds))';
            trackedFeatureIndx = (1:numFeatures(1))';
            numFrames = numel(inds);
            
            %initialize auxiliary matrices for storing information related to tracks
            %fragments
            numTracksWorstCase = round(sum(numFeatures)/10); %arbitrary large number
            if numTracksWorstCase>0;
                
                trackedFeatureIndxAux = zeros(numTracksWorstCase,numFrames);
                rowEnd = numTracksWorstCase; %We'll fill this from the bottom up
                
                for i=1:numFrames-1
                    %get indices of features in 2nd frame that are connected to features in 1st frame
                    %indx1C - indexes in frame 1 that are linked to frame 2
                    %indx2C - indexes in frame 2 that are linked to indx1C in frame 1
                    
                    
                    numFeaturesFrame1 = numFeatures(i);
                    numFeaturesFrame2 = numFeatures(i+1);
                    [indx2C,indx1C] = find(WellCells{i}.Link21Mat);
                    
                    %%
                    %find existing tracks that are not connected to features in 2nd frame
                    numExistTracks = size(trackedFeatureIndx,1);
                    indx1U = setdiff(1:numExistTracks,indx1C); %features in 1 not connected to 2
                    numRows = length(indx1U);
                    %%
                    %determine where to store these tracks in auxiliary matrix
                    %extend auxiliary matrices if necessary
                    rowStart = rowEnd - numRows + 1;
                    while rowStart <= 1
                        trackedFeatureIndxAux = [zeros(numTracksWorstCase,numFrames); ...
                            trackedFeatureIndxAux];
                        rowEnd = rowEnd + numTracksWorstCase;
                        rowStart = rowStart + numTracksWorstCase;
                    end
                    
                    %% move rows of tracks that are not connected to points in
                    %2nd frame to auxilary matrix
                    trackedFeatureIndxAux(rowStart:rowEnd,1:i) = trackedFeatureIndx(indx1U,:);
                    
                    %%
                    %assign space for new connectivity matrix
                    tmp = zeros(numFeaturesFrame2,i+1);
                    %fill in the feature numbers in 2nd frame
                    tmp(1:numFeaturesFrame2,i+1) = (1:numFeaturesFrame2)';
                    %shuffle existing tracks to get the correct connectivity with 2nd frame
                    tmp(indx2C,1:i) = trackedFeatureIndx(indx1C,:);
                    %update the connectivity matrix "trackedFeatureIndx"
                    trackedFeatureIndx = tmp;
                    
                    %update rowEnd to indicate until which row the auxiliary
                    %matrices are ampty
                    rowEnd = rowStart - 1;
                end
                %
                %add information from last frame to auxiliary matrices
                numRows = size(trackedFeatureIndx,1);
                rowStart = rowEnd - numRows + 1;
                if rowStart <= 1
                    trackedFeatureIndxAux = [zeros(numRows,numFrames); ...
                        trackedFeatureIndxAux];
                    rowEnd = rowEnd + numRows;
                    rowStart = rowStart + numRows;
                end
                trackedFeatureIndxAux(rowStart:rowEnd,:) = trackedFeatureIndx;
                
                %remove all empty rows
                trackedFeatureIndx = trackedFeatureIndxAux(rowStart:end,:);
                clear trackedFeatureIndxAux
                
                % get total number of tracks
                numTracks = size(trackedFeatureIndx,1);
                
                %find the frame where each track begins and then sort the vector
                frameStart = zeros(numTracks,1);
                for i=1:numTracks
                    frameStart(i) = find((trackedFeatureIndx(i,:)~=0),1,'first');
                end
                [frameStart,indx] = sort(frameStart);
                
                %rearrange "trackedFeatureIndx" such that tracks are sorted in ascending order by their
                %starting point. Note that this ends up also arranging tracks starting at the
                %same time in descending order from longest to shortest.
                trackedFeatureIndx = trackedFeatureIndx(indx,:);
                
                
                %% Filter short tracks
                movieInfo = [];
                for i=inds'
                    n = WellCells{i}.num;
                    if n~=0
                        movieInfo(i).xCoord = [WellCells{i}.Centroids(:,1) zeros(n,1)];
                        movieInfo(i).yCoord = [WellCells{i}.Centroids(:,2) zeros(n,1)];
                        movieInfo(i).zCoord = [WellCells{i}.Int90Prctile{indtrckChnl}  zeros(n,1)];
                        movieInfo(i).amp = [WellCells{i}.Int90Prctile{indNucChnl}  zeros(n,1)];
                        movieInfo(i).num = WellCells{i}.num;
                    else
                        movieInfo(i).xCoord = [];
                        movieInfo(i).yCoord = [];
                        movieInfo(i).zCoord = [];
                        movieInfo(i).amp = [];
                        movieInfo(i).num = WellCells{i}.num;
                    end
                end
                probDim = 3;
                trackedFeatureInfo = coordAmpMatFromIndicesSparse(trackedFeatureIndx,movieInfo,...
                    numFrames,probDim);
                trackSEL = getTrackSEL(trackedFeatureInfo);
                
                
                
                
                %remove stubs whose length is less than minTrackLen
                minTrackLen=10;
                
                indxKeep = find(trackSEL(:,3) >= minTrackLen);
                trackSEL = trackSEL(indxKeep,:);
                trackedFeatureIndx = trackedFeatureIndx(indxKeep,:);
                trackedFeatureInfo = trackedFeatureInfo(indxKeep,:);
                numTracks = size(trackSEL,1)
                clear movieInfo;
                
                
                %% Close gaps with merging/splitting
                
                %if there are gaps to close (i.e. if there are tracks that start after the
                %first frame and tracks that end before the last frame) ...
                
                
                
                numTracksLink = size(trackedFeatureIndx,1);
                
                
                %% Find all possible links based on thresholds

                trackStartTime = trackSEL(:,1);
                trackEndTime   = trackSEL(:,2);
                
                CentroidsStarts = cell2mat(arrayfun(@(a) full(trackedFeatureInfo(a,8*(trackStartTime(a)-1)+1:8*(trackStartTime(a)-1)+2))',1:numTracks,'UniformOutput',false))';
                CentroidsEnds = cell2mat(arrayfun(@(a) full(trackedFeatureInfo(a,8*(trackEndTime(a)-1)+1:8*(trackEndTime(a)-1)+2))',1:numTracks,'UniformOutput',false))';
                
                
                AmpStarts = cell2mat(arrayfun(@(a) full(trackedFeatureInfo(a,8*(trackStartTime(a)-1)+4:8*(trackStartTime(a)-1)+4))',1:numTracks,'UniformOutput',false))';
                AmpEnds = cell2mat(arrayfun(@(a) full(trackedFeatureInfo(a,8*(trackEndTime(a)-1)+4:8*(trackEndTime(a)-1)+4))',1:numTracks,'UniformOutput',false))';

                
                indStarts = cell2mat(arrayfun(@(a) trackedFeatureIndx(a,trackStartTime(a)),1:numTracks,'UniformOutput',false))';
                indEnds = cell2mat(arrayfun(@(a) trackedFeatureIndx(a,trackEndTime(a)),1:numTracks,'UniformOutput',false))';
  
                
                indx1=[];
                indx2=[];
                Dists=[];
                Skips = [];
                
                f = waitbar(0,'Calculating ditances of possible links...');
                for ind1 = 1:numTracks
                    if ~rem(ind1,100)
                        waitbar(ind1/numTracks,f,'Calculating ditances of possible links...');
                    end
          
                    for ind2=1: numTracks

                               
                        dT = trackStartTime(ind2)-trackEndTime(ind1);

                        if dT>0 && dT<=maxTimeJump;%condition on time
                            dR = norm(CentroidsStarts(ind2,:)-CentroidsEnds(ind1,:));
                            if dR<=sqrt(dT)*maxStep;%condition on space
                                if abs(AmpStarts(ind2)-AmpEnds(ind1))<maxAmpDiff %condition on amp
                                    indx1=[indx1 ind1];
                                    indx2=[indx2 ind2];
                                    Dists=[Dists dR];
                                    Skips=[Skips dT];
                                end
                            end
                        end
                    end
                end
                close(f)
                %% calculate cost matrix
                f = waitbar(0,'Gap closing...');
                if nnz(trackedFeatureInfo)==0
                    longtracksFinal = [];
                else
                    [trackStats,statsRelChange,errFlag] = getTrackStats(trackedFeatureInfo,0.70,maxTimeJump);
                    dispSqTheta = trackStats.dispSqTheta;
                    dispSqR = trackStats.dispSqR;
                    addConst =   - dispSqR.*log(dispSqTheta) + log(gamma(dispSqR));%log(ampDiffStd)
                    
                    costs = [];
                    
                    for i=1 : numel(indx1)
                        if ~rem(i,100)
                            waitbar(i/numel(indx1),f,'Calculating cost matrix for gap closing...');
                        end
                        dta = Skips(i);
                        distA = Dists(i);
                        dispSqT = dispSqTheta(dta);
                        dispSqAR = dispSqR(dta);
                        addCons = addConst(dta);
                        costs = [costs, dispSqT.*distA.^2-(dispSqAR-1).*log(max(distA.^2,realmin))+addCons];
                    end
                    close(f)
                    
                    %list the tracks that start and end in each frame
                    tracksPerFrame = repmat(struct('starts',[],'ends',[]),numFrames,1);
                    for iFrame = 1 : numFrames
                        tracksPerFrame(iFrame).starts = find(trackStartTime == iFrame); %starts
                        tracksPerFrame(iFrame).ends = find(trackEndTime == iFrame); %ends
                    end
                    
                    %            %todo - add costs for cell splitting
                    numSplit  =  0; %index counting splitting events
                    indxSplit = []; %vector storing splitting track number
                    altCostSplit = []; %vector storing alternative costs of not splitting
                    maxDispAllowed = 75;
                    %
                    %go over all tracks
                    [numTracksLink,numFrames] = size(trackedFeatureIndx);
                    trackCenter = zeros(numTracksLink,probDim);
                    trackMeanDispMag = NaN(numTracksLink,1);
                    for iTrack = 1 : numTracksLink
                        
                        %get current track's coordinates
                        currentTrack = (reshape(trackedFeatureInfo(iTrack,:)',8,[]))';
                        currentTrack = currentTrack(:,1:probDim);
                        currentTrack = full(currentTrack(trackStartTime(iTrack):trackEndTime(iTrack),:));
                        
                        %calculate the track's center of mass
                        trackCenter(iTrack,:) = mean(currentTrack);
                        
                        %calculate the track's mean displacement
                        if size(currentTrack,1) > 1
                            trackMeanDispMag(iTrack) = mean(sqrt(sum(diff(currentTrack,1,1).^2,2)));
                        end
                        
                    end
                    
                    %costs of splitting
                    if mergeSplit
                        
                        %go over all track starting times
                        for startTime = 2 : numFrames
                            
                            %find tracks that start in this frame
                            startsToConsider = tracksPerFrame(startTime).starts;
                            
                            %find tracks that start before this frame and end after or in this frame
                            splitsToConsider = intersect(vertcat(tracksPerFrame(1:startTime-1).starts),...
                                vertcat(tracksPerFrame(startTime:end).ends));
                            
                            %get index indicating time of splitting
                            timeIndx  = (startTime-2)*8;
                            
                            %calculate displacement between track starts and other tracks in the
                            %previous frame
                            dispMat2 = createDistanceMatrix(CentroidsStarts(startsToConsider,:), ...
                                full(trackedFeatureInfo(splitsToConsider,timeIndx+1:timeIndx+size(CentroidsStarts,2))));
                            
                            %find possible pairs
                            [indxStart2,indxSplit2] = find(dispMat2 <= maxDispAllowed);
                            numPairs = length(indxStart2);
                            
                            %clear memory
                            clear dispMat2
                            
                            %map from indices to track indices
                            indxStart2 = startsToConsider(indxStart2);
                            indxSplit2 = splitsToConsider(indxSplit2);
                            
                            %reserve memory for cost vectors and related vectors
                            indx1MS   = zeros(numPairs,1);
                            indx2MS   = zeros(numPairs,1);
                            costMS    = zeros(numPairs,1);
                            altCostMS = zeros(numPairs,1);
                            indxMSMS  = zeros(numPairs,1);
                            
                            %go over all possible pairs
                            for iPair = 1 : numPairs
                                
                                %get indices of starting track and track it might have split from
                                iStart = indxStart2(iPair);
                                iSplit = indxSplit2(iPair);
                                
                                %calculate the vector connecting the end of track iStart to the
                                %point of splitting and compute its magnitude
                                dispVec = CentroidsStarts(iStart,:) - full(trackedFeatureInfo(iSplit,...
                                    timeIndx+1:timeIndx+size(CentroidsStarts,2)));
                                dispVecMag = sqrt(dispVec * dispVec');
                                
                                %get the amplitude of track iStart after its start - take
                                %the first nTpMS+1 points
                                nTpMS=2;
                                indxAfter = 8*(startTime-1)+3 + 8*(0:nTpMS);
                                indxAfter = indxAfter(indxAfter < 8*numFrames);
                                ampS = full(trackedFeatureInfo(iStart,indxAfter));
                                ampS = mean(ampS(ampS~=0));
                                
                                %get the amplitude of the splitting track after and before
                                %splitting - take nTpMS+1 points on each side
                                ampSp1 = full(trackedFeatureInfo(iSplit,indxAfter)); %after splitting
                                ampSp1 = mean(ampSp1(ampSp1~=0));
                                indxBefore = 8*(startTime-2)+3 - 8*(0:nTpMS);
                                indxBefore = indxBefore(indxBefore > 1);
                                ampSp = full(trackedFeatureInfo(iSplit,indxBefore)); %before splitting
                                ampSp = mean(ampSp(ampSp~=0));
                                
                                
                                %look at displacement and amplitude ratio only (no
                                %directionality)
                                possibleLink = dispVecMag <= maxDispAllowed; %&& ...
                                
                                
                                %if this is a possible link ...
                                if possibleLink
                                    meanDispAllTracks =  nanmean(trackMeanDispMag);
                                    %calculate the cost of linking
                                    dispVecMag2 = dispVecMag ^ 2; %due to displacement
                                    
                                    ampDiffIndSpS = abs(ampSp - ampS);
                                    ampDiffIndSpSp1 = abs(ampSp - ampSp1);
                                    
                                    %ampCost = ampRatio; %due to amplitude
                                    %ampCost(ampCost<1) = ampCost(ampCost<1) ^ (-2); %punishment harsher when intensity of splitting feature < sum of intensities of features after splitting
                                    ampCost = 10^3*(ampDiffIndSpS+ampDiffIndSpSp1);
                                    meanDisp2Tracks = trackMeanDispMag(iStart); %for displacement scaling
                                    if isnan(meanDisp2Tracks)
                                        meanDisp2Tracks = nanmean(meanDispAllTracks);
                                    end
                                    cost12 = ampCost + dispVecMag2 / (meanDisp2Tracks^2);
                                    
                                    
                                    
                                    
                                    %if the lifetime consideration does not make this link impossible
                                    if ~isinf(cost12)
                                        
                                        %add this cost to the list of costs
                                        costMS(iPair) = cost12;
                                        
                                        %check whether the track being split from has had something
                                        %possibly splitting from it in this same frame
                                        prevAppearance = find(indxMSMS == iSplit);
                                        
                                        %if this track in this frame did not appear before ...
                                        %THIS SECTION IS OUTDATED
                                        if isempty(prevAppearance)
                                            
                                            %increase the "split index" by one
                                            numSplit = numSplit + 1;
                                            
                                            %save the splitting track's number
                                            indxMSMS(iPair) = iSplit;
                                            
                                            %store the location of this pair in the cost matrix
                                            indx1MS(iPair) = numSplit+numTracks; %row number
                                            indx2MS(iPair) = iStart; %column number
                                            
                                            %calculate the alternative cost of not splitting for the
                                            %track that the start is possibly splitting from
                                            
                                            %get the average square displacement in this track
                                            trackCoord = trackedFeatureInfo(indxMSMS(iPair),:);
                                            trackCoord = reshape(trackCoord',8,[]);
                                            if issparse(trackCoord)
                                                trackCoord = full(trackCoord);
                                                trackCoord(trackCoord==0) = NaN;
                                                %if probDim == 2
                                                trackCoord(3,:) = 0;
                                                trackCoord(7,:) = 0;
                                                %end
                                            end
                                            dispVecMag2 = (diff(trackCoord,1,2)).^2;
                                            dispVecMag2 = nanmean(dispVecMag2,2);
                                            dispVecMag2 = sum(dispVecMag2(1:probDim));
                                            
                                            %if the average square displacement is smaller
                                            %than resLimit^2, then expand it
                                            %  dispVecMag2 = max([dispVecMag2 resLimit^2]);
                                            
                                            %calculate intensity cost if no split happens
                                            ampCost = ampSp / ampSp1;
                                            ampCost(ampCost<1) = ampCost(ampCost<1) ^ (-2);
                                            
                                            %get track's mean displacement
                                            meanDisp1Track = trackMeanDispMag(indxMSMS(iPair));
                                            meanDisp1Track(isnan(meanDisp1Track)) = ...
                                                meanDispAllTracks;
                                            
                                            %calculate alternative cost
                                            cost12 = dispVecMag2 * ampCost ...
                                                / (meanDisp1Track^2);
                                            
                                            %although alternative cost is still calculated,
                                            %it is actually not used any more in the end
                                            
                                            %add this cost to the list of alternative costs
                                            altCostMS(iPair) = cost12;
                                            
                                        else %if this track in this frame appeared before
                                            
                                            %do not increase the "split index" or save the
                                            %splitting track's index (they are already saved)
                                            
                                            %store the location of this pair in the cost matrix
                                            indx1MS(iPair) = indx1MS(prevAppearance); %row number
                                            indx2MS(iPair) = iStart; %column number
                                            
                                            %no need to calculate and save the alternative cost
                                            %since that is already saved from previous appearance
                                            
                                        end %(if isempty(prevAppearance))
                                        
                                    end %(if ~isinf(cost12))
                                    
                                end %(if possibleLink)
                                
                            end %(for for iPair = 1 : numPairs)
                            
                            %keep only pairs that turned out to be possible
                            possiblePairs = find(indx1MS ~= 0);
                            indx1MS   = indx1MS(possiblePairs);
                            indx2MS   = indx2MS(possiblePairs);
                            costMS    = costMS(possiblePairs);
                            possibleSplits = find(indxMSMS ~= 0);
                            altCostMS = altCostMS(possibleSplits);
                            indxMSMS  = indxMSMS(possibleSplits);
                            clear possiblePairs possibleSplits
                            
                            %append these vectors to overall cost and related vectors
                            indx1 = [indx1, indx1MS'];
                            indx2 = [indx2, indx2MS'];
                            costs  = [costs, costMS'];
                            altCostSplit = [altCostSplit, altCostMS'];
                            indxSplit = [indxSplit, indxMSMS'];
                            
                        end %(for startTime = 2 : numFrames)
                        
                        
                    end
                    
                    
                    
                    
                    
                    costMat = sparse(indx1,indx2,costs,numTracks+numSplit,numTracks);
                    
                    costMat(isnan(costMat))=0;
                    costMat(isinf(costMat))=0;
                    
                    
                    
                    %link tracks based on this cost matrix, allowing for birth and death
                    [link12,link21] = lap(costMat,[],[], 1,100);
                    link12 = double(link12);
                    link21 = double(link21);
                    
                    %put the indices of all tracks from linking in one vector
                    tracks2Link = (1:numTracksLink)';
                    tracksRemaining = tracks2Link;
                    
                    %reserve memory space for matrix showing track connectivity
                    compoundTrack = zeros(numTracksLink,600);
                    
                    %initialize compTrackIndx
                    compTrackIndx = 0;
                    %
                    while ~isempty(tracksRemaining)
                        
                        %update compound track index by 1
                        compTrackIndx = compTrackIndx + 1;
                        
                        %take first track as a seed to build a compound track with
                        %closed gaps and merges/splits
                        trackSeed = tracksRemaining(1);
                        seedLength = 1;
                        seedLengthOld = 0; %dummy just to get into the while loop
                        
                        %while current seed contains more tracks than previous seed, i.e.
                        %whie new track segments are still being added to the compound
                        %track
                        while seedLength > seedLengthOld
                            
                            %store current seed for later comparison
                            seedLengthOld = seedLength;
                            
                            %find tracks connected to ends of seed tracks
                            tmpTracks = link12(trackSeed);
                            trackLink2End = tmpTracks(tmpTracks <= numTracksLink); %starts linked to ends
                            trackMerge = [];
                            %if mergeSplit
                            %    trackMerge = indxMerge(tmpTracks(tmpTracks > numTracksLink & ...
                            %        tmpTracks <= numTracksLink+numMerge) - numTracksLink); %tracks that ends merge with
                            %end
                            
                            %find tracks connected to starts of seed tracks
                            tmpTracks = link21(trackSeed);
                            trackLink2Start = tmpTracks(tmpTracks <= numTracksLink); %ends linked to starts
                            trackSplit = [];
                            if mergeSplit
                                trackSplit = indxSplit(tmpTracks(tmpTracks > numTracksLink & ...
                                    tmpTracks <= numTracksLink+numSplit) - numTracksLink); %tracks that starts split from
                            end
                            
                            %put all tracks together as the new seed
                            trackSeed = [trackSeed; trackLink2End; trackLink2Start; ...
                                trackMerge'; trackSplit'];
                            
                            %remove repetitions and arrange tracks in ascending order
                            trackSeed = unique(trackSeed);
                            
                            %get number of tracks in new seed
                            seedLength = length(trackSeed);
                            
                            %expand new seed if merging/splitting are allowed
                            if mergeSplit
                                
                                %variables storing merge/split seed tracks
                                mergeSeed = [];
                                splitSeed = [];
                                
                                %go over all seed tracks
                                for iSeed = 1 : seedLength
                                    
                                    %get the location(s) of this track in indxMerge
                                    %mergeSeed = [mergeSeed; find(indxMerge == trackSeed(iSeed))];
                                    
                                    %get the location(s) of this track in indxSplit
                                    splitSeed = [splitSeed; find(indxSplit == trackSeed(iSeed))'];
                                    
                                end
                                
                                %add numTracksLink to mergeSeed and splitSeed to determine
                                %their location in the cost matrix
                                %mergeSeed = mergeSeed + numTracksLink;
                                splitSeed = splitSeed + numTracksLink;
                                
                                %find tracks merging with seed tracks
                                trackMerge = [];
                                %for iSeed = 1 : length(mergeSeed)
                                %    trackMerge = [trackMerge; find(link12(1:numTracksLink)==mergeSeed(iSeed))];
                                %end
                                
                                %find tracks splitting from seed tracks
                                trackSplit = [];
                                for iSeed = 1 : length(splitSeed)
                                    trackSplit = [trackSplit; find(link21(1:numTracksLink)==splitSeed(iSeed))];
                                end
                                
                                %add these track to the seed
                                trackSeed = [trackSeed; trackMerge; trackSplit];
                                
                                %remove repetitions and arrange tracks in ascending order
                                trackSeed = unique(trackSeed);
                                
                                %get number of tracks in new seed
                                seedLength = length(trackSeed);
                                
                            end %(if mergeSplit)
                            
                        end %(while length(trackSeed) > length(trackSeedOld))
                        
                        %expand trackSeed to reserve memory for connetivity information
                        trackSeedConnect = [trackSeed zeros(seedLength,2)];
                        
                        %store the tracks that the ends of the seed tracks are linked to,
                        %and indicate whether it's an end-to-start link (+ve) or a merge (-ve)
                        tmpTracks = link12(trackSeed);
                        if mergeSplit
                            %tmpTracks(tmpTracks > numTracksLink & tmpTracks <= ...
                            %    numTracksLink+numMerge) = -indxMerge(tmpTracks(tmpTracks > ...
                            %    numTracksLink & tmpTracks <= numTracksLink+numMerge) - numTracksLink);
                        end
                        tmpTracks(tmpTracks > numTracksLink) = NaN;
                        trackSeedConnect(:,2) = tmpTracks;
                        
                        %store the tracks that the starts of the seed tracks are linked to,
                        %and indicate whether it's a start-to-end link (+ve) or a split (-ve)
                        tmpTracks = link21(trackSeed);
                        if mergeSplit
                            tmpTracks(tmpTracks > numTracksLink & tmpTracks <= ...
                                numTracksLink+numSplit) = -indxSplit(tmpTracks(tmpTracks > ...
                                numTracksLink & tmpTracks <= numTracksLink+numSplit) - numTracksLink);
                        end
                        tmpTracks(tmpTracks > numTracksLink) = NaN;
                        trackSeedConnect(:,3) = tmpTracks;
                        
                        %store tracks making up this compound track and their connectivity
                        compoundTrack(compTrackIndx,1:3*seedLength) = reshape(...
                            trackSeedConnect,3*seedLength,1)';
                        
                        %in the list of all tracks, indicate that these tracks have
                        %been taken care of by placing NaN instead of their number
                        tracks2Link(trackSeed) = NaN;
                        
                        %retain only tracks that have not been linked to anything yet
                        tracksRemaining = tracks2Link(~isnan(tracks2Link));
                        
                    end %(while ~isempty(tracksRemaining))
                    
                    %remove empty rows
                    maxValue = max(compoundTrack,[],2);
                    compoundTrack = compoundTrack(maxValue > 0,:);
                    
                    %determine number of tracks after gap closing (including merge/split if
                    %specified)
                    numTracksCG = size(compoundTrack,1);
                    
                    %reserve memory for structure storing tracks after gap closing
                    tracksFinal = repmat(struct('tracksFeatIndxCG',[],...
                        'tracksCoordAmpCG',[],'seqOfEvents',[]),numTracksCG,1);
                    
                    f = waitbar(0,'Closing gaps...');
                    
                    %go over all compound tracks
                    for iTrack = 1 : numTracksCG
                        if ~rem(iTrack,100)
                            waitbar(iTrack/numTracksCG,f,'Closing gaps...');
                        end
                        %get indices of tracks from linking making up current compound track
                        %determine their number and connectivity
                        trackSeedConnect = compoundTrack(iTrack,:)';
                        trackSeedConnect = trackSeedConnect(trackSeedConnect ~= 0);
                        seedLength = length(trackSeedConnect)/3; %number of segments making current track
                        trackSeedConnect = reshape(trackSeedConnect,seedLength,3);
                        
                        %get their start times
                        segmentStartTime = trackStartTime(trackSeedConnect(:,1));
                        
                        %arrange segments in ascending order of their start times
                        [segmentStartTime,indxOrder] = sort(segmentStartTime);
                        trackSeedConnect = trackSeedConnect(indxOrder,:);
                        
                        %get the segments' end times
                        segmentEndTime = trackEndTime(trackSeedConnect(:,1));
                        
                        %calculate the segments' positions in the matrix of coordinates and
                        %amplitudes
                        segmentStartTime8 = 8 * (segmentStartTime - 1) + 1;
                        segmentEndTime8   = 8 * segmentEndTime;
                        
                        %instead of having the connectivity in terms of the original track
                        %indices, have it in terms of the indices of this subset of tracks
                        %(which are arranged in ascending order of their start times)
                        for iSeed = 1 : seedLength
                            value = trackSeedConnect(iSeed,2);
                            if value > 0
                                trackSeedConnect(iSeed,2) = find(trackSeedConnect(:,1) == ...
                                    value);
                            elseif value < 0
                                trackSeedConnect(iSeed,2) = -find(trackSeedConnect(:,1) == ...
                                    -value);
                            end
                            value = trackSeedConnect(iSeed,3);
                            if value > 0
                                trackSeedConnect(iSeed,3) = find(trackSeedConnect(:,1) == ...
                                    value);
                            elseif value < 0
                                trackSeedConnect(iSeed,3) = -find(trackSeedConnect(:,1) == ...
                                    -value);
                            end
                        end
                        
                        %get track information from the matrices storing linking information
                        tracksFeatIndxCG = trackedFeatureIndx(trackSeedConnect(:,1),:);
                        tracksCoordAmpCG = trackedFeatureInfo(trackSeedConnect(:,1),:);
                        
                        %convert zeros to NaNs where approriate for the case of sparse
                        %matrices
                        if issparse(tracksCoordAmpCG)
                            
                            %convert sparse to full
                            tracksCoordAmpCG = full(tracksCoordAmpCG);
                            
                            %go over all the rows in this compound track
                            for iRow = 1 : size(tracksCoordAmpCG,1)
                                
                                %find all the zero entries
                                colZero = find(tracksCoordAmpCG(iRow,:)==0);
                                colZero = colZero(:)';
                                
                                %find the columns of the x-coordinates corresponding to
                                %the zero columns
                                xCoordCol = colZero - mod(colZero-1,8*ones(size(colZero)));
                                
                                %keep only the columns whose x-coordinate is zero as
                                %well
                                colZero = colZero(tracksCoordAmpCG(iRow,xCoordCol)==0);
                                
                                %replace zero with NaN in the surviving columns
                                tracksCoordAmpCG(iRow,colZero) = NaN;
                                
                            end
                            
                        end
                        
                        %perform all gap closing links and modify connectivity accordingly
                        %go over all starts in reverse order
                        for iSeed = seedLength : -1 : 2
                            
                            %find the track this track might be connected to
                            track2Append = trackSeedConnect(iSeed,3);
                            
                            %if there is a track (which is not a split)
                            if track2Append > 0
                                
                                %put track information in the relevant row
                                tracksFeatIndxCG(track2Append,segmentStartTime(iSeed):...
                                    segmentEndTime(iSeed)) = tracksFeatIndxCG(iSeed,...
                                    segmentStartTime(iSeed):segmentEndTime(iSeed));
                                tracksFeatIndxCG(iSeed,:) = 0;
                                tracksCoordAmpCG(track2Append,segmentStartTime8(iSeed):...
                                    segmentEndTime8(iSeed)) = tracksCoordAmpCG(iSeed,...
                                    segmentStartTime8(iSeed):segmentEndTime8(iSeed));
                                tracksCoordAmpCG(iSeed,:) = NaN;
                                
                                %update segment information
                                segmentEndTime(track2Append) = segmentEndTime(iSeed);
                                segmentEndTime8(track2Append) = segmentEndTime8(iSeed);
                                segmentEndTime(iSeed) = NaN;
                                segmentEndTime8(iSeed) = NaN;
                                segmentStartTime(iSeed) = NaN;
                                segmentStartTime8(iSeed) = NaN;
                                
                                %update connectivity
                                trackSeedConnect(track2Append,2) = trackSeedConnect(iSeed,2);
                                trackSeedConnect(trackSeedConnect(:,2) == iSeed,2) = track2Append;
                                trackSeedConnect(trackSeedConnect(:,3) == iSeed,3) = track2Append;
                                trackSeedConnect(trackSeedConnect(:,2) == -iSeed,2) = -track2Append;
                                trackSeedConnect(trackSeedConnect(:,3) == -iSeed,3) = -track2Append;
                                
                            end %(if track2Append > 0)
                            
                        end %(for iSeed = seedLength : -1 : 2)
                        
                        %find rows that are not empty
                        maxValue = max(tracksFeatIndxCG,[],2);
                        rowsNotEmpty = find(maxValue > 0);
                        
                        %remove empty rows
                        tracksFeatIndxCG = tracksFeatIndxCG(rowsNotEmpty,:);
                        tracksCoordAmpCG = tracksCoordAmpCG(rowsNotEmpty,:);
                        segmentEndTime   = segmentEndTime(rowsNotEmpty);
                        segmentStartTime = segmentStartTime(rowsNotEmpty);
                        trackSeedConnect = trackSeedConnect(rowsNotEmpty,:);
                        
                        %update connectivity accordingly
                        %by now, only merges and splits are left - thus no need for minus
                        %sign to distinguish them from closed gaps
                        for iSeed = 1 : length(rowsNotEmpty)
                            trackSeedConnect(trackSeedConnect(:,2) == -rowsNotEmpty(...
                                iSeed),2) = iSeed;
                            trackSeedConnect(trackSeedConnect(:,3) == -rowsNotEmpty(...
                                iSeed),3) = iSeed;
                        end
                        
                        %determine new "seedLength"
                        seedLength = length(rowsNotEmpty);
                        
                        %store the sequence of events of this track
                        seqOfEvents = [segmentStartTime ones(seedLength,1) ...
                            (1:seedLength)' trackSeedConnect(:,3); ...
                            segmentEndTime 2*ones(seedLength,1) ...
                            (1:seedLength)' trackSeedConnect(:,2)];
                        
                        %sort sequence of events in ascending order of time
                        [tmp,indxOrder] = sort(seqOfEvents(:,1));
                        seqOfEvents = seqOfEvents(indxOrder,:);
                        
                        %add 1 to the times of merges
                        indx = find(~isnan(seqOfEvents(:,4)) & seqOfEvents(:,2) == 2);
                        seqOfEvents(indx,1) = seqOfEvents(indx,1) + 1;
                        
                        %find the frame where the compound track starts and the frames
                        %where it ends
                        frameStart = seqOfEvents(1,1);
                        frameEnd   = seqOfEvents(end,1);
                        
                        %store final tracks, removing frames before anything happens and
                        %after everything happens
                        tracksFinal(iTrack).tracksFeatIndxCG = tracksFeatIndxCG(:,...
                            frameStart:frameEnd);
                        tracksFinal(iTrack).tracksCoordAmpCG = tracksCoordAmpCG(:,...
                            8*(frameStart-1)+1:8*frameEnd);
                        tracksFinal(iTrack).seqOfEvents = seqOfEvents;
                        
                    end %(for iTrack = 1 : numTracksCG)
                    close(f)
                    %%
                    
                    a = arrayfun(@(x) size(x.tracksFeatIndxCG,2),tracksFinal);
                    longtracksFinal = tracksFinal(a>25);% remove all tracks shorter than 10 frames
                end
            else
                longtracksFinal = [];
            end
            
            R.setTracks(longtracksFinal,pos);
        end
        
        
        
        
        function R = merge(Rvec,varargin)
            % TODO: extend support for merging object with different
            % TimeVec structures, for now I'm assuming they are the same!
            arg.prefix='%g_';
            arg = parseVarargin(varargin,arg);
            arg.Results = MultiPositionSingleCellWoundResults;
            R = merge@MultiPositionResults(Rvec,arg);
            
            %R.PIVlbl=PIVLbl.empty(1,0);
            %for i=1:numel(Rvec)
            %    R.PIVlbl = [R.PIVlbl Rvec(i).PIVlbl];
            %end
            
            allfields = {};
            for i=1:numel(Rvec)
                allfields = union(allfields,fieldnames(Rvec(i).TimeVecs));
            end
            for i=1:numel(Rvec)
                T=Rvec(i).TimeVecs;
                missingflds = setdiff(allfields,fieldnames(T));
                for j=1:numel(missingflds)
                    T(1).(missingflds{j})=[];
                end
                T=orderfields(T,allfields);
                if i==1
                    TimeVec=T;
                else
                    TimeVec=[TimeVec T];  %#ok<AGROW>
                end
            end
            R.TimeVecs=TimeVec;
        end
        
        
        
        function h = plotCellVsTime(R,pos)
            WellCells = R.getWellsLbl(pos);
            plot((0:(numel(R.WellLbls{1})-1))/2,cellfun(@(x) x.num, WellCells));
            ylabel('#Cells')
            xlabel('Time(h)')
            shg
        end
        
        
        
        
        %plotting
        %         function HeatMapData(R, dataname,pos, frame, varargin)
        %             arg.clims = [-15, 15];
        %             arg = parseVarargin(varargin,arg);
        %
        %             a = R.getData(dataname,pos);
        %             a = a(:,frame);
        %                       if min(arg.clims)>=0
        %                           colormap(magma())
        %                       else
        %             colormap(makeColorMap([0.6 0 0.6],[0 0 0],[0.8 0.8 0]))
        %                       end;
        %             imagesc(unique(R.PIVlbl{1}.X), unique(R.PIVlbl{1}.Y), reshape(a,numel(unique(R.PIVlbl{1}.X)),numel(unique(R.PIVlbl{1}.Y)))',arg.clims);
        %             set(gcf,'color','w');
        %             axis equal
        %             title(dataname)
        %             colorbar
        %             set(gca,'ydir','normal')
        %             shg;
        %         end
        
        
        function h = PlotDisp(R,pos, i, dt,varargin)
            
            %function to plot displacement bw frame i and i+dt
            WellCells = R.getWellsLbl(pos);
            CentroidsD = WellCells{i}.Centroids;
            CentroidsA = WellCells{i+dt}.Centroids;
            
            
            CM = double(WellCells{i}.Link12Mat);%continuous connectivity map
            for j=i+1:i+dt-1
                CM=CM*double(WellCells{j}.Link12Mat);
            end
            [indx1, indx2] = find(CM);
            
            [indx1, J] = sort(indx1);
            
            idxToPlot=J;
            
            
            indx2 = indx2(idxToPlot);
            CentroidsD = CentroidsD(indx1,:);
            CentroidsA = CentroidsA(indx2,:);
            
            range = 1:length(indx1);
            if nargin>4
                range = varargin{1};
                if range(end)>indx1(end)
                    range = range(1):indx1(end);
                    disp('range end out of bounds, drawing all points up to the total # of cells.')
                end
                if range(1)>indx1(end)
                    range = 1:length(indx1);
                    disp('range fully out of bounds, drawing all points.')
                end
            end;
            
            %tzeva = viridis(length(indx1));
            %scatter3(CentroidsD(range,1),CentroidsD(range,2),-CentroidsD(range,3),[],tzeva(range,:),'*');
            %hold on
            %scatter3(CentroidsA(range,1),CentroidsA(range,2),-CentroidsA(range,3),[],tzeva(range,:));
            dx=CentroidsA(range,1)-CentroidsD(range,1);
            dy=CentroidsA(range,2)-CentroidsD(range,2);
            dz=CentroidsA(range,3)-CentroidsD(range,3);
            
            h = quiver3(CentroidsD(range,1),CentroidsD(range,2),-CentroidsD(range,3),dx,dy,-dz, 0);
            
            
            set(gca,'xlim',[-200 3000],'ylim',[-200 3000],'CameraPositionMode','manual','CameraPosition',[-2.6364e+03 -1.4045e+04 2.1889e+03])
            %hold on
            shg
        end
        
        
        
        %% Some sanity checks
        function allTracks = allTrackMatrix(R,pos)
            longtracksFinal = R.getTracks(pos);
            frames = R.Frames;
            allTracks = zeros(size(longtracksFinal,1),numel(frames));
            nTracks = 1;
            for i=1:size(longtracksFinal,1);
                trackStart = longtracksFinal(i).seqOfEvents(1,1);
                trackEnd = longtracksFinal(i).seqOfEvents(end,1);
                for j=1:size(longtracksFinal(i).tracksFeatIndxCG,1)
                    allTracks(nTracks,trackStart:trackEnd) = longtracksFinal(i).tracksFeatIndxCG(j,:);
                    nTracks = nTracks+1;
                end
            end
        end
        
        function allTracks = heatmapTracks(R,pos, dataname)
            longtracksFinal = R.getTracks(pos);
            frames = R.Frames;
            allTracks = zeros(size(longtracksFinal,1),numel(frames));
            for i=1:size(longtracksFinal,1);
                trackStart = longtracksFinal(i).seqOfEvents(1,1);
                trackEnd = longtracksFinal(i).seqOfEvents(end,1);
                allTracks(i,trackStart:trackEnd) = nanfill(longtracksFinal(i).(dataname));
            end
        end
        
        
        
        % %         function  plotInfectedOverT(R,pos)
        %             WellCells = R.getWellsLbl(R.PosNames{WellNum});
        %             ThreshVirus = median(WellCells{1}.VirusIntensities)+3*std(WellCells{1}.VirusIntensities);
        %
        %                 {1}.VirusIntensities
        %
        %         end
        
        function h = plotFracInTracks(R,pos)
            h = plot(sum(allTrackMatrix(R,pos)>0)'./cellfun(@(x) x.num, R.getWellsLbl(pos)));
            xlabel('Frame')
            ylabel({'Fraction of detected cells' 'in an active track'})
            shg
        end
        
        function histTrackLength(R,pos)
            hist(sum(allTrackMatrix(R,pos)>0,2),100);
            xlabel('Length (frames)')
            ylabel('#')
            shg
        end
        
        function plotCompTrack(R,pos,trackNum)
            
            PL = getTracks(R,pos);
            plotCompTrack(PL(trackNum))
            
            Ch = get(gcf,'Children');
            
            Ch(1).YLabel.String = 'Hoechst (a.u.)';
            Ch(2).YLabel.String = 'HSV (a.u.)';
            Ch(1).YLim = [0 .8];
            Ch(2).YLim = [.05 .15];
            
        end
        function J = TracksThatGetInfected(R,pos)
            PL = getTracks(R,pos);
            J = find([PL.CellsGetInfected]);
        end
        
        function J = TracksThatDie(R,pos)
            PL = getTracks(R,pos);
            J = find([PL.CellDies]);
        end
        
        function paramMat = putInMat(R,Tracks,Param)
            frames = R.Frames;
            assert(isfield(Tracks,Param),[Param ' is not a string in Tracks'])
            dt = (numel(Tracks(1).T)-numel(Tracks(1).(Param)))/2;
            paramMat = nan(numel(Tracks),numel(frames));
            for i=1:numel(Tracks)
                paramMat(i,Tracks(i).T(dt+1:end-dt)) = Tracks(i).(Param);
            end
        end
        
        function plotTrackInIJ(R,pos,i)
            T = R.getTracks(pos);
            Track = T(i);
            Xs = Track.tracksCoordAmpCG(1,1:8:end);
            Ys = Track.tracksCoordAmpCG(1,2:8:end);
            
            v = [Xs(1:end-1); Ys(1:end-1); Xs(2:end); Ys(2:end)];
            plotLinesInIJ(v(:))
        end
    end
    
    
    
end