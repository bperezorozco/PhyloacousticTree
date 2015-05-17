classdef WAVRec < handle
    % WAVRec < handle
    %
    % Uses histcn.m as provided by Bruno Luong in Mathworkd File Exchange
    % (see http://www.mathworks.co.uk/matlabcentral/fileexchange/23897-n-dimensional-histogram)
    %
    % TODO: (Related to NOTE: in computeSpectralFeatures method) Check what happens if conditioned spectrogram contains
    % NaNs (because the InterestingFrameFraction property of the used FEATHISTParams object is set to lower than 1).
    %
    % TODO: check if 6D histogram computation is a bottlenec and if so consider changing the use of histcn method with
    % the ndhistc code provided by Kangwon Lee in Mathworkd File Exchange
    % (see http://www.mathworks.com/matlabcentral/fileexchange/3957-ndhistc)
    
    %% WAVRec PlantUML Class Diagram
    %{
    @startuml
    skinparam classAttributeIconSize 0
    class WAVRec {
        #+  WavDataPath;
        #+  SamplingFrequency;
        #+  NumofSamplesPerChannel;
        #+  NumofChannels;
        #+  DurationSec;
        #+  BitDepth;
        +   WAVRec(_) : this
        -   addAudioInfoProperties(this)
        +   computeSpectrogram(FH)
        +   computeFeatures(this,_,FH)
        +   computeSlidingWindowSFM(this,FH,_)
        +   computeSlidingWindowKurtosis(this,FH)
        +   computeSlidingWindowNegEntropy(this,FH)
        +   checkAudioReadConsistency(this)
        +   <<Static>> computeHistogram(featurearray,FH)
        +   <<Static>> conditionSpectrogram(_,_,FH)
        -   <<Static>> histcn(_,varargin)
    }
    @enduml
    %}
    %% Properties
    properties (SetAccess = protected, GetAccess = public)
        WavDataPath;
        SamplingFrequency;
        NumofSamplesPerChannel;
        NumofChannels;
        DurationSec;
        BitDepth;
    end % properties
    %% Methods
    methods %Constructor
        function this = WAVRec(wavdatapath)
            this.WavDataPath = wavdatapath;
            this.addAudioInfoProperties;
        end % WAVRec
    end
    methods (Access = private)
        function this = addAudioInfoProperties(this)
            temp = audioinfo(this.WavDataPath);
            this.SamplingFrequency = temp.SampleRate;
            this.NumofSamplesPerChannel = temp.TotalSamples;
            this.NumofChannels = temp.NumChannels;
            this.DurationSec = temp.Duration;
            if isfield(temp,'BitsPerSample')
                this.BitDepth = temp.BitsPerSample;
            else
                this.BitDepth = NaN;
            end
        end % addAudioInfoProperties
    end
    methods
        function [spec, F, T, framelimssamp] = computeSpectrogram(this, paramsobj)
            %Computes unnormalised spectrogram (no frame scaling to pmf) for frequencies from 0 to Fs/2
            %
            % [spec, F, T, framelimssamp] = WAVRec.computeSpectrogram(paramsobj)
            %
            % The first output is the absolute value of the spectrogram with increasing frequencies in increasing row
            % indices and increasing frame time in increasing column indices
            % The second and third outputs are the vectors of frequencies (in Hz) and frame center times (in sec) as
            % specified in the documentation of the spectrogram Matlab function 
            % The fourth outpus is a two columnn array containing the limits of the frames in samples
            %
            % The input is a handle object of the class FEATHISTParams
            
            %%
            % Read wav file
            temp = audioread(this.WavDataPath);
            % Get the channel and time interval that are specified in paramsobj
            startind = 1 + floor( paramsobj.AudioLimits(1) * size(temp,1) );
            endind = round( paramsobj.AudioLimits(2) * size(temp,1) );
            x = temp( startind:endind, paramsobj.AudioChannel).';
            % Prepare spectrogram computation parameters
            tempwindow = window(paramsobj.WindowType, round(paramsobj.WindowLengthms*this.SamplingFrequency/1000)).';
            tempnoverlap = round(paramsobj.OverlapFraction*length(tempwindow));
            switch round(this.SamplingFrequency)
                case 96000
                    tempnfft = paramsobj.PointsFFT96k;
                case 48000
                    tempnfft = paramsobj.PointsFFT48k;
                case 44100
                    tempnfft = paramsobj.PointsFFT44k;
                case 32000
                    tempnfft = paramsobj.PointsFFT32k;
                case 22050
                    tempnfft = paramsobj.PointsFFT22k;
                case 16000
                    tempnfft = paramsobj.PointsFFT16k;
                otherwise
                    error('problem in the determination of nfft based on the sampling frequency')
            end
            % Compute spectrogram for all frequencies from 0 to Fs/2
            [tempspec, F, T] = spectrogram(x, tempwindow, tempnoverlap, tempnfft, this.SamplingFrequency);
            spec = abs(tempspec);
            %% Determine frame limits in samples
            % (also check that that this is consistent with the spec variable)
            framelensamp = round(paramsobj.WindowLengthms*this.SamplingFrequency/1000);
            framestepsamp = round((1 - paramsobj.OverlapFraction-eps) * framelensamp);
            frameendsamp = ( framelensamp : framestepsamp : length(x) );
            framestartsamp = frameendsamp - framelensamp + 1;
            framelimssamp = [framestartsamp(:) frameendsamp(:)];
            if size(spec,2) ~= size(framelimssamp,1)
                % error('problem in reverse engineering of spectrogram computation frame limits')
                warning('WAVRec_computeSpectrogram:NumLimsFramesDiscrepancy', ...
                    ['discrepancy in #frames/frame limits for mdid ' num2str(this.MetadataID)]);
%                 warning('WAVRec:frames', 'somemsg')
            end
        end % computeSpectrogram
        function [featurearray, framelimssamp, framelimsms] = computeFeatures(this, paramsobj)
            %Computes features
            %
            % featurearray = WAVRec.computeFeatures(paramsobj)
            %
            % First input argument is a handle object of the class FEATHISTParams
            %
            % The first output is an array with as many rows as there are feature types in .FeaturesType property of the
            % FEATHISTArray object input argument and as many columns as there are frames (determined according to the
            % .WindowLengthms and .OverlapFraction properties of the FEATHISTArray object in the second input argument
            % and Matlab's spectorgram function.
            %
            % The second and third outputs are the limits of the frames (start : end) in samples and ms correspondingly
            
            %% Compute spectogram and determine frame limits in samples and ms
            % NOTE: This computation is done even if it is not necessary, i.e. if the check
            % any(ismember(paramsobj.FeaturesType,{'mean', 'std', 'median', 'IQR', 'normIQR', 'normIQRold',
            % 'MutualInformation', 'SFM', 'SFM_Madhu', 'SFM_Madhu_2_8'})) evaluates to 0. This is done to get the time
            % and frequency indices and could be replicated by including here a reverse engineering of the frame slicing
            % and frequnecy points determination of Matlab's spectrogram
            [spec, F, T, framelimssamp] = this.computeSpectrogram(paramsobj);
            [X, xf] = this.conditionSpectrogram(spec, F, paramsobj);
            notallnaninds = ~all(isnan(X));
            xf = xf(:).';
            framelimsms = repmat(T(:)*1000,[1,2]) + ...
                repmat([-paramsobj.WindowLengthms/2 paramsobj.WindowLengthms/2],[length(T),1]);
            framelensamp = framelimssamp(1,2)-framelimssamp(1,1)+1;
            framestepsamp = framelimssamp(2,1)-framelimssamp(1,1);
            %% Compute features
            % Initialise featurearray
            featurearray = zeros( length(paramsobj.FeaturesType), length(T) );
            % Compute each line of featurearray (replicated computations for some features, could include extra 'if'
            % clauses to compute what is necessary only once)
            if any(ismember({'median','IQR','normIQRold'},paramsobj.FeaturesType))
                tempcdf = cumsum(X);
                nonnantempcdf = tempcdf(:,notallnaninds);
                nonnanfeaturearray = zeros(1,size(nonnantempcdf,2));
                
            end
            if any(ismember({'IQR','normIQRold'},paramsobj.FeaturesType))
                tempfeatarr1 = nan( 1, size(nonnantempcdf,2) );
                tempfeatarr2 = nan( 1, size(nonnantempcdf,2) );
            end
            if any(ismember({'HighestPeakFrequency','HighestPeak3dBWidth','SecondPeakFrequency''SecondPeak3dBWidth'},...
                    paramsobj.FeaturesType))
                [~,tempsortXind] = sort(X,1,'descend');
            end
            for ii = 1:length(paramsobj.FeaturesType)
                switch paramsobj.FeaturesType{ii}
                    case 'mean'
                        featurearray(ii,:) = xf*X;
                    case 'std'
                        meantemp = xf*X;
                        featurearray(ii,:) = sqrt(sum(...
                            ((repmat(xf.',1,size(X,2))-repmat(meantemp,length(xf),1)).^2).*X...
                            ));
                    case 'var'
                        meantemp = xf*X;
                        featurearray(ii,:) = sum(...
                            ((repmat(xf.',1,size(X,2))-repmat(meantemp,length(xf),1)).^2).*X...
                            );
                    case 'skew'
                        meantemp = xf*X;
                        featurearray(ii,:) = sum(...
                            ((repmat(xf.',1,size(X,2))-repmat(meantemp,length(xf),1)).^3).*X...
                            );
                    case 'kurt'
                        meantemp = xf*X;
                        featurearray(ii,:) = sqrt(sum(...
                            ((repmat(xf.',1,size(X,2))-repmat(meantemp,length(xf),1)).^4).*X...
                            ));
                    case 'median'
                        percval = 0.5;
                        boundinds = nan(2,size(nonnantempcdf,2));
                        for jj = 1:size(nonnantempcdf,2)
                            boundinds(:,jj) = [find(nonnantempcdf(:,jj)-percval<=0,1,'last') ; ...
                                               find(nonnantempcdf(:,jj)-percval>=0,1,'first')];
                        end
                        interpinds = boundinds(2,:)-boundinds(1,:)~=0;
                        noninterpinds = boundinds(2,:)-boundinds(1,:)==0;
                        templinind = (1:size(nonnantempcdf,2));
                        tempnum = (xf(boundinds(2,interpinds)) - xf(boundinds(1,interpinds)));
                        tempdenom1 = nonnantempcdf(sub2ind(size(nonnantempcdf),boundinds(1,interpinds),templinind(interpinds)));
                        tempdenom2 = nonnantempcdf(sub2ind(size(nonnantempcdf),boundinds(2,interpinds),templinind(interpinds)));
                        tempdenom = tempdenom2-tempdenom1;
                        nonnanfeaturearray(interpinds) = xf(boundinds(1,interpinds))+...
                            ((percval-tempdenom1).*(tempnum./tempdenom));
                        nonnanfeaturearray(noninterpinds) = xf(boundinds(1,noninterpinds));
                        featurearray(ii,:) = nan;
                        featurearray(ii,notallnaninds) = nonnanfeaturearray;
                    case 'IQR'
                        percval = 0.25;
                        boundinds = nan(2,size(nonnantempcdf,2));
                        for jj = 1:size(nonnantempcdf,2)
                            boundinds(:,jj) = [find(nonnantempcdf(:,jj)-percval<=0,1,'last') ; ...
                                               find(nonnantempcdf(:,jj)-percval>=0,1,'first')];
                        end
                        interpinds = boundinds(2,:)-boundinds(1,:)~=0;
                        noninterpinds = boundinds(2,:)-boundinds(1,:)==0;
                        templinind = (1:size(nonnantempcdf,2));
                        tempnum = (xf(boundinds(2,interpinds)) - xf(boundinds(1,interpinds)));
                        tempdenom1 = nonnantempcdf(sub2ind(size(nonnantempcdf),boundinds(1,interpinds),templinind(interpinds)));
                        tempdenom2 = nonnantempcdf(sub2ind(size(nonnantempcdf),boundinds(2,interpinds),templinind(interpinds)));
                        tempdenom = tempdenom2-tempdenom1;
                        tempfeatarr1(1,interpinds) = xf(boundinds(1,interpinds))+...
                            ((percval-tempdenom1).*(tempnum./tempdenom));
                        tempfeatarr1(1,noninterpinds) = xf(boundinds(1,noninterpinds));
                        percval = 0.75;
                        boundinds = nan(2,size(nonnantempcdf,2));
                        for jj = 1:size(nonnantempcdf,2)
                            boundinds(:,jj) = [find(nonnantempcdf(:,jj)-percval<=0,1,'last') ; ...
                                               find(nonnantempcdf(:,jj)-percval>=0,1,'first')];
                        end
                        interpinds = boundinds(2,:)-boundinds(1,:)~=0;
                        noninterpinds = boundinds(2,:)-boundinds(1,:)==0;
                        templinind = (1:size(nonnantempcdf,2));
                        tempnum = (xf(boundinds(2,interpinds)) - xf(boundinds(1,interpinds)));
                        tempdenom1 = nonnantempcdf(sub2ind(size(nonnantempcdf),boundinds(1,interpinds),templinind(interpinds)));
                        tempdenom2 = nonnantempcdf(sub2ind(size(nonnantempcdf),boundinds(2,interpinds),templinind(interpinds)));
                        tempdenom = tempdenom2-tempdenom1;
                        tempfeatarr2(1,interpinds) = xf(boundinds(1,interpinds))+...
                            ((percval-tempdenom1).*(tempnum./tempdenom));
                        tempfeatarr2(1,noninterpinds) = xf(boundinds(1,noninterpinds));
                        featurearray(ii,:) = nan;
                        featurearray(ii,notallnaninds) = tempfeatarr2 - tempfeatarr1;
                    case 'normIQRold'
                        for jj = 1:size(X,2)
                            temppctl = interp1(tempcdf(:,jj),xf,[0.25 0.5 0.75],'linear','extrap');
                            featurearray(ii,jj) = (temppctl(3) - temppctl(1))/temppctl(2);
                        end
                    case 'normIQR'
                        if all(ismember({'median','IQR'},paramsobj.FeaturesType)) && ...
                                (find(strcmp('median',paramsobj.FeaturesType)) < ...
                                find(strcmp('normIQR',paramsobj.FeaturesType))) && ...
                                (find(strcmp('IQR',paramsobj.FeaturesType)) < ...
                                find(strcmp('normIQR',paramsobj.FeaturesType)))
                            indmedian = strcmp('median',paramsobj.FeaturesType);
                            indIQR = strcmp('IQR',paramsobj.FeaturesType);
                            featurearray(ii,:) = featurearray(indIQR,:)./featurearray(indmedian,:);
                        else
                            error('include ''median'' and ''IQR'' in FEATHISTParams.FeaturesType before ''normIQR''')
                        end
                    case 'SFM' %standard
                        % NOTE: spectrogram normalised to pmf is used (becasue of invariability to mult.factor)
                        tempcondspec = X;
                        tempcondspec(tempcondspec==0) = eps;
                        featurearray(ii,:) = ...
                            (prod(tempcondspec,1).^(1/size(tempcondspec,1)))./(sum(tempcondspec,1)/size(tempcondspec,1));
                    case 'SFM_Madhu' %madhu
                        % NOTE: spectrogram normalised to pmf is used (becasue this is how the measure is described in
                        % the Madhu paper)
                        tempcondspec = X;
                        tempcondspec(tempcondspec==0) = eps;
                        featurearray(ii,:) = 2.^(-sum(tempcondspec.*log2(tempcondspec),1)/log2(size(tempcondspec,1)))-1;
                    case 'SFM_Madhu_2_8' %madhu
                        % NOTE: spectrogram normalised to pmf is used (becasue this is how the measure is described in
                        % the Madhu paper)
%                         tempcondspec = X;
                        paramsobjtemp = paramsobj.deepcopy;
                        paramsobjtemp.setFrequencyLimitsParams([2e3 8e3]);
                        [Xtemp, ~] = this.conditionSpectrogram(spec, F, paramsobjtemp);
                        Xtemp(Xtemp==0) = eps;
                        featurearray(ii,:) = 2.^(-sum(Xtemp.*log2(Xtemp),1)/log2(size(Xtemp,1)))-1;    
                    case 'MutualInformation'
                        error('not implemented yet')
                    case 'Kurtosis'
                        wavdata = audioread(this.WavDataPath);
                        wavdata = buffer(wavdata, framelensamp, framelensamp-framestepsamp, 'nodelay');
                        wavdata = wavdata(:, 1:length(T));
                        featurearray(ii,:) = kurtosis(wavdata,0);
                    case 'NegEntropy'
                        wavdata = audioread(this.WavDataPath);
                        wavdata = buffer(wavdata, framelensamp, framelensamp-framestepsamp, 'nodelay');
                        wavdata = wavdata(:, 1:length(T));
                        % Compute NegEntropy (this could be vectorised if it becomes a speed bottleneck)
                        negent = nan(1,length(T));
                        for jj = 1:length(negent)
                            framemean = mean(wavdata(:,jj));
                            framestd = std(wavdata(:,jj));
                            histlims = [-10*framestd 10*framestd] + framemean;
                            binwidth = .01;
                            wavdatahist = hist(wavdata(:,jj),(histlims(1):binwidth:histlims(2)));
                            wavdatahist(wavdatahist == 0) = eps;
                            wavdatahist = wavdatahist/sum(wavdatahist.*(binwidth*ones(1,length(wavdatahist))));
                            wavdataentr = -sum(wavdatahist.*log(wavdatahist).*(binwidth*ones(1,length(wavdatahist))));
                            normcomp = framestd * randn(framelensamp,1) + framemean;
                            normhist = hist(normcomp,(histlims(1):binwidth:histlims(2)));
                            normhist(normhist == 0) = eps;
                            normhist = normhist/sum(normhist.*(binwidth*ones(1,length(normhist))));
                            normentr = -sum(normhist.*log(normhist).*(binwidth*ones(1,length(normhist))));
                            negent(jj) = normentr - wavdataentr;
                        end
                        featurearray(ii,:) = negent;
                    case 'HighestPeakFrequency'
                        featurearray(ii,:) = nan;
                        featurearray(ii,notallnaninds) = xf(tempsortXind(1,notallnaninds));
                    case 'HighestPeak3dBWidth'
%                         above3dB = X > repmat(0.5*xf(tempsortf(1,:)),[size(X,1),1]);
                        for jj = 1:size(X,2)
                            above3dB = X(:,jj) > 0.5*X(tempsortXind(1,jj));
                            templowind = find(~above3dB(1:tempsortXind(1,jj)),1,'last');
                            if ~isempty(templowind)
                                templowfreq = xf(templowind);
                            else
                                templowfreq = xf(1);
                            end
                            temphighind = tempsortXind(1,jj)-1+find(~above3dB(tempsortXind(1,jj):end),1,'first');
                            if ~isempty(temphighind)
                                templhighfreq = xf(temphighind);
                            else
                                templhighfreq = xf(end);
                            end
                            featurearray(ii,jj) = templhighfreq-templowfreq;
                        end
                end
            end
        end % computeFeatures
        function [SFM, slidewinlimssamp, slidewinlimsms] = computeSlidingWindowSFM(this, paramsobj, sfmtype)
            % Computes SFM of wavrecobj for sliding windows as specified by the .WindowLengthms and .OverlapFraction
            % properties of the provided FEATHISTParams object (first input). Specifying these properties equal to the
            % FEATHISTParams object used for the computation of the spectrogram will thus give 'per frame' SFM. Otherwise
            % these properties can copy those specified in the 'sliding window' histogram computation of the
            % .appendSlidingWindow method of the HISTArray class (with proper conversion of overlap fraction to stepms).
            % The second input must be a string equal to either 'madhu' or 'starndard'
            %
            % The first output is a vector of SFM values per sliding window
            % The second output is a two column matrix of sliding window limits in samples
            % The third output is a two column matrix of sliding window limits in ms
            
            if ~sum(strcmp(sfmtype,{'standard','madhu'}))
                error('the optional third argument should be either ''standard'' or ''madhu''')
            end
            spec = this.computeSpectrogram(paramsobj);
            condspec = this.conditionSpectrogram(spec, paramsobj);
            switch sfmtype
                case 'madhu'
                    SFM = 2.^(-sum(condspec.*log2(condspec),1)/log2(size(condspec,1)))-1; % madhu
                case 'standard'
                    condspec = condspec + eps;
                    SFM = (prod(condspec,1).^(1/size(condspec,1)))./(sum(condspec,1)/size(condspec,1)); %standard
            end
            % Determine sliding window limits in samples
            slidewinlensamp = round(paramsobj.WindowLengthms*this.SamplingFrequency/1000);
            slidewinstepsamp = round((1 - paramsobj.OverlapFraction) * slidewinlensamp);
            slidewinendsamp = ( slidewinlensamp : slidewinstepsamp : this.NumofSamplesPerChannel );
            slidewinstartsamp = slidewinendsamp - slidewinlensamp + 1;
            % Convert sliding windows limits from samples to ms            
            slidewinstartms = slidewinstartsamp/this.SamplingFrequency*1000;
            slidewinendms = slidewinendsamp/this.SamplingFrequency*1000;
            slidewinlimssamp = [slidewinstartsamp(:) slidewinendsamp(:)];
            slidewinlimsms = [slidewinstartms(:) slidewinendms(:)];
        end % computeSlidingWindowSFM
        function [kurt, slidewinlimssamp, slidewinlimsms] = computeSlidingWindowKurtosis(this, paramsobj)
            % Computes Kurtosis of wavrecobj for sliding windows as specified by the .WindowLengthms and .OverlapFraction
            % properties of the provided FEATHISTParams object (first input). Specifying these properties equal to the
            % FEATHISTParams object used for the computation of the spectrogram will thus give 'per frame' Kurtosis.
            % Otherwise these properties can copy those specified in the 'sliding window' histogram computation of the
            % .appendSlidingWindow method of the HISTArray class to give same time windows as those of the sliding window
            % histograms.
            %
            % The first output is a vector of Kurtosis values per sliding window
            % The second output is a two column matrix of sliding window limits in samples
            % The third output is a two column matrix of sliding window limits in ms
            
            % Determine sliding window limits in samples (in the same way as is done in the .appendSlidingWindow method
            % of the HISTArray class)
            slidewinlensamp = round(paramsobj.WindowLengthms*this.SamplingFrequency/1000);
            slidewinstepsamp = round((1 - paramsobj.OverlapFraction) * slidewinlensamp);
            slidewinendsamp = ( slidewinlensamp : slidewinstepsamp : this.NumofSamplesPerChannel );
            slidewinstartsamp = slidewinendsamp - slidewinlensamp + 1;
            % Create array of overlapping sliding windows from audioread data
            wavdata = audioread(this.WavDataPath);
            wavdata = buffer(wavdata, slidewinlensamp, slidewinlensamp-slidewinstepsamp, 'nodelay');
            wavdata = wavdata(:, 1:length(slidewinendsamp));
            % Compute Kurtosis
            kurt = kurtosis(wavdata);
%             % Determine sliding window limits in samples
%             slidewinlensamp = round(paramsobj.WindowLengthms*this.SamplingFrequency/1000);
%             slidewinstepsamp = round((1 - paramsobj.OverlapFraction) * slidewinlensamp);
%             slidewinendsamp = ( slidewinlensamp : slidewinstepsamp : this.NumofSamplesPerChannel );
%             slidewinstartsamp = slidewinendsamp - slidewinlensamp + 1;
            % Convert sliding windows limits from samples to ms            
            slidewinstartms = slidewinstartsamp/this.SamplingFrequency*1000;
            slidewinendms = slidewinendsamp/this.SamplingFrequency*1000;
            slidewinlimssamp = [slidewinstartsamp(:) slidewinendsamp(:)];
            slidewinlimsms = [slidewinstartms(:) slidewinendms(:)];
        end % computeSlidingWindowKurtosis
        function [negent, slidewinlimssamp, slidewinlimsms] = computeSlidingWindowNegEntropy(this, paramsobj)
            % Computes NegEntropy of wavrecobj for sliding windows as specified by the .WindowLengthms and
            % .OverlapFraction properties of the provided FEATHISTParams object (first input). Specifying these
            % properties equal to the FEATHISTParams object used for the computation of the spectrogram will thus give
            % 'per frame' NegEntropy. Otherwise these properties can copy those specified in the 'sliding window'
            % histogram computation of the .appendSlidingWindow method of the HISTArray class to give same time windows
            % as those of the sliding window histograms.
            %
            % The first output is a vector of NegEntropy values per sliding window
            % The second output is a two column matrix of sliding window limits in samples
            % The third output is a two column matrix of sliding window limits in ms
            
            % Determine sliding window limits in samples (in the same way as is done in the .appendSlidingWindow method
            % of the HISTArray class)
            slidewinlensamp = round(paramsobj.WindowLengthms*this.SamplingFrequency/1000);
            slidewinstepsamp = round((1 - paramsobj.OverlapFraction) * slidewinlensamp);
            slidewinendsamp = ( slidewinlensamp : slidewinstepsamp : this.NumofSamplesPerChannel );
            slidewinstartsamp = slidewinendsamp - slidewinlensamp + 1;
            % Create array of overlapping sliding windows from audioread data
            wavdata = audioread(this.WavDataPath);
            wavdata = buffer(wavdata, slidewinlensamp, slidewinlensamp-slidewinstepsamp, 'nodelay');
            wavdata = wavdata(:, 1:length(slidewinendsamp));
            % Compute NegEntropy (this could be vectorised if it becomes a speed bottleneck)
            negent = nan(1,length(slidewinendsamp));
            for ii = 1:length(negent)
                slidewinmean = mean(wavdata(:,ii));
                slidewinstd = std(wavdata(:,ii));
                histlims = [-10*slidewinstd 10*slidewinstd] + slidewinmean;
                binwidth = .01;
                wavdatahist = hist(wavdata(:,ii),(histlims(1):binwidth:histlims(2)));
                wavdatahist(wavdatahist == 0) = eps;
                wavdatahist = wavdatahist/sum(wavdatahist.*(binwidth*ones(1,length(wavdatahist))));
                wavdataentr = -sum(wavdatahist.*log(wavdatahist).*(binwidth*ones(1,length(wavdatahist))));
                normcomp = slidewinstd * randn(slidewinlensamp,1) + slidewinmean;
                normhist = hist(normcomp,(histlims(1):binwidth:histlims(2)));
                normhist(normhist == 0) = eps;
                normhist = normhist/sum(normhist.*(binwidth*ones(1,length(normhist))));
                normentr = -sum(normhist.*log(normhist).*(binwidth*ones(1,length(normhist))));
                negent(ii) = normentr - wavdataentr;
            end
%             % Determine sliding window limits in samples
%             slidewinlensamp = round(paramsobj.WindowLengthms*this.SamplingFrequency/1000);
%             slidewinstepsamp = round((1 - paramsobj.OverlapFraction) * slidewinlensamp);
%             slidewinendsamp = ( slidewinlensamp : slidewinstepsamp : this.NumofSamplesPerChannel );
%             slidewinstartsamp = slidewinendsamp - slidewinlensamp + 1;
            % Convert sliding windows limits from samples to ms            
            slidewinstartms = slidewinstartsamp/this.SamplingFrequency*1000;
            slidewinendms = slidewinendsamp/this.SamplingFrequency*1000;
            slidewinlimssamp = [slidewinstartsamp(:) slidewinendsamp(:)];
            slidewinlimsms = [slidewinstartms(:) slidewinendms(:)];
        end % computeSlidingWindowNegEntropy
        function check_flag = checkAudioReadConsistency(this)
            % Compares the #samples read by audioread and those returned by .TotalSamples of audioinfo and returns 1 if
            % they are equal or 0 if not (as happens for .mp3 and other non-wav formats in Windows and Linux platforms)
            temp = audioread(this.WavDataPath);
            if size(temp,1) == this.NumofSamplesPerChannel
                check_flag = 1;
            else
                check_flag = 0;
            end
        end % checkAudioReadConsistency
    end
    methods (Static, Access = public)
        function [spec, F] = conditionSpectrogram(spec, F, paramsobj)
            %Limits spectrogram to specified freq lims, normalises frames to pmfs and sets low sum frames to NaNs
            %
            % spec = WAVRec.conditionSpectrogram(spec, paramsobj)
            %
            % First input argument is spectrogram as computed by this.computeSpectrogram
            %
            % Second input is a vector specifying the frequencies (in Hz) corresponding to the rows of the spectrogram in
            % the first input argument (this is the same as returned in the second output argument of the
            % WAVRec.computeSpectrogram method)
            %
            % Third input is a handle object of the class FEATHISTParams
            %
            % First output is an array with same #cols as first input but #rows truncated to the frequency limits
            % specified in the paramsobj.FrequencyLimitsHz property and with each column either normalised to a pmf or
            % set to NaN
            %
            % Second output is a vector of the frequencies (in Hz) kept after the truncation of frequencies outside the
            % frequency limits specified in the FEATHISTParams object input
            
            [~,tempindlow] = min(abs( F - paramsobj.FrequencyLimitsHz(1) ));
            [~,tempindhigh] = min(abs( F - paramsobj.FrequencyLimitsHz(2) ));
            spec = spec(tempindlow:tempindhigh,:);
            zeroind = (sum(spec)==0);
            if any(zeroind)
                spec(:,zeroind) = eps*abs(randn(size(spec,1),length(find(zeroind))));
            end
            F = F(tempindlow:tempindhigh);
            %% Find frames within percentile of interesting frames
            frameinterest = sum(spec);
            [~, tempind] = sort(frameinterest,'descend');
            cutoffpoint = ceil(paramsobj.InterestingFrameFraction*length(frameinterest));
            %% Normalise interesting frames and set non-interesting frames to NaN
            % Normalising factor
            temp = repmat(frameinterest(tempind(1:cutoffpoint)), size(spec,1), 1);
            % Normalise chosen frames
            spec(:,tempind(1:cutoffpoint)) = spec(:,tempind(1:cutoffpoint))./temp;
            % Set not-chosen frames to NaN
            spec(:,tempind(cutoffpoint+1:end)) = nan;
        end % conditionSpectrogram
        function histogram = computeHistogram(featurearray, paramsobj)
            %Computes histogram of given features array (with/out Laplace smoothing and normalisation to pdf)
            %
            % histogram = WAVRec.computeHistogram(featurearray, paramsobj)
            % 
            % The first input is the 6-row features array are returned by the conditionSpectrogram method.
            % The second input is a handle object of the class FEATHISTParams
            % 
            % The output is the a histogram array with as many dimensions as the non-NaN elements of the HistogramEdges
            % property of paramsobj (which is an object of the class FEATHISTParams)
            % 
            % This method uses histcn.m as provided by Bruno Luong in Mathworkd File Exchange (see
            % http://www.mathworks.co.uk/matlabcentral/fileexchange/23897-n-dimensional-histogram)
            
            %% Remove NaN entries from paramsobj.HistogramEdges cell array
            tempedges = {};
            tempfeat = [];
            for ii = 1:length(paramsobj.HistogramEdges);
                if ~isnan(paramsobj.HistogramEdges{ii});
                    tempedges{end+1} = paramsobj.HistogramEdges{ii}; %#ok<AGROW>
                    tempfeat(end+1,:) = featurearray(ii,:); %#ok<AGROW>
                end;
            end;
            %% Compute Histogram
            histogram = WAVRec.histcn(tempfeat.', tempedges{:});
            if paramsobj.HistogramSmoothing == 1
                histogram = 1 + histogram;
            elseif paramsobj.HistogramSmoothing == 0
                histogram = eps + histogram;
            else
                error('The HistogramSmoothing parameter of the provided FEATHISTParams object must be either 0 or 1');
            end
            if paramsobj.HistogramNorm == 1
                temp = histogram;
                for ii = 1:ndims(histogram)
                    temp = sum(temp);
                end
                histogram = histogram/temp;
            end

        end % computeHistogram
    end
    methods (Static, Access = private)
        function [count, edges, mid, loc] = histcn(X, varargin)
            % function [count edges mid loc] = histcn(X, edge1, edge2, ..., edgeN)
            %
            % NOTE: This is a slightly modified version of the code provided in Mathworkd File Exchange
            % (see http://www.mathworks.co.uk/matlabcentral/fileexchange/23897-n-dimensional-histogram)
            %
            % Purpose: compute n-dimensional histogram
            %
            % INPUT
            %   - X: is (M x N) array, represents M data points in R^N
            %   - edgek: are the bin vectors on dimension k, k=1...N.
            %     If it is a scalar (Nk), the bins will be the linear subdivision of
            %     the data on the range [min(X(:,k)), max(X(:,k))] into Nk
            %     sub-intervals
            %     If it's empty, a default of 32 subdivions will be used
            %
            % OUTPUT
            %   - count: n-dimensional array count of X on the bins, i.e.,
            %         count(i1,i2,...,iN) = cardinal of X such that
            %                  edge1(i1) <= X(:,i1) < edge1(i1)+1 and
            %                       ...
            %                  edgeN(iN) <= X(:,iN) < edgeN(iN)+1
            %   - edges: (1 x N) cell, each provides the effective edges used in the
            %     respective dimension
            %   - mid: (1 x N) cell, provides the mid points of the cellpatch used in
            %     the respective dimension
            %   - loc: (M x N) array, index location of X in the bins. Points have out
            %     of range coordinates will have zero at the corresponding dimension.
            %
            % DATA ACCUMULATE SYNTAX:
            %   [ ... ] = histcn(..., 'AccumData', VAL);
            %   where VAL is M x 1 array. Each VAL(k) corresponds to position X(k,:)
            %   will be accumulated in the cell containing X. The accumulate result
            %   is returned in COUNT.
            %   NOTE: Calling without 'AccumData' is similar to having VAL = ones(M,1)
            %
            %   [ ... ] = histcn(..., 'AccumData', VAL, 'FUN', FUN);
            %     applies the function FUN to each subset of elements of VAL.  FUN is
            %     a function that accepts a column vector and returns
            %     a numeric, logical, or char scalar, or a scalar cell.  A has the same class
            %     as the values returned by FUN.  FUN is @SUM by default.  Specify FUN as []
            %     for the default behavior.
            %
            % Usage examples:
            %   M = 1e5;
            %   N = 3;
            %   X = randn(M,N);
            %   [N edges mid loc] = histcn(X);
            %   imagesc(mid{1:2},N(:,:,ceil(end/2)))
            %
            % % Compute the mean on rectangular patch from scattered data
            %   DataSize = 1e5;
            %   Lat = rand(1,DataSize)*180;
            %   Lon = rand(1,DataSize)*360;
            %   Data = randn(1,DataSize);
            %   lat_edge = 0:1:180;
            %   lon_edge = 0:1:360;
            %   meanData = histcn([Lat(:) Lon(:)], lat_edge, lon_edge, 'AccumData', Data, 'Fun', @mean);
            %
            % See also: HIST, ACCUMARRAY
            %
            % Bruno Luong: <brunoluong@yahoo.com>
            % Last update: 25/August/2011
            
            if ndims(X)>2 %#ok<ISMAT>
                error('histcn: X requires to be an (M x N) array of M points in R^N');
            end
            DEFAULT_NBINS = 32;
            
            AccumData = [];
            Fun = {};
            
            % Looks for where optional parameters start
            % For now only 'AccumData' is valid
            split = find(cellfun('isclass', varargin, 'char'), 1, 'first');
            if ~isempty(split)
                for k = split:2:length(varargin)
                    if strcmpi(varargin{k},'AccumData')
                        AccumData = varargin{k+1}(:);
                    elseif strcmpi(varargin{k},'Fun')
                        Fun = varargin(k+1); % 1x1 cell
                    end
                end
                varargin = varargin(1:split-1);
            end
            
            % Get the dimension
            nd = size(X,2);
            edges = varargin;
            if nd<length(edges)
                nd = length(edges); % wasting CPU time warranty
            else
                edges(end+1:nd) = {DEFAULT_NBINS};
            end
            
            % Allocation of array loc: index location of X in the bins
            loc = zeros(size(X));
            sz = zeros(1,nd);
            % Loop in the dimension
            for d=1:nd
                ed = edges{d};
                Xd = X(:,d);
                if isempty(ed)
                    ed = DEFAULT_NBINS;
                end
                if isscalar(ed) % automatic linear subdivision
                    ed = linspace(min(Xd),max(Xd),ed+1);
                end
                edges{d} = ed;
                % Call histc on this dimension
                [~, loc(:,d)] = histc(Xd, ed, 1);
                % Use sz(d) = length(ed); to create consistent number of bins
                sz(d) = length(ed)-1;
            end % for-loop
            
            % Clean
            clear dummy
            
            % This is need for seldome points that hit the right border
            sz = max([sz; max(loc,[],1)]);
            
            % Compute the mid points
            mid = cellfun(@(e) 0.5*(e(1:end-1)+e(2:end)), edges, ...
                'UniformOutput', false);
            
            % Count for points where all coordinates are falling in a corresponding
            % bins
            if nd==1
                sz = [sz 1]; % Matlab doesn't know what is one-dimensional array!
            end
            
            hasdata = all(loc>0, 2);
            if ~isempty(AccumData)
                count = accumarray(loc(hasdata,:), AccumData(hasdata), sz, Fun{:});
            else
                count = accumarray(loc(hasdata,:), 1, sz);
            end
            
            return
            
        end % histcn
    end
end