classdef FEATHISTParams < handle
    %FEATParams properties are the params for determining histogram from .wav & methods are for setting these params
    %
    % TODO: write validateattributes for the set methods
    
    %% FEATParams PlantUML Class Diagram
    %{
    @startuml
    skinparam classAttributeIconSize 0
    class FEATHISTParams {
        -+  AudioChannel
        -+  AudioLimits
        -+  FrequencyLimitsHz
        -+  PointsFFT96k
        -+  PointsFFT48k
        -+  PointsFFT44k
        -+  PointsFFT32k
        -+  PointsFFT22k
        -+  PointsFFT16k
        -+  WindowLengthms
        -+  WindowType
        -+  OverlapFraction
        -+  InterestingFrameFraction
        -+  FeaturesType
        -+  HistogramEdges
        -+  HistogramSmoothing
        -+  HistogramNorm
        -+  ParametersSet
        +   FEATHISTParams() : this
        +   setBriggsParams(this)
        +   setDefaultParams(this)
        +   setAudioParams(this,_,_)
        +   setFrequencyLimitsParams(this,_)
        +   setPointsFFTParams(this,_)
        +   setWindowParams(this,_,_,_)
        +   setSpectrogramConditionParams(this,_)
        +   setFeaturesType(this,_)
        +   setHistogramParams(this,_,_,_)
        +   getParams(this)
        +   deepcopy(this)
    }
    @enduml
    %}
    %% Properties
    properties (SetAccess = private, GetAccess = public)
        AudioChannel;
        % AudioChannel is a scalar taking the value 1 or 2 to indicate if the left (the default) or right channel of a
        % stereo recording is used
        AudioLimits;
        % AudioLimits is a two element vector with positive real elements taking values between (and including) 0 and 1
        % (the second being not smaller than the first) indicating which part of the audio recording is used. If this is
        % set to [0 1] (the default) then all the recording is used.
        FrequencyLimitsHz;
        % FrequencyLimitsHz is a 2-element vector with frequency limits in Hz, the second higher than the first and both
        % lower than half the sampling frequency of the wav file to be used for feature extraction and histogram
        % computation
        PointsFFT96k;
        PointsFFT48k;
        PointsFFT44k;
        PointsFFT32k;
        PointsFFT22k;
        PointsFFT16k;
        WindowLengthms;
        WindowType;
        % WindowType is a Window function gallery function handle (see www.mathworks.co.uk/help/signal/ref/window.html)
        OverlapFraction;
        % OverlapFraction is a real scalar in [0 1)
        InterestingFrameFraction;
        % InterestingFrameFraction is a positive real scalar in [0 1] specifying the fraction of all frames that are kept
        % and normalised to pdf (These amount to a number of frames equal to
        % ceil(size(specgram,2)*(interestingframefraction), the ones that rank highest in their sum(frame) value). The
        % remaining frames are set to NaN. If this second argument is set to 1 then all frames are normalised and none is
        % set to NaN.
        FeaturesType;
        % SpectraclFeaturesType is  a string equal to either 'mean_std' or 'median_IQR'.
        HistogramEdges;
        % HistogramEdges is a 6x1 cell each element of which is either the vector of edges to be used in the computation
        % of the (up to 6 dimensional) histogram for the corresponding feature or NaN to specify that the corresponding
        % feature is not included in the histogram (which will then be of less than 6 dimensions).
        % The 6 elements of the HistogramEdges cell correspond to spectral features as follows:
        % Egdes vector in cell element #1 corresponds to Mean Frequency fc as defined in Briggs09 section 2.1.2
        % Egdes vector in cell element #2 corresponds to Bandwidth BW as defined in Briggs09 section 2.1.2
        % Egdes vector in cell element #3 corresponds to D of Mean Frequency fc
        % Egdes vector in cell element #4 corresponds to D of Bandwidth BW
        % Egdes vector in cell element #5 corresponds to DD of Mean Frequency fc
        % Egdes vector in cell element #6 corresponds to DD of Bandwidth BW
        HistogramSmoothing;
        % HistogramSmoothing is a scalar equal to either 1 (Pascal smoothing to be applied, default) or 0 (no smoothing).
        HistogramNorm;
        % HistogramNorm is a scalar equal to either 1 (histogram normalised to pmf, default) or 0 (no normalisation).
        ParametersSet;
    end % properties
   %% Methods
   methods (Access = public)
       function this = FEATHISTParams % Constructor
           this.ParametersSet = 0;
       end % FEATHISTParams
       function setBriggsParams(this)
           %Sets all parameters to a virtual replica of what is described in Briggs09
           this.AudioChannel = 1;
           this.AudioLimits = [0 1];
           this.FrequencyLimitsHz = [1378 10852];
           this.PointsFFT44k = 256;
           this.PointsFFT96k = round(this.PointsFFT44k*96/44.1);
           this.PointsFFT48k = round(this.PointsFFT44k*48/44.1);
           this.PointsFFT32k = round(this.PointsFFT44k*32/44.1);
           this.PointsFFT22k = round(this.PointsFFT44k*22.05/44.1);
           this.PointsFFT16k = round(this.PointsFFT44k*16/44.1);
           this.WindowLengthms = 256/44.1;
           this.WindowType = @rectwin;
           this.OverlapFraction = 0.5;
           this.InterestingFrameFraction = 0.1;
           this.FeaturesType = {'mean';'std'};
           this.HistogramEdges = {[-Inf (100:100:9900) Inf]
                                  [-Inf (100:100:4900) Inf]
                                  NaN
                                  NaN
                                  NaN
                                  NaN};
           this.HistogramSmoothing = 1;
           this.HistogramNorm = 1;
           this.ParametersSet = 1;
       end % setBriggsParams
       function setDefaultParams(this)
           %Sets all parameters to my default
           this.AudioChannel = 1;
           this.AudioLimits = [0 1];
           this.FrequencyLimitsHz = [500 5000];
           this.PointsFFT96k = 1024;
           this.PointsFFT48k = round(this.PointsFFT96k*48/96);
           this.PointsFFT44k = round(this.PointsFFT96k*44.1/96);
           this.PointsFFT32k = round(this.PointsFFT96k*32/96);
           this.PointsFFT22k = round(this.PointsFFT96k*22.05/96);
           this.PointsFFT16k = round(this.PointsFFT96k*16/96);
           this.WindowLengthms = 5;
           this.WindowType = @rectwin;
           this.OverlapFraction = 0.5;
           this.InterestingFrameFraction = 1;
           this.FeaturesType = {'mean';'std'};
           this.HistogramEdges = {[-Inf (500:500:7500) Inf]
                                  [-Inf (500:500:3000) Inf]
                                  NaN
                                  NaN
                                  NaN
                                  NaN};
           this.HistogramSmoothing = 1;
           this.HistogramNorm = 1;
           this.ParametersSet = 1;
       end % setDefaultParams
       function setAudioParams(this,channel,timelimits)
           this.AudioChannel = channel;
           this.AudioLimits = timelimits;
       end % setAudioParams
       function setFrequencyLimitsParams(this,freqlims)
           %Sets the FrequencyLimitsHz properties (in Hz) as specified in the input.
           %
           % The input must be a 2-element vector with frequency limits in Hz, the second higher than the first and both
           % lower than half the sampling frequency of the wav file to be used for feature extraction and histogram
           % computation
           this.FrequencyLimitsHz = freqlims;
       end % setFrequencyLimitsParams
       function setPointsFFTParams(this,nfft)
           %Sets the #fftpoints properties for 96k, 48k, 44.1k, 32k, 22.05k and 16k Fs as in first to sixth elements of
           %input
           this.PointsFFT96k = nfft(1);
           this.PointsFFT48k = nfft(2);
           this.PointsFFT44k = nfft(3);
           this.PointsFFT32k = nfft(4);
           this.PointsFFT22k = nfft(5);
           this.PointsFFT16k = nfft(6);
       end % setPointsFFTParams
       function setWindowParams(this,windowlenms,windowtype,overlap)
           %Sets the WindowLengthms, WindowType and OverlapFraction properties as specified in the three inputs
           %
           % First input must be a positive real scalar (duration in ms)
           % Second input must be a function handle from the Window function gallery (see
           % http://www.mathworks.co.uk/help/signal/ref/window.html)
           % Third input must be a real scalar in [0 1)
           this.WindowLengthms = windowlenms;
           this.WindowType = windowtype;
           this.OverlapFraction = overlap;
       end % setWindowParams
       function setSpectrogramConditionParams(this,intframfrac)
           %Sets the InterestingFrameFraction property as specified in its input
           %
           % Input must be a positive real scalar in [0 1] specifying the fraction of all frames that are kept and
           % normalised to pdf (These amount to a number of frames equal to
           % ceil(size(specgram,2)*(interestingframefraction), the ones that rank highest in their sum(frame) value). The
           % remaining frames are set to NaN. If this second argument is set to 1 then all frames are normalised and none
           % is set to NaN.
           this.InterestingFrameFraction = intframfrac;
       end % setSpectrogramConditionParams
       function setFeaturesType(this,feattypes)
           % The input must be a column cell with string elements being any (non repeating) of
           % 'mean'
           % 'std'
           % 'median'
           % 'IQR'
           % 'normIQR'
           % 'MutualInformation'
           % 'SFM'
           % 'SFM_Madhu'
           % 'NegEntropy'
           % 'Kurtosis'
           
           % TODO: validateattributes for input being:
                % column cell
                % of strings
                % being any non repeating of the ones listed above
           
           this.FeaturesType = feattypes;
       end
       function setHistogramParams(this,histedges,histsmooth,histnorm)
           %Sets the HistogramEgdes, HistogramSmoothing and HistogramNorm properties equal to the three inputs
           %
           % The first input must be a 6x1 cell each element of which is either the vector of edges to be used in the
           % computation of the (up to 6 dimensional) histogram for the corresponding feature or NaN to specify that the
           % corresponding feature is not included in the histogram (which will then be of less than 6 dimensions).
           % The 6 elements of the HistogramEdges cell correspond to spectral features as follows:
           % Egdes vector in cell element #1 corresponds to Mean Frequency fc as defined in Briggs09 section 2.1.2
           % Egdes vector in cell element #2 corresponds to Bandwidth BW as defined in Briggs09 section 2.1.2
           % Egdes vector in cell element #3 corresponds to D of Mean Frequency fc
           % Egdes vector in cell element #4 corresponds to D of Bandwidth BW
           % Egdes vector in cell element #5 corresponds to DD of Mean Frequency fc
           % Egdes vector in cell element #6 corresponds to DD of Bandwidth BW
           % The second input must be a scalar equal to either 1 (Pascal smoothing to be applied, default) or 0 (no
           % smoothing).
           % The third input must be a scalar equal to either 1 (histogram normalised to pmf, default) or 0 (no
           % normalisation).
           this.HistogramEdges = histedges;
           this.HistogramSmoothing = histsmooth;
           this.HistogramNorm = histnorm;
       end % setHistogramParams
       function params = getParams(this)
           %Returns a structure copying the FEATHISTParams object properties
           s = warning('off', 'MATLAB:structOnObject');
           params = struct(this);
           warning(s.state, 'MATLAB:structOnObject');
       end % getParams
       function clone = deepcopy(this)
           if this.ParametersSet ~= 1
               error('must give FEATHISTParams object with set properties')
           end
           clone = FEATHISTParams;
           clone.AudioChannel = this.AudioChannel;
           clone.AudioLimits = this.AudioLimits;
           clone.FrequencyLimitsHz = this.FrequencyLimitsHz;
           clone.PointsFFT96k = this.PointsFFT96k;
           clone.PointsFFT48k = this.PointsFFT48k;
           clone.PointsFFT44k = this.PointsFFT44k;
           clone.PointsFFT32k = this.PointsFFT32k;
           clone.PointsFFT22k = this.PointsFFT22k;
           clone.PointsFFT16k = this.PointsFFT16k;
           clone.WindowLengthms = this.WindowLengthms;
           clone.WindowType = this.WindowType;
           clone.OverlapFraction = this.OverlapFraction;
           clone.InterestingFrameFraction = this.InterestingFrameFraction;
           clone.FeaturesType = this.FeaturesType;
           clone.HistogramEdges = this.HistogramEdges;
           clone.HistogramSmoothing = this.HistogramSmoothing;
           clone.HistogramNorm = this.HistogramNorm;
           clone.ParametersSet = this.ParametersSet;
       end % deepcopy
   end % methods
end % FEATParams





























