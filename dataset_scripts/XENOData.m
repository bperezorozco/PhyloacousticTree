classdef XENOData < handle
    % XENOData < handle
    
    %% XENOData PlantUML Class Diagram
    %{
    @startuml
    skinparam classAttributeIconSize 0
    class XENOData {
        i+  DataPath
        i+  DataTable
        i+  MetadataTable
        i+  ColCatalogueNumber
        i+  ColSpeciesLat
        i+  ColSpeciesEng
        i+  ColSongType (removed)
        i+  ColBackLatin (removed)
        i+  SpeciesStruct (removed)
        i+  RecSpecList (removed)
        +   XENOData(varargin) : this
        +   convMDIDstoSpeciesIDs(this,_)
        +   convCatNumstoMDIDs(this,_)
        +   convMDIDstoCatNums(this,_)
        +   getPerSpeciesMetadataIDs(this,_)
        +   convMDIDstoLatinNames(this,_)
        +   convMDIDstoEnglishNames(this,_)
        +   getPathfromMDIDs(this,_)
    }
    @enduml
    %}
    %% Properties
    properties (SetAccess = immutable, GetAccess = public)
        DataPath;
        DataTable;
        MetadataTable;
        ColMetadataID;
        ColCatalogueNumber;
        ColSpeciesLat;
        ColSpeciesEng;
%         ColSongType;
%         ColBackLatin;
%         SpeciesStruct;
%         RecSpecList;
    end
    %% Methods
    methods
        function this = XENOData(xcdatapath)  % Constructor
            % xcdatapath must be path to downloaded xeno-canto data folder
            try
                if xcdatapath(end) ~= filesep
                    this.DataPath = [xcdatapath filesep];
                else
                    this.DataPath = xcdatapath;
                end
                temp = load([this.DataPath 'Tables.mat']);
                this.DataTable = temp.DataTableCell;
                this.MetadataTable = temp.MetadataTableCell;
                this.ColMetadataID = find(strcmpi('MetadataID',this.MetadataTable(1,:)));
                this.ColSpeciesLat = [find(strcmpi('gen',this.MetadataTable(1,:))) ...
                    find(strcmpi('sp',this.MetadataTable(1,:)))];
                this.ColCatalogueNumber = find(strcmpi('id',this.MetadataTable(1,:)));
                this.ColSpeciesEng = find(strcmpi('en',this.MetadataTable(1,:)));
            catch
                error('Tables.mat containing DataTableCell and MetadataTableCell vars not in path specified');
            end
            % Convert Catalogue Numbers from string (as returned by XC api) to nums
            temp = this.MetadataTable(2:end,this.ColCatalogueNumber);
            temp = cellfun(@str2num,temp);
            this.MetadataTable(2:end,this.ColCatalogueNumber) = num2cell(temp);
            clear temp
        end % XENOData
        function specids = convMDIDstoSpeciesIDs(this, mdids)
            % For a given (integer vector) input corresponding to a collection of metadata IDs (as listed in
            % [XENOData.MetadataTable{2:end,XENOData.ColMetadataID}]) it returns the species number of each recording (as
            % listed in [XENOData.DataTable{2:end,1}])
            specids = nan(size(mdids));
            for ii = 1:length(mdids)
                temp = [this.MetadataTable{mdids(ii)+1,this.ColSpeciesLat(1)} ' ' ...
                    this.MetadataTable{mdids(ii)+1,this.ColSpeciesLat(2)}];
                specids(ii) = this.DataTable{strcmpi(temp,this.DataTable(:,2)),1};
            end
        end % convMDIDstoSpeciesIDs
        function mdids = convCatNumstoMDIDs(this,catnums)
            % For an input vector containing catalogue numbers of XC recordings it returns the corresponding metadata
            % IDs (as listed in [XENOData.MetadataTable{2:end,XENOData.ColMetadataID}]).
            % Input must be a nonempty vector of positive integers.
            % For any element in the input vector that does not correspond to a XC recording catalogue numbers in
            % the dataset, a NaN is returned
            % TODO: could be vectorised with arrayfun for better speed?
            
            validateattributes(catnums,{'numeric'},{'nonempty','vector','integer','positive'});
            catnumsall = [this.MetadataTable{2:end,this.ColCatalogueNumber}];
            mdids = nan(size(catnums));
            for ii = 1:length(catnums)
                temp = find(catnums(ii) == catnumsall);
                if any(temp)
                    mdids(ii) = this.MetadataTable{temp+1,this.ColMetadataID};
                end
            end
        end % convCatNumstoMDIDs
        function catnums = convMDIDstoCatNums(this,mdids)
            % For an input vector containing metadata IDs (as listed in
            % [XENOData.MetadataTable{2:end,XENOData.ColMetadataID}]) it returns the corresponding catalogue numbers of
            % XC recordings.
            % Input must be a nonempty vector of positive integers.
            % For any element in the input vector that does not correspond to a metadataid a NaN is returned
            % TODO: can be vectorised?
            
            validateattributes(mdids,{'numeric'},{'nonempty','vector','integer','positive','nonnan'});
            mdidsall = [this.MetadataTable{2:end,this.ColMetadataID}];
            mdidsall = [nan ; mdidsall(:)];
            catnums = nan(size(mdids));
            for ii = 1:length(mdids)
                if ismember(mdids(ii),mdidsall)
                    catnums(ii) = this.MetadataTable{mdids(ii) == mdidsall ,this.ColCatalogueNumber};
                end
            end
        end % convMDIDstoCatNums
        function mdids = getPerSpeciesMetadataIDs(this,specid)
            % For a sclar input corresponding to a Species ID (as listed in [XENOData.DataTable{2:end,1}]) it returns
            % the corresponding metadata IDs (as listed in [XENOData.MetadataTable{2:end,XENOData.ColMetadataID}]).
            % Input must be a nonempty positive integer scalar.
            % If the input does not correspond to a Species ID an emptry array is returned
            
            validateattributes(specid,{'numeric'},{'nonempty','scalar','integer','positive','nonnan'});
            specidsall = this.DataTable(:,1);
            specidsall{1} = nan;
            specidsall = [specidsall{:}];
            try
                mdids = this.DataTable{ specid == specidsall ,  3 };
            catch
                mdids = [];
            end
        end % getPerSpeciesMetadataIDs
        function latnamescell = convMDIDstoLatinNames(this,mdids)            
            % For an input vector containing metadata IDs (as listed in
            % [XENOData.MetadataTable{2:end,XENOData.ColMetadataID}]) it returns a cell containing the corresponding
            % Latin names.
            % Input must be a nonempty vector of positive integers.
            % If any element in the input vector does not correspond to a metadataid an error is tripped
            % TODO: can be vectorised?
            
            validateattributes(mdids,{'numeric'},{'nonempty','vector','integer','positive','nonnan'});
            mdidsall = [this.MetadataTable{2:end,this.ColMetadataID}];
            mdidsall = [nan ; mdidsall(:)];
            latnamescell = cell(size(mdids));
            try
                for ii = 1:length(mdids)
                    tempind = mdids(ii) == mdidsall;
                    latnamescell{ii} = [this.MetadataTable{tempind,this.ColSpeciesLat(1)} ' ' ...
                        this.MetadataTable{tempind,this.ColSpeciesLat(2)}];
                end
            catch
                error('input not corresponding to mdid')
            end
        end % convMDIDstoLatinNames
        function engnamescell = convMDIDstoEnglishNames(this,mdids)            
            % For an input vector containing metadata IDs (as listed in
            % [XENOData.MetadataTable{2:end,XENOData.ColMetadataID}]) it returns a cell containing the corresponding
            % English names.
            % Input must be a nonempty vector of positive integers.
            % If any element in the input vector does not correspond to a metadataid an error is tripped
            % TODO: can be vectorised?
            
            validateattributes(mdids,{'numeric'},{'nonempty','vector','integer','positive','nonnan'});
            mdidsall = [this.MetadataTable{2:end,this.ColMetadataID}];
            mdidsall = [nan ; mdidsall(:)];
            engnamescell = cell(size(mdids));
            try
                for ii = 1:length(mdids)
                    tempind = mdids(ii) == mdidsall;
                    engnamescell{ii} = this.MetadataTable{tempind,this.ColSpeciesEng(1)};
                end
            catch
                error('input not corresponding to mdid')
            end
        end % convMDIDstoEnglishNames
        function recpathscell  = getPathfromMDIDs(this,mdids)
            % For an input vector containing metadata IDs (as listed in
            % [XENOData.MetadataTable{2:end,XENOData.ColMetadataID}]) it returns a cell containing the corresponding
            % recordings datapath and filenames.
            % Input must be a nonempty vector of positive integers.
            % If any element in the input vector does not correspond to a metadataid an error is tripped
            % TODO: can be vectorised?
            
            validateattributes(mdids,{'numeric'},{'nonempty','vector','integer','positive','nonnan'});
            mdidsall = [this.MetadataTable{2:end,this.ColMetadataID}];
            mdidsall = [nan ; mdidsall(:)];
            recpathscell = cell(size(mdids));
            try
                for ii = 1:length(mdids)
                    tempind = mdids(ii) == mdidsall;
                    tempfname = [num2str(this.MetadataTable{tempind,this.ColCatalogueNumber}) '.xcrec'];
                    tempdir = [this.MetadataTable{tempind,this.ColSpeciesLat(1)} ' ' ...
                        this.MetadataTable{tempind,this.ColSpeciesLat(2)}];
                    recpathscell{ii} = [this.DataPath tempdir filesep tempfname];
                end
            catch
                error('input not corresponding to mdid')
            end
        end % getPathfromMDIDs
%         function specid = getSpeciesIDforLatName(this,speclatname)
%             
%         end % getSpeciesIDforLatName
    end
    methods (Static)
        function [DataTableCell, MetadataTableCell] = ripData(speccell, qualstr, backgrstr)
            % speccell has to be a cell with strings containing latin names (i.e. 'Cyanistes caeruleus')
            % qualstr has to be of the form 'q:A' or 'q<:C' or 'q>:C'
            % backgrstr has to be either 'none' or 'all' (in the first case only recs without other species in the
            % background are downloaded
            
            %% Specify the wave data weblink and get the save path from user
            save_path = uigetdir(cd, 'Select directory to save the wave data files');
            if ~save_path
                disp('No directory specified, no data saved')
                return
            else
                save_path = [save_path filesep];
            end
            
            %% Download data
            speccellhttp = cellfun(@(x)strrep(x,' ','%20'),speccell,'UniformOutput',0);
            allrecscell = cell(length(speccell),1);
            for jj = 1:length(speccell)
                mkdir([save_path speccell{jj}]);
                urlwrite(['http://www.xeno-canto.org/api/2/recordings?query=' speccellhttp{jj} '%20' qualstr],...
                    'temp.json');
                datori = loadjson('temp.json');
                allrecscell{jj} = cell(length(datori.recordings),1);
                switch backgrstr
                    case 'none'
                        kk = 1;
                        for ii = 1:length(datori.recordings)
                            teststr = urlread(datori.recordings{ii}.url);
                            if ~isempty(strfind(teststr,'<tr><td>Background</td><td valign=''top''>none</td></tr>'))
                                disp(['downloading rec #' num2str(kk)])
                                temppathfname = [save_path speccell{jj} filesep datori.recordings{ii}.id '.xcrec'];
                                urlwrite(datori.recordings{ii}.file,temppathfname);
                                try
                                    tempaudioinfo = audioinfo(temppathfname); %#ok<NASGU>
                                    allrecscell{jj}(kk) = datori.recordings(ii);
                                    kk = kk + 1;
                                catch
                                    warning(['audioinfo problem in ' temppathfname])
                                    pause(10)
                                    warning(['deleting ' temppathfname])
                                    delete(temppathfname)
                                    if exist(temppathfname,'file')
                                        error('could not delete file')
                                    end
                                end
                            end
                        end
                        allrecscell{jj} = allrecscell{jj}(1:kk-1);
                    case 'all'
                        allrecscell{jj} = datori.recordings;
                        for ii = 1:length(datori.recordings)
                            disp(['downloading rec #' num2str(ii)])
                            urlwrite(datori.recordings{ii}.file,...
                                [save_path speccell{jj} filesep datori.recordings{ii}.id '.xcrec']);
                        end
                    otherwise
                        error('third input must be either ''none'' or ''all''')
                end
                delete('temp.json')
            end
            
            %% Create DataTable and MetadataTable
            DataTableCell{1,1} = 'Species ID';
            DataTableCell(2:length(speccell)+1,1) = num2cell((1:length(speccell))');
            DataTableCell{1,2} = 'Latin Name';
            DataTableCell(2:end,2) = speccell(:);
            DataTableCell{1,3} = 'Rec MDIDs';
            for ii = 1:length(speccell)
                DataTableCell{ii+1,3} = nan(length(allrecscell{ii}),1);
                if isempty(allrecscell{ii})
                    warning(['No records where found for species #' num2str(ii)])
                end
            end
            temp1 = cellfun(@length,DataTableCell(2:end,3));
            temp = find(temp1 ~= 0, 1, 'first');
            if ~isempty(temp)
                MetadataTableCell(1,:) = [{'MetadataID'} fieldnames(allrecscell{temp}{1})'];
                MetadataTableCell{sum(temp1),end} = 'Problem with MetadataTable creation';
                kk = 1;
                for ii = 1:length(allrecscell)
                    for jj = 1:length(allrecscell{ii})
                        MetadataTableCell(kk+1,:) = [kk struct2cell(allrecscell{ii}{jj})'];
                        DataTableCell{ii+1,3}(jj) = kk;
                        kk = kk + 1;
                    end
                end
                save([save_path 'Tables.mat'],'DataTableCell','MetadataTableCell');
                try
                    xlswrite('DataTable,xlsx',DataTableCell);
                    xlswrite('MetadataTable.xlsx',MetadataTableCell);
                catch
                    warning('No xlsx files written, probably because of xlswrite problem in OSX')
                end
            else
                warning('No records where found for any of the specified species')
                DataTableCell = {};
                MetadataTableCell = {};
            end
            
        end % ripData
    end
end