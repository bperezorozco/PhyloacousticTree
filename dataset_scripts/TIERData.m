classdef TIERData < handle
    % TIERData < handle
    
    %% TIERData PlantUML Class Diagram
    %{
    @startuml
    skinparam classAttributeIconSize 0
    class TIERData {
        i+  DataPath
        i+  DataTable
        i+  MetadataTable
        i+  ColPatternID
        i+  ColSpeciesLat
        i+  ColSpeciesEng
        i+  ColType
        i+  ColCallType
        i+  ColDuration
        i+  ColClassID
        i+  ColIndividualID
        i+  ClassBounds
        i+  RecSpecList
        +   TIERData(varargin) : this
        +   getClassBoundIndices(this)
        +   getMetadataTableRow(this,_)
        +   getSameIndIDs(this,_)
        +   getSameSpecies(this,_)
        +   getSpeciesNumber(this,_)
        +   findDiscrepancies(this)
        +   getMDIDsFromLatName(this,_)
        +	<<Static>> ripData()
        +   <<Static>> ripMetadata()
    }
    @enduml
    %}
    %% Properties
    properties (SetAccess = immutable, GetAccess = public)
        DataPath;
        DataTable;
        MetadataTable;
        ColPatternID;
        ColSpeciesLat;
        ColSpeciesEng;
        ColType;
        ColCallType;
        ColDuration
        ColClassID
        ColIndividualID;
        ClassBounds;
        RecSpecList;
    end
    %% Methods
    methods
        function this = TIERData(varargin)  % Constructor
            % Either give one argument (path to khf data folder) or none (ui will open)
            switch length(varargin)
                case 1
                    try
                        if varargin{1}(end) ~= filesep
                            this.DataPath = [varargin{1} filesep];
                        else
                            this.DataPath = varargin{1};
                        end
                        temp = load([this.DataPath 'data_tables.mat']);
                        this.DataTable = temp.data_table;
                        this.MetadataTable = temp.metadata_table;
                    catch
                        error(['could not find a ''data_tables.mat'' file containing ''data_table'' and '...
                            '''metadata_table'' variables in the path specified']);
                    end
                case 0
                    datapath = uigetdir([getenv('USERPROFILE') filesep ...
                        'Documents\current data\oxford biosound project'], ...
                        'Select directory with TIER data');
                    if ~datapath
                        disp('No directory specified, no object created')
                        return
                    else
                        this.DataPath = [datapath filesep];
                        try
                            temp = load([this.DataPath 'data_tables.mat']);
                            this.DataTable = temp.data_table;
                            this.MetadataTable = temp.metadata_table;
                        catch
                            error(['could not find a ''data_tables.mat'' file containing ''data_table'' and '...
                                '''metadata_table'' variables in the path specified']);
                        end
                    end
                otherwise
                    error('Please specify either a valid path with the ripped TIER data or no input argument')
            end
            this.ColPatternID = find(strcmp('Pattern ID',this.MetadataTable(1,:)));
            this.ColSpeciesLat = find(strcmp('Species (lat)',this.MetadataTable(1,:)));
            this.ColSpeciesEng = find(strcmp('Species (en)',this.MetadataTable(1,:)));
            this.ColType = find(strcmp('Type',this.MetadataTable(1,:)));
            this.ColCallType = find(strcmp('Call Type',this.MetadataTable(1,:)));
            this.ColDuration = find(strcmp('Duration [s]',this.MetadataTable(1,:)));
            this.ColClassID = find(strcmp('Class ID',this.MetadataTable(1,:)));
            this.ColIndividualID = find(strcmp('Individual ID',this.MetadataTable(1,:)));
            this.ClassBounds = this.getClassBoundIndices;
            this.RecSpecList = this.getSpeciesNumber( 1:size(this.MetadataTable,1)-1 );
        end
        function boundscell = getClassBoundIndices(this)
            % Returns cell of dims #species x 3. Column 1 is the Latin name and columns 2-3 are the first and last
            % number of recording of the corresponding species
            [C,ia,~] = unique(this.MetadataTable(2:end,this.ColSpeciesLat),'stable');
            temp = [ia(2:end);size(this.MetadataTable,1)]-1;
            boundscell = [C num2cell([ia temp])];
        end
        function metadatarow = getMetadataTableRow(this, patternid)
            %For a given (string) input equal to a recording's Pattern ID it returns the recording's ID
            % (i.e its row number in the MetadataTable - 1)
            metadatarow = find(strcmp( patternid, this.MetadataTable{2:end,this.ColPatternID} ));
        end
        function sameindids = getSameIndIDs(this, recid)
            %For a given (integer scalar) input corresponding to a recording's MetadataID (this is row number in the
            %MetadataTable - 1) it returns all recording IDs (row numbers in the MetadataTable - 1) with the same IndID
            allIndIDs = [this.MetadataTable{2:end,this.ColIndividualID}];
            recIndID = this.MetadataTable{recid+1,this.ColIndividualID};
            sameindids = find( recIndID == allIndIDs );
        end
        function specids = getSameSpecies(this, recid)
            % For a given (integer scalar) input corresponding to a MetadataID
            % (i.e. row number in the MetadataTable - 1) it returns the MetadataIDs with the same species
            boundscell = this.ClassBounds;
%             temp = find(recid >= [boundscell{:,2}] & recid <= [boundscell{:,3}]);            
            temp = this.RecSpecList(recid);
            specids = (boundscell{temp,2}:boundscell{temp,3});
        end
        function specids = getSpeciesNumber(this, recids)
            % For a given (integer vector) input corresponding to a collection of recording IDs
            % (i.e. row numbers in the MetadataTable - 1) it returns the species number of each recording
            boundscell = this.ClassBounds;
            specids = zeros(length(recids),1);
            for ii = 1:length(recids)
                temp = recids(ii) >= [boundscell{:,2}] & recids(ii) <= [boundscell{:,3}];
                specids(ii) = find(temp);
            end
        end
        function flag = findDiscrepancies(this)
            % Checks for discrepancies in downloaded data
            %% Initialise to all ok
            flag = 'All Ok';
            %% List of [dir -> filenames] checks with [metadatatable -> boundscell -> Pattern IDs]
            specdirslist = dir(this.DataPath);
            boundscell = this.ClassBounds;
            for ii = 1:length(specdirslist)
                if specdirslist(ii).isdir && ~sum(strcmp(specdirslist(ii).name,{'.','..'}))
                    perspecieswavslist = dir([this.DataPath specdirslist(ii).name]);
                    boundspecies = strcmp(specdirslist(ii).name,boundscell(:,1));
                    try
                        temp1 = {perspecieswavslist(3:end).name};
                        temp1 = cellfun(@(x) x(1:end-4), temp1, 'UniformOutput',false)';
                        temp2 = this.MetadataTable(boundscell...
                            {boundspecies,2}+1:boundscell{boundspecies,3}+1,this.ColPatternID);
                        if strcmp(temp1,temp2) ~= true(length(temp1),1)
                            error('')
                        end
                    catch
                        flag = ['discrepancy in folder ' specdirslist(ii).name];
                        return
                    end
                end
            end
        end
        function mdids = getMDIDsFromLatName(this,latname)
            %For given Latin Name (string) it returns vector of mdids (equal to row numbers in MetadataTable - 1)
            validateattributes(latname,{'char'},{'vector'})
            specnum = find(strcmp(this.ClassBounds(:,1),latname));
            if isempty(specnum)
                warning('provided Latin Name is not in the database')
                mdids = [];
            else
                mdids = (this.ClassBounds{specnum,2}:this.ClassBounds{specnum,3});
                mdids = mdids(:);
            end
        end
        function mdids = getOpenAccess(this,mdids)
            %For input mdids vector (row numbers in MetadataTable - 1) it returns (subset) vector of OpenAccess mdids
            validateattributes(mdids,{'numeric'},{'vector','integer','positive','<=',size(this.MetadataTable,1)-1})
            mdids = mdids(strcmp(this.MetadataTable(mdids + 1,strcmp(this.MetadataTable(1,:),'Open Access')),'yes'));
        end
        
    end
    methods (Static)
        function ripData
            %Downloads data (only wave data not metadata) from the Animal Sound Archive - Museum für Naturkunde Berlin
            %
            % The script creates directories per species with the .wav data unzipped in them
            %
            % The data ripping is based on finding the species names as enclosed between the expressions '&Species=' and
            % '&Page=' in the str output of the url 'http://www.animalsoundarchive.org/RefSys/Statistics.php' and on then
            % downloading the corresponding .zip files from the addresses
            % http://www.animalsoundarchive.org/RefSys/Species.php?&Species= [species name] &Page=1
            %
            % It is possible that for the function to run ok the user must first open a browser and log in to the
            % 'Register' section of the website
            %
            % For more information and links about the data see
            % http://www.animalsoundarchive.org/RefSys/ProjectDescription.php The password is 'bubobubo'
            
            %% Specify the wave data weblink and get the save path from user
            stats_url = 'http://www.animalsoundarchive.org/RefSys/Statistics.php';
            download_url = 'http://www.animalsoundarchive.org/RefSys/downloadSpecies.php?&Species=';
            save_path = uigetdir(cd, 'Select directory to save the wave data files');
            if ~save_path
                disp('No directory specified, no data saved')
                return
            else
                save_path = [save_path '\'];
            end
            %% Download data
            urlstr = urlread(stats_url);
            indstart = strfind(urlstr,'&Species=');
            indend = strfind(urlstr,'&Page=');
            species_names = cell(length(indstart),1);
            for ii = 1:length(species_names)
                species_names{ii} = urlstr(indstart(ii)+9:indend(ii)-1);
                urlwrite([download_url strrep(species_names{ii}, ' ', '%20')],[save_path species_names{ii} '.zip']);
                mkdir([save_path species_names{ii}]);
                unzip([save_path species_names{ii} '.zip'],[save_path species_names{ii}]);
                delete([save_path species_names{ii} '.zip'])
            end
        end
        function ripMetadata
            % Guides the user through a series of copy paste actions during which metadata from the Animal Sound Archive
            % of the Museum für Naturkunde Berlin are downloaded and saved in 2 xlsx files, 'data_table.xlsx' and
            % 'metadata_table.xlsx'.
            
            %% Web addresses and filenames
            species_url = 'http://www.animalsoundarchive.org/RefSys/Statistics.php';
            overall_url = 'http://www.animalsoundarchive.org/RefSys/Preview.php';
            password = 'bubobubo';
            species_table_xlsx_filename = 'data_table.xlsx';
            metadata_table_xlsx_filename = 'metadata_table.xlsx';
            %% Get path from user
            save_path = uigetdir(cd, 'Select directory to save the metadata files');
            if ~save_path
                disp('No directory specified, no metadata saved')
            return
            else
                save_path = [save_path filesep];
            end
            %% Rip 'Statistics' table for per-species downloadable data
            xlswrite([save_path species_table_xlsx_filename],NaN);
            % Copy paste the species table from the website to the xlxs file
            web(species_url,'-browser');
            winopen([save_path species_table_xlsx_filename])
            disp('copy-paste the table from the webpage that just opened in Matlab''s Web Browser')
            disp('(including the header and sum lines)')
            disp(['to the ' species_table_xlsx_filename ' file that just opened'])
            disp('and save/close the xlsx file w/o changing its name.')
            disp(['if needed the website can be accessed through ' overall_url])
            disp(['and the password is ' password])
            disp('when done execute ''return'' in the prompt to exit keyboard mode')
            keyboard            
            %% Download and convert 'Metadata' table for per-recording metadata
            % Download .csv file
            disp('Go to http://www.animalsoundarchive.org/RefSys/MetadataToCsv1FormPage.php')
            disp('Or to http://www.animalsoundarchive.org/RefSys/Metadata.php ''CSV Export'' link')
            disp('Set options to ''All'', ''All'', ''as selected'', ''Semicolon'', ''None'', ''Windows''')
            disp(['and download .csv file in the folder ' save_path ' with its original name'])
            disp('when done execute ''return'' in the prompt to exit keyboard mode')
            keyboard
            % Convert to .xlsx
            disp(['Create a ' save_path metadata_table_xlsx_filename ' excel file'])
            disp('and in it open the just downloaded .csv file with delimiter set to '';''')
            disp(['Then save it in ' save_path])
            disp(['as ' metadata_table_xlsx_filename])
            disp('when done execute ''return'' in the prompt to exit keyboard mode')
            keyboard
            %% Save data_table and metadata_table variables in data_tables.mat
            % Unwrap interleaved first column of the 'Statistics' table as copy pasted in data_table.xlsx 
            [~, ~, temp] = xlsread([save_path 'data_table.xlsx']);
            data_table(1,:) = {temp{1,1} 'Common Name' temp{1,2:end}};
            data_table{(size(temp,1))/2+1,end} = NaN;
            for ii = 1:size(temp,1)/2-1
                data_table{ii+1,1} = temp{ii*2,1};
                data_table{ii+1,2} = temp{ii*2+1,1};
                data_table(ii+1,3:end) = temp(ii*2,2:end);
            end
            data_table(end,1:end-2) = temp(end,[1 1 2:end-2]);
            data_table{end,end-1} = NaN; %#ok<NASGU>
            % Get the table downloaded in metadata_table.xlsx
            [~, ~, metadata_table] = xlsread([save_path 'metadata_table.xlsx']); %#ok<NASGU>
            save([save_path 'data_tables'], 'data_table', 'metadata_table');
            
        end
    end
end