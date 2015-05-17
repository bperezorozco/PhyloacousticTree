classdef TIERRec < WAVRec
    % TIERRec < WAVRec
    
    %% TIERRec PlantUML Class Diagram
    %{
    @startuml
    skinparam classAttributeIconSize 0
    class TIERRec {
        -+  MetadataID
        -+  RecID
        -+  LatinName
        -+  CommonName
        -+  IndID
        +   TIERRec(_,TD)
    }
    @enduml
    %}
    %% Properties
    properties (SetAccess = private, GetAccess = public)
        MetadataID; %This is equal to the row number in the TIERData.MetadataTable cell excluding the header first row
                    %i.e. for the recording described in TIERData.MetadataTable(x,:) it will be equal to x-1;
        RecID;
        LatinName;
        CommonName;
        IndID;
    end
    %% Methods
    methods
        function this = TIERRec(metadataid, tierdataobj) %Constructor
            templatinname = tierdataobj.MetadataTable{metadataid+1,tierdataobj.ColSpeciesLat};
            temprecid = tierdataobj.MetadataTable{metadataid+1,tierdataobj.ColPatternID};
            this = this@WAVRec([tierdataobj.DataPath templatinname filesep temprecid '.wav']);
            this.MetadataID = metadataid;
            this.RecID = tierdataobj.MetadataTable{this.MetadataID+1,tierdataobj.ColPatternID};
            this.LatinName = tierdataobj.MetadataTable{this.MetadataID+1,tierdataobj.ColSpeciesLat};
            this.CommonName = tierdataobj.MetadataTable{this.MetadataID+1,tierdataobj.ColSpeciesEng};
            this.IndID = tierdataobj.MetadataTable{this.MetadataID+1,tierdataobj.ColIndividualID};
        end
    end
end