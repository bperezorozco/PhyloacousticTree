classdef XENORec < WAVRec
    % XENORec < WAVRec
    
    %% XENORec PlantUML Class Diagram
    %{
    @startuml
    skinparam classAttributeIconSize 0
    class XENORec {
        -+  MetadataID
        -+  CatalogueNumber
        -+  LatinName
        -+  CommonName
        -+  SpeciesID
    }
    @enduml
    %}
    %% Properties
    properties (SetAccess = private, GetAccess = public)
        MetadataID;
        CatalogueNumber;
        LatinName;
        CommonName;
        SpeciesID;
    end
    %% Methods
    methods
        function this = XENORec(metadataid, xenodataobj) %Constructor
            tempcell = xenodataobj.getPathfromMDIDs(metadataid);
            this  = this@WAVRec(tempcell{1});
            clear tempcell;
            this.MetadataID = metadataid;
            this.CatalogueNumber = xenodataobj.convMDIDstoCatNums(metadataid);
            this.SpeciesID = xenodataobj.convMDIDstoSpeciesIDs(metadataid);
            tempcell = xenodataobj.convMDIDstoLatinNames(metadataid);
            this.LatinName = tempcell{1};
            clear tempcell;
            tempcell = xenodataobj.convMDIDstoEnglishNames(metadataid);
            this.CommonName = tempcell{1};
            clear tempcell;
        end
    end
end