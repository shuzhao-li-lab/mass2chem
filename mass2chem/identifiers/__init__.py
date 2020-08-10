"""

JSON data in the Metabolic Atlas repo 
https://github.com/SysBioChalmers/Human-GEM
retrieved on 2020-03-31 


# humanGEMMetAssoc.JSON  humanGEMRxnAssoc.JSON

​import json

jm = json.load(open('Human-GEM/data/annotation/humanGEMMetAssoc.JSON'))

​print(jm.keys())

dict_keys(['mets', 'metsNoComp', 'metBiGGID', 'metKEGGID', 'metHMDBID', 'metChEBIID', 'metPubChemID', 'metLipidMapsID', 'metEHMNID', 'metHepatoNET1ID', 'metRecon3DID', 'metMetaNetXID'])



RefMet is from Metabolomics Workbench
Retrived 2020-03-29 


"""

