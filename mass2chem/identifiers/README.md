
# Mapping compound identifiers 

We use an internal Azimuth identifier. 
This will compile information from HMDB, PubChem, KEGG, chem_xref and refmet.

HMDB can have multiple accession numbers. 

	  <accession>HMDB0000031</accession>
	  <status>quantified</status>
	  <secondary_accessions>
	    <accession>HMDB00031</accession>
	  </secondary_accessions>
	  <name>Androsterone</name>
	  
	  <kegg_id>C00523</kegg_id>
	  <foodb_id>FDB021881</foodb_id>
	  <drugbank_id/>
	  <chemspider_id>5668</chemspider_id>
	  <pubchem_compound_id>5879</pubchem_compound_id>
	  <pdb_id/>
	  <chebi_id>16032</chebi_id>
	  <biocyc_id/>
	  <wikipedia_id>Androsterone</wikipedia_id>
	  <knapsack_id/>
	  <phenol_explorer_compound_id/>
	  <bigg_id>35244</bigg_id>
	  <metlin_id>2797</metlin_id>
	  
As long as information is easily retrieved from a database, there's no need to replicate too much from it.

Parsing code is under mass2chem/utils/, and examples are given in /docs/.
 