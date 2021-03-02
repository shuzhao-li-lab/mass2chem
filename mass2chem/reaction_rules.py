'''
mass differences attributed to common biochemical reactions 
(in biological systems, to distinguish from in source reactions)






'''

source_file = "edited_reaction_data.txt"
# this file was 'reverse engineered' from the MyCompoundID paper.
# next step is to calculate those from genome scale models

class RuleOfChemicalReaction:

    def __init__(self, ):

        self.rule = ''              #e.g. glucuronidation, hydrolysis
        self.mass_diff = 0
        self.moved_chem_group = ''  # e.g. NH3

        self.isotopes = []



def get_list_RuleOfChemicalReaction(source=source_file):
    '''
    To do
    
    '''
    return []