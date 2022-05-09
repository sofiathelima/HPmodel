
amino_acid_dict = [{
    'name':'Glycine',
    'single_letter':'G',
    'three_letter':'Gly',
    'Hydrophobic':'H' # technically neutral
},{
    'name':'Alanine',
    'single_letter':'A',
    'three_letter':'Ala',
    'Hydrophobic':'H'
},{
    'name':'Leucine',
    'single_letter':'L',
    'three_letter':'Leu',
    'Hydrophobic':'H'
},{
    'name':'Methionine',
    'single_letter':'M',
    'three_letter':'Met',
    'Hydrophobic':'H'
},{
    'name':'Phenylalanine',
    'single_letter':'F',
    'three_letter':'Phe',
    'Hydrophobic':'H'
},{
    'name':'Tryptophan',
    'single_letter':'W',
    'three_letter':'Trp',
    'Hydrophobic':'H'
},{
    'name':'Lysine',
    'single_letter':'K',
    'three_letter':'Lys',
    'Hydrophobic':'P'
},{
    'name':'Glutamine',
    'single_letter':'Q',
    'three_letter':'Gln',
    'Hydrophobic':'P'
},{
    'name':'Glutamic Acid',
    'single_letter':'E',
    'three_letter':'Glu',
    'Hydrophobic':'P'
},{
    'name':'Serine',
    'single_letter':'S',
    'three_letter':'Ser',
    'Hydrophobic':'H' # technically neutral
},{
    'name':'Proline',
    'single_letter':'P',
    'three_letter':'Pro',
    'Hydrophobic':'P' # technically neutral
},{
    'name':'Valine',
    'single_letter':'V',
    'three_letter':'Val',
    'Hydrophobic':'H'
},{
    'name':'Isoleucine',
    'single_letter':'I',
    'three_letter':'Ile',
    'Hydrophobic':'H'
},{
    'name':'Cysteine',
    'single_letter':'C',
    'three_letter':'Cys',
    'Hydrophobic':'H'
},{
    'name':'Tyrosine',
    'single_letter':'Y',
    'three_letter':'Tyr',
    'Hydrophobic':'P' # technically neutral
},{
    'name':'Histidine',
    'single_letter':'H',
    'three_letter':'His',
    'Hydrophobic':'P' # technically neutral
},{
    'name':'Arginine',
    'single_letter':'R',
    'three_letter':'Arg',
    'Hydrophobic':'P'
},{
    'name':'Asparagine',
    'single_letter':'N',
    'three_letter':'Asn',
    'Hydrophobic':'P'
},{
    'name':'Aspartic Acid',
    'single_letter':'D',
    'three_letter':'Asp',
    'Hydrophobic':'P'
},{
    'name':'Threonine',
    'single_letter':'T',
    'three_letter':'Thr',
    'Hydrophobic':'H'
}]

# takes amino acid dictionary, amino acid chain, and boolean telling if the chain is single or 3 letter coding
def HP_chain(amino_acid_dict, amino_acids, is_single_letter):
    HP_chain = []
    if is_single_letter == True:
        for i in amino_acids:
            for j in range(len(amino_acid_dict)):
                if amino_acid_dict[j]['single_letter'] == i:
                    HP_chain.append(amino_acid_dict[j]['Hydrophobic'])
                    break
    else:
        for i in amino_acids:
            for j in range(len(amino_acid_dict)):
                if amino_acid_dict[j]['three_letter'] == i:
                    HP_chain.append(amino_acid_dict[j]['Hydrophobic'])
                    break
    
    return HP_chain

