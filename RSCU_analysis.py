import pandas as pd
#import os

def RSCU():
    #path = 'full/path/to/directory'
    #ext = ('codonTable.out')
    #codonTable = pd.read_csv("kluyveromyces_marxianus.final.cds.all_codonTable.out")
    codon_csv_Table = pd.read_csv("kluyveromyces_marxianus.final.cds.all_codonTable.out")
    yeast_table = pd.read_csv("1154yeasts_21outrgoups_info_20220408.csv")

    #sed -i '1i Sequence, Codon, Amino acid, Frequency, Percentage, Codon Count' filename (had to add a header to the files)


    clade = yeast_table["clade"]== "CUG-Ser1 clade"

    aa_dict = {'Met':['ATG'], 'Phe':['TTT', 'TTC'], 'Leu':['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'], 'Cys':['TGT', 'TGC'], 'Tyr':['TAC', 'TAT'], 'Trp':['TGG'], 'Pro':['CCT', 'CCC', 'CCA', 'CCG'], 'His':['CAT', 'CAC'],
    'Gln':['CAA', 'CAG'], 'Arg':['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'], 'Ile':['ATT', 'ATC', 'ATA'], 'Thr':['ACT', 'ACC', 'ACA', 'ACG'],
    'Asn':['AAT', 'AAC'], 'Lys':['AAA', 'AAG'], 'Ser':['AGT', 'AGC', 'TCT', 'TCC', 'TCA', 'TCG'], 'Val':['GTT', 'GTC', 'GTA', 'GTG'],
    'Ala':['GCT', 'GCC', 'GCA', 'GCG'], 'Asp':['GAT', 'GAC'], 'Glu':['GAA', 'GAG'], 'Gly':['GGT', 'GGC', 'GGA', 'GGG'], '*':['TAA','TAG','TGA']}

    codon_counts = aa_dict.values()
    #print(codon_counts)



    lst = []

    for codon in codon_counts:
        he_table = codon_csv_Table.loc[codon_csv_Table["Codon"].isin(codon)]
        he_table["RSCU"] = he_table["Frequency"].apply(lambda x: x*len(codon))
        lst.append(he_table)

    df_all = pd.concat(lst)
    df_all = df_all.drop(["Amino acid","Frequency","Percentage","Codon Count"], axis=1)
    print(df_all.iloc[10500])


        # if i == ['ATG']:
        #     print("p")
        #     #he_table = codon_csv_Table.loc[codon_csv_Table["Codon"].isin(i)]
        # elif i == ['TTT', 'TTC']:
        #     print("y")
        #
        # elif i == ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG']:
        #     print("o")





if __name__ == '__main__':
    RSCU()
