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

    codon_dict = {'1':['ATG'], '2':['TTT', 'TTC'], '6':['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'], '2':['TGT', 'TGC'], '2':['TAC', 'TAT'], '1':['TGG'], '4':['CCT', 'CCC', 'CCA', 'CCG'], '2':['CAT', 'CAC'],
    '2':['CAA', 'CAG'], '6':['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'], '3':['ATT', 'ATC', 'ATA'], '4':['ACT', 'ACC', 'ACA', 'ACG'],
    '2':['AAT', 'AAC'], '2':['AAA', 'AAG'], '6':['AGT', 'AGC', 'TCT', 'TCC', 'TCA', 'TCG'], '4':['GTT', 'GTC', 'GTA', 'GTG'],
    '4':['GCT', 'GCC', 'GCA', 'GCG'], '2':['GAT', 'GAC'], '2':['GAA', 'GAG'], '4':['GGT', 'GGC', 'GGA', 'GGG'], '3':['TAA','TAG','TGA']}

    codon_counts = codon_dict.keys()
    #print(codon_counts)


    # Phe_table = codon_csv_Table.loc[codon_csv_Table["Codon"].isin(["TTT","TTC"])]
    # Phe_table["Frequency"] = Phe_table["Frequency"].apply(lambda x: x*2)
    # print(Phe_table)

    Phe_table = codon_csv_Table.loc[codon_csv_Table["Codon"].isin(["TTT","TTC"])]
    Phe_table["RSCU"] = codon_csv_Table["Frequency"].apply(lambda x: x*2)
    #print(Phe_table)

    Leu_table = codon_csv_Table.loc[codon_csv_Table["Codon"].isin(['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'])]
    Leu_table["RSCU"] = codon_csv_Table["Frequency"].apply(lambda x: x*6)
    #print(Leu_table)

    Cys_table = codon_csv_Table.loc[codon_csv_Table["Codon"].isin(['TGT', 'TGC'])]
    Cys_table["RSCU"] = codon_csv_Table["Frequency"].apply(lambda x: x*2)

    Tyr_table = codon_csv_Table.loc[codon_csv_Table["Codon"].isin(['TAC', 'TAT'])]
    Tyr_table["RSCU"] = codon_csv_Table["Frequency"].apply(lambda x: x*2)

    df_all = pd.concat([Phe_table, Leu_table, Cys_table, Tyr_table])
    print(df_all)
    #print(Phe_table)
    #Leu_table







if __name__ == '__main__':
    RSCU()
