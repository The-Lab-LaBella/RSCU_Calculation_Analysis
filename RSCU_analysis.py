import pandas as pd
import re
import os

def RSCU():
    path = '/scratch/bzavalam/cds_csv'
    ext = ('out.csv')

    #Read yeast csv profile of all the genomes and respective clades
    yeast_table = pd.read_csv("1154yeasts_21outrgoups_info_20220408.csv")
    clade_reassign_lst = ["CUG-Ser1 clade","CUG-Ala clade", "CUG-Ser2 clade"]

    Ser1_clade_table = yeast_table.loc[yeast_table["clade"].isin(["CUG-Ser1 clade"])]
    Ser1_lst = Ser1_clade_table["assembly_fullID"].values.tolist()

    Ala_clade_table = yeast_table.loc[yeast_table["clade"].isin(["CUG-Ala clade"])]
    Ala_lst = Ala_clade_table["assembly_fullID"].values.tolist()

    Ser2_clade_table = yeast_table.loc[yeast_table["clade"].isin(["CUG-Ser2 clade"])]
    Ser2_lst = Ser2_clade_table["assembly_fullID"].values.tolist()

    yeast_table = yeast_table[yeast_table.clade.isin(["CUG-Ser1 clade"]) == False]
    yeast_table = yeast_table[yeast_table.clade.isin(["CUG-Ser2 clade"]) == False]
    yeast_table = yeast_table[yeast_table.clade.isin(["CUG-Ala clade"]) == False]
    yeast_lst = yeast_table["assembly_fullID"].values.tolist()

    lst = []


    aa_reg_dict = {'Met':['ATG'], 'Phe':['TTT', 'TTC'], 'Leu':['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'], 'Cys':['TGT', 'TGC'], 'Tyr':['TAC', 'TAT'], 'Trp':['TGG'], 'Pro':['CCT', 'CCC', 'CCA', 'CCG'], 'His':['CAT', 'CAC'],
    'Gln':['CAA', 'CAG'], 'Arg':['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'], 'Ile':['ATT', 'ATC', 'ATA'], 'Thr':['ACT', 'ACC', 'ACA', 'ACG'],
    'Asn':['AAT', 'AAC'], 'Lys':['AAA', 'AAG'], 'Ser':['AGT', 'AGC', 'TCT', 'TCC', 'TCA', 'TCG'], 'Val':['GTT', 'GTC', 'GTA', 'GTG'],
    'Ala':['GCT', 'GCC', 'GCA', 'GCG'], 'Asp':['GAT', 'GAC'], 'Glu':['GAA', 'GAG'], 'Gly':['GGT', 'GGC', 'GGA', 'GGG'], '*':['TAA','TAG','TGA']}

    aa_Ser1_dict = {'Leu':['TTA', 'TTG', 'CTT', 'CTC', 'CTA'], 'Ser':['AGT', 'AGC', 'TCT', 'TCC', 'TCA', 'TCG', 'CTG'], 'Met':['ATG'], 'Phe':['TTT', 'TTC'], 'Cys':['TGT', 'TGC'], 'Tyr':['TAC', 'TAT'], 'Trp':['TGG'], 'Pro':['CCT', 'CCC', 'CCA', 'CCG'], 'His':['CAT', 'CAC'],
    'Gln':['CAA', 'CAG'], 'Arg':['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'], 'Ile':['ATT', 'ATC', 'ATA'], 'Thr':['ACT', 'ACC', 'ACA', 'ACG'],
    'Asn':['AAT', 'AAC'], 'Lys':['AAA', 'AAG'], 'Val':['GTT', 'GTC', 'GTA', 'GTG'],
    'Ala':['GCT', 'GCC', 'GCA', 'GCG'], 'Asp':['GAT', 'GAC'], 'Glu':['GAA', 'GAG'], 'Gly':['GGT', 'GGC', 'GGA', 'GGG'], '*':['TAA','TAG','TGA']}

    aa_Ser2_dict = {'Leu':['TTA', 'TTG', 'CTT', 'CTC', 'CTA'], 'Ser':['AGT', 'AGC', 'TCT', 'TCC', 'TCA', 'TCG', 'CTG'], 'Met':['ATG'], 'Phe':['TTT', 'TTC'], 'Cys':['TGT', 'TGC'], 'Tyr':['TAC', 'TAT'], 'Trp':['TGG'], 'Pro':['CCT', 'CCC', 'CCA', 'CCG'], 'His':['CAT', 'CAC'],
    'Gln':['CAA', 'CAG'], 'Arg':['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'], 'Ile':['ATT', 'ATC', 'ATA'], 'Thr':['ACT', 'ACC', 'ACA', 'ACG'],
    'Asn':['AAT', 'AAC'], 'Lys':['AAA', 'AAG'], 'Val':['GTT', 'GTC', 'GTA', 'GTG'],
    'Ala':['GCT', 'GCC', 'GCA', 'GCG'], 'Asp':['GAT', 'GAC'], 'Glu':['GAA', 'GAG'], 'Gly':['GGT', 'GGC', 'GGA', 'GGG'], '*':['TAA','TAG','TGA']}

    aa_ala_dict = {'Leu':['TTA', 'TTG', 'CTT', 'CTC', 'CTA'], 'Ala':['GCT', 'GCC', 'GCA', 'GCG', 'CTG'], 'Met':['ATG'], 'Phe':['TTT', 'TTC'], 'Cys':['TGT', 'TGC'], 'Tyr':['TAC', 'TAT'], 'Trp':['TGG'], 'Pro':['CCT', 'CCC', 'CCA', 'CCG'], 'His':['CAT', 'CAC'],
    'Gln':['CAA', 'CAG'], 'Arg':['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'], 'Ile':['ATT', 'ATC', 'ATA'], 'Thr':['ACT', 'ACC', 'ACA', 'ACG'],
    'Asn':['AAT', 'AAC'], 'Lys':['AAA', 'AAG'], 'Ser':['AGT', 'AGC', 'TCT', 'TCC', 'TCA', 'TCG'], 'Val':['GTT', 'GTC', 'GTA', 'GTG'], 'Asp':['GAT', 'GAC'], 'Glu':['GAA', 'GAG'], 'Gly':['GGT', 'GGC', 'GGA', 'GGG'], '*':['TAA','TAG','TGA']}

    codon_counts = aa_reg_dict.values()
    Ser1_codon_counts = aa_Ser1_dict.values()
    Ser2_codon_counts = aa_Ser2_dict.values()
    ala_codon_counts = aa_ala_dict.values()



    for files in os.listdir(path):
        if files.startswith(tuple(Ser1_lst)) or files.startswith(tuple(Ser2_lst)):
            codon_csv_Table = pd.read_csv(files)
            for codon in Ser1_codon_counts:
                if codon == ['TTA', 'TTG', 'CTT', 'CTC', 'CTA']:
                    test_codon_table = codon_csv_Table.loc[codon_csv_Table["Codon"].isin(codon)]
                    test_codon_table = test_codon_table.drop(["Amino acid","Frequency","Percentage"], axis=1)
                    average_codon_table = test_codon_table.groupby(['Sequence'], as_index=False).agg(Sum = ('Codon Count','sum'))#.sum()

                    merged = test_codon_table.merge(average_codon_table)
                    merged['RSCU'] = (merged['Codon Count'] / merged['Sum']) * len(codon)
                    merged = merged.drop(['Sum',"Codon Count"], axis=1)
                    lst.append(merged)
                elif codon == ['AGT', 'AGC', 'TCT', 'TCC', 'TCA', 'TCG', 'CTG']:
                    test_codon_table = codon_csv_Table.loc[codon_csv_Table["Codon"].isin(codon)]
                    test_codon_table = test_codon_table.drop(["Amino acid","Frequency","Percentage"], axis=1)
                    average_codon_table = test_codon_table.groupby(['Sequence'], as_index=False).agg(Sum = ('Codon Count','sum'))#.sum()

                    merged = test_codon_table.merge(average_codon_table)
                    merged['RSCU'] = (merged['Codon Count'] / merged['Sum']) * len(codon)
                    merged = merged.drop(['Sum',"Codon Count"], axis=1)
                    lst.append(merged)
                else:
                    codon_csv_Table = codon_csv_Table[codon_csv_Table.Codon.isin(['TTA', 'TTG', 'CTT', 'CTC', 'CTA']) == False]
                    codon_csv_Table = codon_csv_Table[codon_csv_Table.Codon.isin(['AGT', 'AGC', 'TCT', 'TCC', 'TCA', 'TCG', 'CTG']) == False]
                    codon_table = codon_csv_Table.loc[codon_csv_Table["Codon"].isin(codon)]
                    codon_table["RSCU"] = codon_table["Frequency"].apply(lambda x: x*len(codon))
                    codon_table = codon_table.drop(["Amino acid","Frequency","Percentage","Codon Count"], axis=1)
                    lst.append(codon_table)

            df_all = pd.concat(lst)
            df_all.to_csv(files + "_RSCU.csv", index = False, header=True)
            lst.clear()
            print(df_all)


        elif files.startswith(tuple(Ala_lst)):
            codon_csv_Table = pd.read_csv(files)
            for codon in ala_codon_counts:
                if codon == ['TTA', 'TTG', 'CTT', 'CTC', 'CTA']:
                    test_codon_table = codon_csv_Table.loc[codon_csv_Table["Codon"].isin(codon)]
                    test_codon_table = test_codon_table.drop(["Amino acid","Frequency","Percentage"], axis=1)
                    average_codon_table = test_codon_table.groupby(['Sequence'], as_index=False).agg(Sum = ('Codon Count','sum'))#.sum()

                    merged = test_codon_table.merge(average_codon_table)
                    merged['RSCU'] = (merged['Codon Count'] / merged['Sum']) * len(codon)
                    merged = merged.drop(['Sum',"Codon Count"], axis=1)
                    lst.append(merged)
                elif codon == ['GCT', 'GCC', 'GCA', 'GCG', 'CTG']:
                    test_codon_table = codon_csv_Table.loc[codon_csv_Table["Codon"].isin(codon)]
                    test_codon_table = test_codon_table.drop(["Amino acid","Frequency","Percentage"], axis=1)
                    average_codon_table = test_codon_table.groupby(['Sequence'], as_index=False).agg(Sum = ('Codon Count','sum'))#.sum()

                    merged = test_codon_table.merge(average_codon_table)
                    merged['RSCU'] = (merged['Codon Count'] / merged['Sum']) * len(codon)
                    merged = merged.drop(['Sum',"Codon Count"], axis=1)
                    lst.append(merged)
                else:
                    codon_csv_Table = codon_csv_Table[codon_csv_Table.Codon.isin(['TTA', 'TTG', 'CTT', 'CTC', 'CTA']) == False]
                    codon_csv_Table = codon_csv_Table[codon_csv_Table.Codon.isin(['GCT', 'GCC', 'GCA', 'GCG', 'CTG']) == False]
                    codon_table = codon_csv_Table.loc[codon_csv_Table["Codon"].isin(codon)]
                    codon_table["RSCU"] = codon_table["Frequency"].apply(lambda x: x*len(codon))
                    codon_table = codon_table.drop(["Amino acid","Frequency","Percentage","Codon Count"], axis=1)
                    lst.append(codon_table)

            df_all = pd.concat(lst)
            df_all.to_csv(files + "_RSCU.csv", index = False, header=True)
            lst.clear()
            print(df_all)

        elif files.startswith(tuple(yeast_lst)):
            codon_csv_Table = pd.read_csv(files)
            for codon in codon_counts:
                codon_table = codon_csv_Table.loc[codon_csv_Table["Codon"].isin(codon)]
                codon_table["RSCU"] = codon_table["Frequency"].apply(lambda x: x*len(codon))
                codon_table = codon_table.drop(["Amino acid","Frequency","Percentage","Codon Count"], axis=1)
                lst.append(codon_table)

            df_all = pd.concat(lst)
            df_all.to_csv(files + "_RSCU.csv", index = False, header=True)
            lst.clear()
            print(df_all)


    #sed -i '1i Sequence, Codon, Amino acid, Frequency, Percentage, Codon Count' filename (had to add a header to the files)


if __name__ == '__main__':
    RSCU()
