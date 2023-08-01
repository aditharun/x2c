import os
import pandas as pd 
import numpy as np
#from IPython.display import clear_output
from multiprocess import Pool, cpu_count

setlen=lambda x:len(set(x)) # Calculate length of set of a list.

datadir = "../data/genomic-data"
outputdir = "../data/braf-data"

df_mutfiles=pd.read_csv(os.path.join(datadir, 'df_mutfiles_specifc.7z'),sep='\t',dtype=str)

df_mutfiles_cystine=pd.read_csv(os.path.join(datadir, 'df_mutfiles_cystine_specific.7z'),sep='\t',dtype=str)

geneset=list(set(df_mutfiles.Hugo_Symbol.values))

histcodes=list(set(df_mutfiles.CODE.sort_values()))

geneset_cys=list(set(df_mutfiles_cystine.Hugo_Symbol.values))

muts_to_keep = ['p.Val600Glu', 'p.Val600Lys']

df = df_mutfiles[df_mutfiles['HGVSp'].isin(muts_to_keep)]

geneset_v600_list = list(set(df["Hugo_Symbol"] + df["HGVSp"]))

df_Out_v600=pd.DataFrame(0,columns=['All']+histcodes,index=geneset_v600_list+['Total'],dtype=int)


def fun_Outdf(inpvec):
    df_Out1=inpvec[1]
    inpdf=inpvec[0]
    for idx,row in inpdf.iterrows():
        igene=row.Hugo_Symbol
        imutsite=row.HGVSp
        icode=row.CODE
        """
        print("Values:")
        print("igene: "+str(igene))
        print("imut: "+str(imutsite))
        print("icode: "+str(icode))
        print([igene+'_'+imutsite,icode])
        print(df_Out1)
        """
        df_Out1.loc[igene+imutsite,icode]=df_Out1.loc[igene+imutsite,icode]+1
        df_Out1.loc[igene+imutsite,'All']=df_Out1.loc[igene+imutsite,'All']+1
    return df_Out1


df_Out_v600=fun_Outdf([df, df_Out_v600])

df['StidSid']=[df.loc[idx,'STUDY_ID']+'_'+df.loc[idx,'Tumor_Sample_Barcode'] for idx in df.index]

# Count number of cases 'Total' within each histology : Perform this action after filtering for curated and sequenced samples

TotalRow=[setlen(df.StidSid.values)]+[setlen(df.StidSid[df.CODE==ihist].values) for ihist in df_Out_v600.columns[1:]]
df_Out_v600.loc['Total']=TotalRow

for item in ['n/a','NAN','NA','na','nan']:
    if item in df_Out_v600.index:
        df_Out_v600=df_Out_v600.drop(index=item)

#remove any histologies for which there are 0 mutations observed
df_Out_v600=df_Out_v600.drop(columns=df_Out_v600.columns[df_Out_v600.loc['Total']==0])

#df_chDegen=pd.read_csv(os.path.join(datadir, 'Genelist_ManyChromosomes.xlsx'))# Only MARCH1 MARCH2 and SEPT15 have any real issues. This is very very likely due to excel errors someone made in the past by copying data incorrectly without realizing.
#degenchlist=df_chDegen[['Hugo_Symbol','Chromosome Locations','Entries_per_Location']].applymap(lambda x:str(x).upper()).values


# Identifying chromosomal degeneracies takes half an hour using 6 cores at >2Ghz. If not updating raw data, uncomment this step to skip next block.
def fun_chromosome(inplist):
    setlen=lambda x:len(set(x))
    df_mutfiles=inplist[0]
    geneset=inplist[1]
    return [igene for igene in geneset if setlen(df_mutfiles[df_mutfiles.Hugo_Symbol==igene].Chromosome)>1]

ncores=cpu_count()-2# number of cores the task can be split into
imarkers=[i*int(len(geneset)/ncores) for i in range(ncores+1)]
imarkers[-1]=len(geneset)
gslist=[[df_mutfiles,geneset[imarkers[i]:imarkers[i+1]]] for i in range(ncores)]
if __name__ == '__main__': #I don't quite understand why this is necessary but it is a part of multiprocessing docs
    po=Pool(ncores) # invoke 6 pooled threads/processes. 
    list_degen=list(po.map(fun_chromosome,gslist)) 
    po.close() 
    po.join()
degenchlist=[elem for row in list_degen for elem in row]
degenchlist=[[igene,set(df_mutfiles[df_mutfiles.Hugo_Symbol==igene].Chromosome)] for igene in degenchlist]
degenchlist=[[igene,ichlist,[sum(df_mutfiles[df_mutfiles.Hugo_Symbol==igene].Chromosome==ich) for ich in ichlist]] for igene,ichlist in degenchlist]
# Check for chromosome assignment issues in the genomic dataset. > THey are present but minimal.
df_chDegen=pd.DataFrame(degenchlist, columns=['Hugo_Symbol','Chromosome Locations','Entries_per_Location'])
df_chDegen.to_csv(os.path.join(outputdir, 'Genelist_ManyChromosomes.xlsx'))# Only MARCH1 MARCH2 and SEPT15 have any real issues. This is very very likely due to excel errors someone made in the past by copying data incorrectly without realizing.

#Import alternate gene nomenclature file from cbioportal: https://docs.cbioportal.org/3.-cbioportal-maintenance/updating-gene-and-gene_alias-tables
#Homo_sapien.gene_info.gz ftp://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz
dfGeneNames=pd.read_csv(os.path.join(datadir, 'Homo_sapiens_gene_info_GM.txt'),sep='\t',dtype=str)
dfGeneNames=dfGeneNames.astype(str).applymap(lambda x:x.upper())
dfGeneNames.Synonyms=[str(row).split('|') for row in dfGeneNames.Synonyms.values]
listSynonyms=[elem for row in dfGeneNames.Synonyms.values for elem in row]
#print(dfGeneNames[[igene in row for row in dfGeneNames.Synonyms]].Symbol.values[0])
def FindGeneName(igene):
    retgene=np.nan
    if (igene in set(dfGeneNames.Symbol)) or (igene not in set(listSynonyms)) or (igene in [row[0] for row in degenchlist]):
        return igene
    else:
        retgene=dfGeneNames[[igene in row for row in dfGeneNames.Synonyms]].Symbol.values[0]
    # proceed to rename if the chromosome no is same.    
    chno_retgene=str(dfGeneNames[dfGeneNames.Symbol==retgene].chromosome.values[0])
    chno_igene=str(df_mutfiles[df_mutfiles.Hugo_Symbol==igene].Chromosome.values[0])# to cover simple renaming situations
    return retgene if ((chno_retgene==chno_igene) and (igene in geneset)) else igene



dfC=df_Out_v600[:].copy(deep=True)
genesrenamed=[]
genesadded=[]
#dfGeneNames = dfGeneNames[~dfGeneNames['Symbol'].isin(cystineGeneList)]
#retgene=dfGeneNames[[igene in row for row in dfGeneNames.Synonyms]].Symbol.values[0]
for igene in dfC.index:
    newgene=FindGeneName(igene)
    if (newgene != igene):
        if (newgene in dfC.index.values):
            dfC.loc[newgene]=dfC.loc[newgene].copy()+dfC.loc[igene].copy()
            dfC.drop(igene,inplace=True)
            genesadded=genesadded+[igene]
        else:
            dfC.loc[newgene]=dfC.loc[igene].copy()
            dfC.drop(igene,inplace=True)
            genesrenamed=genesrenamed+[igene]

idx1=list(dfC.index)
idx1=sorted(idx1)
idx1.remove('Total')
idx1=idx1+['Total']
dfC=dfC.loc[idx1]

colist=dfC.columns.sort_values()
colist=[colist[-1]]+list(colist[:-1])
dfC=dfC[colist]

dfC.to_csv(os.path.join(outputdir, 'Genomics_Output_Processed_specific_Braf.txt'),header=True,sep='\t',index_label='Hugo_Symbol')
print("Done")




























