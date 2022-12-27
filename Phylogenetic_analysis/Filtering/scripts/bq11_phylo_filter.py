#required packages
import pandas as pd


def subsampling(lindf, nonconvacc, convsampsize, nonconvsize, outfile, minconv):
    
    """
    - lindf: GISAID metadata dataframe with all filtered accessions of a sublineage (e.g. BA.1, BA.2, BA.4, BA.5)
    - nonconvacc: list of accessions to consider from the 'background' set of sequences (non-converging)
    - convsampsize: number of sequences to sample from each 'converging' PANGO lineage
    - nonconvsize: number of sequences to sample from the full 'background' set
    - outfile: the output file name
    - minconv: the minimum number of key Spike substitutions for considering a sequence 'converging' 
    """
    
    #filter background sequences by list 
    newacc = [x for x in list(pd.read_csv(nonconvacc, header=None)[0])]
            
    thedf_f_nonconv = lindf[lindf.ID.isin(newacc)]
    
    #get converging sequences
    thedf_conv = lindf[lindf['conv-sites']>=minconv]
    
    thedf_conv_sub1 = pd.DataFrame()
    
    #for each converging PANGO lineage get defined random sample
    for p in set(thedf_conv.pango):
        if len(thedf_conv[thedf_conv.pango == p]) > convsampsize:
            thedf_conv_sub1 = pd.concat([thedf_conv_sub1, thedf_conv[thedf_conv.pango == p].sample(n=convsampsize, random_state=1) ])
        else:
            thedf_conv_sub1 = pd.concat([thedf_conv_sub1, thedf_conv[thedf_conv.pango == p] ])    
    
    #for background accessions weigth the sampling by the frequency of each PANGO lineage
    thedf_f_nonconv_sub1 = thedf_f_nonconv.sample(n = nonconvsize, weights = thedf_f_nonconv.groupby('pango')['pango'].transform('count'))
    
    #put together
    thedf_f_all_sub1 = pd.DataFrame()
    
    thedf_f_all_sub1 = pd.concat([thedf_conv_sub1, thedf_f_nonconv_sub1])

    
    print(len(thedf_f_all_sub1))
    
    #export files    
    thedf_f_all_sub1.to_csv(outfile, sep ='\t', index = False)
    thedf_f_all_sub1dats = thedf_f_all_sub1[['ID', 'c-date']]
    thedf_f_all_sub1dats.columns = ['name', 'date']
    thedf_f_all_sub1dats.to_csv(outfile.replace('.tsv', '_dates.csv'), index = False)
    thedf_f_all_sub1dats[['name']].to_csv(outfile.replace('.tsv', '_IDs.txt'), index = False, header = None)
    
    return thedf_f_all_sub1