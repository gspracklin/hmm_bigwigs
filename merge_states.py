def merge_different_hmmstates(df,state1=0, state2=2, state3=1):
    "merge 3 HMM states into 2 "
    import pandas as pd

    chr_list=[]
    start_list=[]
    end_list=[]


    for item in df['chrom'].unique():
        chrom_df=df[df['chrom']==item]
        start=1
        print('Starting chr '+item)
        for index, row in chrom_df.iterrows():
            if start == 1:
                if df['state'].iloc[index] == state1 or df['state'].iloc[index] == state2:
                    chr_list.append(df['chrom'].iloc[index])
                    start_list.append(df['start'].iloc[index])
                    start = 0
                    continue
                else:
                    continue
            elif df['state'].iloc[index] == state3:
                end_list.append(df['end'].iloc[(index-1)])
                start=1
            else:continue
        if start == 0:
            end_list.append(df['end'].iloc[(index)])

        if len(chr_list) != len(start_list) or len(start_list) != len(end_list):
            print('Wrong Lengths!')
            break

    keys=['chrom','start','end','state']
    values=[chr_list,start_list,end_list]
    dictionary = dict(zip(keys, values))
    df_merge= pd.DataFrame.from_dict(dictionary)
    return df_merge
