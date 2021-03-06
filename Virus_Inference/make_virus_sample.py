import pandas as pd
import numpy as np

df=pd.read_csv('fullset_test.csv',header=None)

sequences=np.asarray(df.iloc[:,1].values)


virus_indices=df.iloc[:,2]==1
non_virus_indices=df.iloc[:,2]==0

sample=list(np.random.choice(sequences[virus_indices],30))+\
list(np.random.choice(sequences[non_virus_indices],90))

sample_df=pd.DataFrame(columns=['sequence'])
sample_df['sequence']=sample

sample_df.to_csv('virus_sample.csv',index=False)
