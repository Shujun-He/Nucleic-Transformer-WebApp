import pandas as pd
import numpy as np

df=pd.read_csv('promoter_sample.csv')

sample_df=df['sequence']

sample_df.to_csv('promoter_sample.csv',index=False)
