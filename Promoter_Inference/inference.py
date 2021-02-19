import os
import torch
import torch.nn as nn
import torch.nn.functional as F
from .Network import NucleicTransformer
from tqdm import tqdm
import pandas as pd
import numpy as np

nt_int={
"A": 0,
"T": 1,
"G": 2,
"C": 3,}

nts=[
"A",
"T",
"G",
"C"]

def nucleatide2int(nt_sequence,target_length=None):
    int_sequence=[]
    for nt in nt_sequence:
        nt=nt.upper()
        if nt in nt_int:
            int_sequence.append(nt_int[nt])
    int_sequence=np.asarray(int_sequence,dtype='int32')
    if target_length:
        int_sequence=np.pad(int_sequence,(0,target_length-len(int_sequence)),constant_values=-1)
    return int_sequence

def int2nucleotide(nt_sequence,target_length=None):
    seq=''
    for nt in nt_sequence:
        seq+=nts[nt]
    return seq

def get_kmers(sequence,k):
    kmers=[]
    for i in range(len(sequence)-k+1):
        kmers.append(sequence[i:i+k])
    return kmers

class Promoter_Inference():
    def __init__(self):
        self.device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        pass


    def load_models(self, path):
        self.models=[]
        for i in range(5):
            # model=NucleicTransformer(opts.ntoken, opts.nclass, opts.ninp, opts.nhead, opts.nhid,
            #                        opts.nlayers, opts.kmer_aggregation, kmers=[opts.kmers],
            #                        dropout=opts.dropout).to(device)

            model=NucleicTransformer(4, 2, 256, 8, 1024, 6, True, kmers=[7],return_aw=True).to(self.device)

            weights_path=f"{path}/fold{i}top1.ckpt"
            checkpoint=torch.load(weights_path)
            model.load_state_dict(checkpoint)
            model.eval()
            self.models.append(model)

    def _geometric_mean(self, preds):
        gmean=np.ones(preds.shape[1:])

        for pred in preds:
            gmean=gmean*pred

        gmean=gmean**(1/len(preds))
        return gmean

    def predict(self, df, batch_size=32, topk=3):

        #df=pd.read_csv(df_path)

        seqs=[]
        for sequence in df.sequence:
            seqs.append(nucleatide2int(sequence))

        seqs=np.asarray(seqs)

        ensemble_predictions=[]
        preds=[]
        softmax=torch.nn.Softmax(1)

        if len(seqs)<batch_size:
            batches=1
        elif len(seqs)//batch_size==0:
            batches=len(seqs)//batch_size
        else:
            batches=len(seqs)//batch_size+1

        top_kmers=[]
        top_kmer_counts=[]

        for j in tqdm(range(batches)):
            for model in self.models:
                outputs=[]
                with torch.no_grad():
                    x=torch.Tensor(seqs[j*batch_size:(j+1)*batch_size]).to(self.device).long()
                    output, attention_weights=model(x,None)
                    outputs.append(softmax(output).cpu().numpy())

                outputs=np.asarray(outputs)
                outputs=self._geometric_mean(outputs)
                pred=np.argmax(outputs,axis=-1)
            for item in pred:
                preds.append(item)

            attention_weights=attention_weights.cpu().numpy().sum(1)

            #kmer extraction here
            for sequence, prediction, attention_vector in zip(x.cpu().numpy(),pred,attention_weights):
                #print(sequence.shape)
                #exit()
                if prediction==1:
                    kmers=get_kmers(int2nucleotide(sequence),7)
                    sorted_indices=np.argsort(attention_vector)
                    #print(kmers[sorted_indices])
                    #exit()
                    top_kmers_sample=np.asarray(kmers)[sorted_indices][-topk:]
                    for kmer in top_kmers_sample:
                        if kmer not in top_kmers:
                            top_kmers.append(kmer)
                            top_kmer_counts.append(1)
                        else:
                            top_kmer_counts[top_kmers.index(kmer)]=top_kmer_counts[top_kmers.index(kmer)]+1
                    # print(top_kmers)
                    # exit()

            #print(attention_weights.shape)


        preds=np.asarray(preds)
        predictions=[]
        for i in preds:
            if i ==0:
                predictions.append('not promoter')
            elif i==1:
                predictions.append('promoter')
        df['predictions']=predictions
        #df.to_csv('out_path',index=False)
        return df, np.asarray(top_kmers), np.asarray(top_kmer_counts)


#code testing here
# path='promoter_small.csv'
# inference_tool=Promoter_Inference()
# inference_tool.load_models()
# inference_tool.predict(path)
