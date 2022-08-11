import os
import torch
import torch.nn as nn
import torch.nn.functional as F
from .Network import NucleicTransformer
from tqdm import tqdm
import pandas as pd
import numpy as np
from .RNA_feature_generation import generate_RNA_features
import matplotlib.pyplot as plt

def preprocess_inputs(sequence, structures, loops):
    tokens = 'ACGU().BEHIMSX'
    input=[]
    for j in range(len(structures)):
        input_seq=np.asarray([tokens.index(s) for s in sequence])
        input_structure=np.asarray([tokens.index(s) for s in structures[j]])
        #input_loop=np.asarray([tokens.index(s) for s in loops[j]])
        input.append(np.stack([input_seq,input_structure],-1))
    input=np.asarray(input).astype('int')
    return input

def get_distance_mask(L):

    m=np.zeros((3,L,L))


    for i in range(L):
        for j in range(L):
            for k in range(3):
                if abs(i-j)>0:
                    m[k,i,j]=1/abs(i-j)**(k+1)
    return m

class RNA_Inference():
    def __init__(self):
        #self.device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        self.device = 'cpu'

    def load_models(self, path):
        self.models=[]
        for i in range(10):
            # model=NucleicTransformer(opts.ntoken, opts.nclass, opts.ninp, opts.nhead, opts.nhid,
            #                        opts.nlayers, opts.kmer_aggregation, kmers=[opts.kmers],
            #                        dropout=opts.dropout).to(device)

            model=NucleicTransformer(15, 5, 256, 32, 1024, 5, True, kmers=[5],return_aw=True)#.to(self.device)

            weights_path=f"{path}/fold{i}top1.ckpt"
            checkpoint=torch.load(weights_path,map_location=self.device)
            model=nn.DataParallel(model)
            model.load_state_dict(checkpoint)
            model=model.module.cpu()
            model.eval()
            self.models.append(model)

    def _geometric_mean(self, preds):
        gmean=np.ones(preds.shape[1:])

        for pred in preds:
            gmean=gmean*pred

        gmean=gmean**(1/len(preds))
        return gmean

    def predict(self, sequence):
        bp_matrix_seq, structures, loops=generate_RNA_features(sequence)

        input_features={'bpps':bp_matrix_seq,
                        'structures':structures,
                        'loops':loops}

        inputs=preprocess_inputs(sequence, structures, loops)
        #print(inputs.shape)
        outputs=[]
        aws=[]
        with torch.no_grad():
            for model in self.models:
                for j in range(len(bp_matrix_seq)):
                    src=torch.tensor(inputs[j]).long()#.to(self.device)
                    bpp=torch.tensor(bp_matrix_seq[j]).float()#.to(self.device)

                    distance_mask=torch.tensor(get_distance_mask(src.shape[0])).float()#.to(self.device)

                    bpp=torch.cat([bpp.unsqueeze(0),distance_mask],0)
                    output,aw=model(src[:,0].unsqueeze(0),bpp.unsqueeze(0))

                    #outputs.append(output.squeeze().cpu().numpy())
                    #aws.append(aw.squeeze().mean(0).cpu().numpy())

                    outputs.append(output.squeeze().numpy())
                    aws.append(aw.squeeze().mean(0).numpy())
                    # plt.subplot(1,2,1)
                    # plt.imshow(bpp.squeeze()[0].cpu().numpy())
                    # plt.subplot(1,2,2)
                    # plt.imshow(aw.squeeze().mean(0).cpu().numpy())
                    # plt.show()
                    # exit()
        #df.to_csv('out_path',index=False)
        outputs=np.stack(outputs,0).mean(0)
        aws=np.stack(aws,0).mean(0)

        # plt.imshow(outputs.transpose(1,0))
        # plt.show()
        # exit()

        # print(outputs.shape)
        # print(aws.shape)
        # exit()

        return outputs, aws, input_features


#code testing here
# path='promoter_small.csv'
# inference_tool=Promoter_Inference()
# inference_tool.load_models()
# inference_tool.predict(path)
