import pandas as pd
import pandas as pd
import numpy as np
import string
import xlsxwriter
import time
import csv 
import sys
sys.path.append("./propy")
import PyPro
from GetProteinFromUniprot import GetProteinSequence

class PAAC_5:

    all_features = []
    
    def seqinp2(self,f):
        start = time.time()
        paac = {}
        error_sequences = []
        seq = [f]
        protein = []
        error_index = []
        field_names=[]
        count = []
        counter=0
        for line in seq:
            if line in protein:
                counter = counter+1
                continue        
            else:
                count.append(counter)
                protein.append(line)
                counter=counter+1
        positive = []
        k=0
        #print(count)
        '''
        for i in inFile['Target']:    
            if k in count:
                positive.append(i)
            #print(k)
            k=k+1
        '''
        print("Length of target: ",len(positive))
        print("Length of protein data: ",len(protein))

        for i in range(1,26):
            temp = 'PAAC'+str(i)
            field_names.append(temp)

        field_names.append('Sequence')
        #field_names.append('Target')
        count = 0
        length_less = []
        error=1
        for prot in protein:
            error = error+1
            try:
                if len(prot)>=6:
                    proteinsequence = prot
                    DesObject = PyPro.GetProDes(proteinsequence)      # construct a GetProDes object
                    paac = DesObject.GetPAAC(lamda=5,weight=0.05)
                    paac['Sequence']=prot
                    #paac['Target']=target
                    self.all_features.append(paac)

                else:
                    length_less.append(error)
            except LookupError:
                count = count+1
                error_sequences.append(prot)
                error_index.append(error)
        print("Number of error sequences: ",count)
        print("Error sequences: ",error_sequences)
        print("Error indices: ",error_index)
        print("Sequence<=21 indices: ",length_less)


        '''
        with open('uniquePAAC_5.csv','w') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames = field_names)
            writer.writeheader()
            writer.writerows(all_features)
        '''

        stop = time.time()
        print(f"Runtime of the program is {stop-start}")
        pc5_df = pd.DataFrame(self.all_features)
        return pc5_df

    def seqinp(self,f):
        start = time.time()
        paac = {}
        error_sequences = []
        inFile = pd.read_csv(f)
        seq = inFile['Sequence']
        protein = []
        error_index = []
        field_names=[]
        count = []
        counter=0
        for line in seq:
            if line in protein:
                counter = counter+1
                continue        
            else:
                count.append(counter)
                protein.append(line)
                counter=counter+1
        positive = []
        k=0
        #print(count)
        '''
        for i in inFile['Target']:    
            if k in count:
                positive.append(i)
            #print(k)
            k=k+1
        '''
        print("Length of target: ",len(positive))
        print("Length of protein data: ",len(protein))

        for i in range(1,26):
            temp = 'PAAC'+str(i)
            field_names.append(temp)

        field_names.append('Sequence')
        #field_names.append('Target')
        count = 0
        length_less = []
        error=1
        for prot in protein:
            error = error+1
            try:
                if len(prot)>=6:
                    proteinsequence = prot
                    DesObject = PyPro.GetProDes(proteinsequence)      # construct a GetProDes object
                    paac = DesObject.GetPAAC(lamda=5,weight=0.05)
                    paac['Sequence']=prot
                    #paac['Target']=target
                    self.all_features.append(paac)

                else:
                    length_less.append(error)
            except LookupError:
                count = count+1
                error_sequences.append(prot)
                error_index.append(error)
        print("Number of error sequences: ",count)
        print("Error sequences: ",error_sequences)
        print("Error indices: ",error_index)
        print("Sequence<=21 indices: ",length_less)


        '''
        with open('uniquePAAC_5.csv','w') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames = field_names)
            writer.writeheader()
            writer.writerows(all_features)
        '''

        stop = time.time()
        print(f"Runtime of the program is {stop-start}")
        pc5_df = pd.DataFrame(self.all_features)
        return pc5_df
