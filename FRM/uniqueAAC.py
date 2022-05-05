import pandas as pd
import numpy as np
import string
import xlsxwriter
import csv
'''
inFile = pd.read_csv("trial_total.csv")
seq = inFile['Sequence']
listLines = []
count = []
counter=0
for line in seq:
    if line in listLines:
        counter = counter+1
        continue
        
    else:
        count.append(counter)
        listLines.append(line)
        counter=counter+1
'''
class AAC:
    
    all_features = []
###############################################################################
    AALetter=["A","R","N","D","C","E","Q","G","H","I","L","K","M","F","P","S","T","W","Y","V"]
   
###############################################################################
    def AminoAcidComp(self,proteinsequence,t):
        """
        Calculate the conjoint triad features from protein sequence.
	
        Useage:
	
        res = CalculateConjointTriad(protein)
	
        Input: protein is a pure protein sequence.
	
        Output is a dict form containing all 343 conjoint triad features.
        """
        field_names = []
        field_names.append('Sequence')
        #field_names.append('Target')
        res={}
        '''
        proteinnum, counter, proteinseq=_Str2Num(proteinsequence)
        for i in range(1,5):
            for j in range(1,5):
                for k in range(1,5):
                    temp=str(i)+str(j)+str(k)
                    field_names.append(temp)
                    res['Sequence']=proteinseq
                    res['Target']=positive[t]
                    res[temp]=proteinnum.count(temp) / counter
    
        t = t+1
        all_features.append(res)
        '''
    
    
        for i in AAC.AALetter:
            temp = str(i)
            field_names.append(temp)
            res['Sequence']=proteinsequence
            #res['Target']=positive[t]
            counter=0    
            for j in proteinsequence:
                if j==i:
                    counter = counter+1    
            res[temp] = counter / len(proteinsequence)
        
        t = t+1
        self.all_features.append(res)    
    
        return res, field_names, self.all_features
 
###############################################################################

    def seqinp2(self,f):
        aac = {}
        seq = [f]
        protein = []
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
        t=0
        test=0
        for i in protein:
            aac, field_names, self.all_features = AAC.AminoAcidComp(self,i,t)
            test=test+1
        aac_df = pd.DataFrame(self.all_features)
        return aac_df
        '''
        with open('uniqueAAC.csv','w') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames = field_names)
            writer.writeheader()
            writer.writerows(all_features)
        '''
 
    def seqinp(self,f):
        aac = {}
        inFile = pd.read_csv(f)
        seq = inFile['Sequence']
        protein = []
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
        t=0
        test=0
        for i in protein:
            aac, field_names, self.all_features = AAC.AminoAcidComp(self,i,t)
            test=test+1
        aac_df = pd.DataFrame(self.all_features)
        return aac_df
        '''
        with open('uniqueAAC.csv','w') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames = field_names)
            writer.writeheader()
            writer.writerows(all_features)
        '''
