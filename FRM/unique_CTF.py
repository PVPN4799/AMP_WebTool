import pandas as pd
import pandas as pd
import numpy as np
import string
import xlsxwriter
import csv
import re
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
class CTF:
###############################################################################
    all_features = []
    #Amino Acid Symbols
    AALetter=["A","R","N","D","C","E","Q","G","H","I","L","K","M","F","P","S","T","W","Y","V"]
 
    #Defining groups
    _repmat={1:["A",'G','V'],2:['I','L','F','P'],3:['Y','M','T','S'],4:['H','N','Q','W'],5:['R','K'],6:['D','E'],7:['C']}
 
###############################################################################
 
    def _Str2Num(self,proteinsequence):
        """
        substitute protein sequences by their corresponding group numbers.
        
        """
        repmat={}
        counter = len(proteinsequence)-2
        print("Counter: ",counter)
        
        for i in CTF._repmat:          #create a reversed dictionery for amino acids and corresp. group numbers
            for j in CTF._repmat[i]:
                repmat[j]=i
                
                
        res=proteinsequence
        
        for i in repmat:
            res=res.replace(i,str(repmat[i]))
        
        return res, counter, proteinsequence
     
     
###############################################################################
    def CalculateConjointTriad(self,proteinsequence,t):
        """
        Calculate the conjoint triad features from protein sequence.
        
        Useage:
        
        res = CalculateConjointTriad(protein)
        
        Input: protein is a pure protein sequence.
        
        Output is a dict form containing all 343 conjoint triad features.
        """
        field_names = []                   #records feature set names to write to csv file
        field_names.append('Sequence')
        #field_names.append('Target')
        res={}                             #dictionary to store feature set values
        
        proteinnum, counter, proteinseq=CTF._Str2Num(self,proteinsequence)     #subst. protein seq. by group numbers
        for i in range(1,8):
            for j in range(1,8):
                for k in range(1,8):
                    temp=str(i)+str(j)+str(k)
                    field_names.append(temp)
                    res['Sequence']=proteinseq
                    #res['Target']=positive[t]
                    res[temp]=proteinnum.count(temp) / counter
                    a = sum(1 for i in range(len(proteinnum)) if proteinnum.startswith(temp, i)) / counter
                    b=proteinnum.count(temp) / counter
                
                
        print(a)
        print(b)
        t = t+1
        print(t)
        self.all_features.append(res)
    
        return res, field_names, self.all_features, t
 
###############################################################################
 
    def seqinp(self,f):
        ctf = {}
        inFile = pd.read_csv(f)       #read protein sequences
        seq = inFile['Sequence']
        protein = []
        count = []
        counter=0
        for line in seq:                     #check for repetitive sequences
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
        for i in inFile['Target']:          #keep a record of the target value of non-repeating sequences
            if k in count:
                positive.append(i)
            #print(k)
            k=k+1
        '''
        t=0
        test=0
        for i in protein:                   #generate feature sets
            ctf, field_names, self.all_features, t = CTF.CalculateConjointTriad(self,i,t)
            test=test+1
        
        ctf_df = pd.DataFrame(self.all_features)
        return ctf_df
        
        '''
        with open('unique_CTF.csv','w') as csvfile:          #write a csv file with all features
            writer = csv.DictWriter(csvfile, fieldnames = field_names)
            writer.writeheader()
            writer.writerows(all_features)
        '''
