import pandas as pd
import numpy as np
import string
import re
from itertools import combinations

class CTD:
    
    '''
    
    CTD calculation for 7 different physical/chemical properties
    The 20 amino acids have been grouped into 3 groups for each property
    This is known specifically as the CTD-3 feature representation method
    Data for the different groupings was collected from Govindan et al (2011)
    
    '''
    
    all_features = []
    #Amino Acid Symbols
    AALetter=["A","R","N","D","C","E","Q","G","H","I","L","K","M","F","P","S","T","W","Y","V"]
    
    #Defining groups
    parameters = {
    
    #based on charge: 1 - neutral; 2 - negative; 3 - positive
    'chg':{1:['A','C','F','G','H','I','L','M','N','P','Q','S','T','V','W','Y'],2:['D','E'],3:['K','R']},
    
    #based on hydrophobicity: 1 - hydrophobic; 2 - neutral; 3- polar
    'hph':{1:['C','F','I','L','M','V','W'],2:['A','G','H','P','S','T','Y'],3:['D','E','K','N','Q','R']},
    
    #based on normalized VDW volume: 1-3; 1 low, 3 high
    'vdw':{1:['A','C','D','G','P','S','T'],2:['E','I','L','N','Q','V'],3:['F','H','K','M','R','W','Y']},
    
    #based on polarity: 1-3; 1 low, 3 high
    'pol1':{1:['C','F','I','L','M','V','W','Y'],2:['A','G','P','S','T'],3:['D','E','H','K','N','Q','R']},
    
    #based on polarizability: 1-3; 1 low, 3 high
    'pol2':{1:['A','D','G','S','T'],2:['C','E','I','L','N','Q','P','V',],3:['F','H','M','K','R','W','Y']},
    
    #based on secondary structure: 1-coil, 2-helix;3-strand
    'sec':{1:['D','G','N','P','S'],2:['A','E','H','K','L','M','Q','R'],3:['C','F','I','T','V','W','Y']},
    
    #based on solvent accessibility: 1-buried, 2-intermediate, 3-exposed
    'solv':{1:['A','C','F','G','I','L','V','W'],2:['H','M','P','S','T','Y'],3:['D','E','K','N','R','Q']}
    }
 
###############################################################################
 
    def Str2Num(self,proteinsequence):
        
        """
        substitute protein sequences by their corresponding group numbers.
        Returns a dict containing numerical sequences for 7 properties
        
        """
        pram = {'chg':{},'hph':{},'vdw':{},'pol1':{},'pol2':{},'sec':{},'solv':{}}
        nseq = {'chg':'','hph':'','vdw':'','pol1':'','pol2':'','sec':'','solv':''}
        dd = {}
        for i in CTD.parameters:        #create a reversed dictionery for amino acids and corresp. group numbers
            for j in CTD.parameters[i]:
                for l in CTD.parameters[i][j]:
                    pram[i][l]=j
                
        res=proteinsequence
        
        for i in pram:
            for j in pram[i]:
                res=res.replace(j,str(pram[i][j]))
            nseq[i]=res
            res = proteinsequence
        
        return nseq, proteinsequence
    
    def ctd_calc(self,nseq,proteinsequence,sr):
        
        ind = ['1','2','3']
        res = {"Sequence": proteinsequence}
        comb = list(combinations(ind,2))
        
        '''
        Composition calculation - produces 3 features for each property:
        C_1: composition of amino acids belonging to group 1
        C_2: composition of amino acids belonging to group 2
        C_3: composition of amino acids belonging to group 3
        '''
        
        for j in nseq:
            for i in ind:
                counter = 0
                for k in nseq[j]:
                    if k == i:
                        counter = counter + 1
                res[j+'_C'+i] = counter/len(nseq[j])
        
        '''
        Transition calculation - produces 3 features for each property:
        T_12: % frequency of sequence transitions from group 1 to group 2 OR from group 2 to group 1
        T_13: % frequency of sequence transitions from group 1 to group 3 OR from group 3 to group 1
        T_23: % frequency of sequence transitions from group 2 to group 3 OR from group 3 to group 2
        '''
        
        for j in nseq:
            for i in comb:
                counter = 0
                for k in range(len(nseq[j])-1):
                    if nseq[j][k] == i[0] and nseq[j][k+1] == i[1]:
                        counter = counter + 1
                    if nseq[j][k] == i[1] and nseq[j][k+1] == i[0]:
                        counter = counter + 1
                    res[j+'_T'+i[0]+i[1]] = counter/(len(nseq[j])-1)
        '''
        Distribution calculation - similar to calculating quartiles for a range of data
        DQ0: Position percentage of first occurrence of the group
        DQ1: Position percentage at which 25% elements of the group are covered
        DQ2: Position percentage at which 50% elements of the group are covered
        DQ3: Position percentage at which 75% elements of the group are covered
        DQ4: Position percentage at which 100% elements of the group are covered
        
        Most of the actual values are in decimals so they have been rounded off to get sequence positions
        '''
        
        loc = {0:0}
        for j in nseq:
            for i in ind:
                counter = 0
                for k in range(len(nseq[j])):
                    if nseq[j][k] == i:
                        counter = counter + 1
                        loc[counter] = k+1
                if counter == 0:
                    q0 = 0
                    q1 = 0
                    q2 = 0
                    q3 = 0
                    q4 = 0
                
                elif counter == 1:
                    q0 = 100*(nseq[j].find(i)/len(nseq[j]))
                    q1 = 100*(nseq[j].find(i)/len(nseq[j]))
                    q2 = 100*(nseq[j].find(i)/len(nseq[j]))
                    q3 = 100*(nseq[j].find(i)/len(nseq[j]))
                    q4 = 100*(nseq[j].find(i)/len(nseq[j]))
                
                else:
                    q0 = 100*(nseq[j].find(i)/len(nseq[j]))
                    l1 = round(0.25*counter)
                    l2 = round(0.5*counter)
                    l3 = round(0.75*counter)
                    l4 = counter
                    q1 = 100*(loc[l1]/len(nseq[j]))
                    q2 = 100*(loc[l2]/len(nseq[j]))
                    q3 = 100*(loc[l3]/len(nseq[j]))
                    q4 = 100*(loc[l4]/len(nseq[j]))
                res[j+'_DQ0_'+i] = q0
                res[j+'_DQ1_'+i] = q1
                res[j+'_DQ2_'+i] = q2
                res[j+'_DQ3_'+i] = q3
                res[j+'_DQ4_'+i] = q4
        
        self.all_features.append(res)
        return self.all_features
    
    def seqinp(self,f):
        file = pd.read_csv(f)
        nseq = {}
        for i in file.index:
            seq=file.loc[i,'Sequence']
            nseq, i = CTD.Str2Num(self,seq)
            self.all_features = CTD.ctd_calc(self,nseq,seq,i)
        ctd_df=pd.DataFrame(self.all_features)
        return ctd_df
    
    def csv_out(self):
        #same as seqinp, but outputs csv file with the AMP data as input
        file = pd.read_csv("AMP_Data_Cleaned.csv")
        nseq = {}
        for i in file.index:
            seq=file.loc[i,'Sequence']
            nseq, i = CTD.Str2Num(self,seq)
            self.all_features = CTD.ctd_calc(self,nseq,seq,i)
        ctd_df=pd.DataFrame(self.all_features)
        ctd_df.to_csv("CTD.csv")
                    
                        
        