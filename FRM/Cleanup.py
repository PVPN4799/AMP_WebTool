import pandas as pd
import numpy as np

class Cleanup:
	
	'''
	Cleaning up input dataset prior to using Pse-AAC-5 and Pse-AAC-10 feature representation
	methods on it
	'''
	def reminp(self,f,n):
		
		db = pd.read_csv(f)
		db1 = pd.DataFrame(columns=db.columns)
		db_old = db
		for i in db.index:
			j = db['Sequence'][i]
			if len(j) < n+1:
				db1 = db1.append(db_old.iloc[i])
				db = db.drop(index=i)
			else:
				continue
		#db.reset_index(inplace=True)
		name1 = "DB"+str(n)+".csv"
		db.to_csv(name1)
		db1.to_csv("other.csv")