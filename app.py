import pandas as pd
import numpy as np
import time
from sklearn.model_selection import train_test_split
from sklearn.model_selection import cross_validate
from sklearn.model_selection import GridSearchCV
from sklearn import svm
import sklearn.metrics as metrics
from skmultilearn.adapt import MLkNN
from sklearn.metrics import classification_report

import matplotlib.pyplot as plt
from statistics import mean

from FRM import uniqueAAC
from FRM import uniquePAAC_5
from FRM import uniquePAAC_10
from FRM import Cleanup

from FRM import MLC_metric as mlc
from FRM.H_S import hamming_score

from joblib import dump, load

from flask import Flask, request, jsonify, render_template
import __main__
__main__.hamming_score = hamming_score

app = Flask(__name__)


model11 = load("Lvl_1_PC10.joblib")
model12 = load("Lvl_1_PC5.joblib")
model13 = load("Lvl_1_AAC.joblib")
model21 = load("Lvl_2_PC10.joblib")
model22 = load("Lvl_2_PC5.joblib")
model23 = load("Lvl_2_AAC.joblib")
model31 = load("Lvl_3_PC10.joblib")
model32 = load("Lvl_3_AAC.joblib")


@app.route('/')
def home():
	return render_template('index.html')

@app.route('/predict',methods=['POST'])
def predict():
	seq = [x for x in request.form.values()][0]
	flag = 0
	output = ""
	AALetter = list("ARNDCEQGHILKMFPSTWYV")
	for i in seq:
		if i in AALetter:
			continue
		else:
			flag = 1
	if flag == 0:
	
		if len(seq) > 10:
			obj7 = uniquePAAC_10.PAAC_10()
			obj7.all_features=[]
			amp_aac = obj7.seqinp2(seq)
			X7 = amp_aac.drop(columns='Sequence')
			pred = model11.predict(X7)
		elif len(seq) > 5:
			obj7 = uniquePAAC_5.PAAC_5()
			obj7.all_features=[]
			amp_aac = obj7.seqinp2(seq)
			X7 = amp_aac.drop(columns='Sequence')
			pred = model12.predict(X7)
		else:
			obj7 = uniqueAAC.AAC()
			obj7.all_features=[]
			amp_aac = obj7.seqinp2(seq)
			X7 = amp_aac.drop(columns='Sequence')
			pred = model13.predict(X7)
		
		if pred[0] == 0:
			output = "Non-AMP"
		else:
			output = "AMP \n"
			if len(seq) > 10:
				obj7 = uniquePAAC_10.PAAC_10()
				obj7.all_features=[]
				amp_aac = obj7.seqinp2(seq)
				X7 = amp_aac.drop(columns='Sequence')
				pred2 = model21.predict(X7)
			elif len(seq) > 5:
				obj7 = uniquePAAC_5.PAAC_5()
				obj7.all_features=[]
				amp_aac = obj7.seqinp2(seq)
				X7 = amp_aac.drop(columns='Sequence')
				pred2 = model22.predict(X7)
			else:
				obj7 = uniqueAAC.AAC()
				obj7.all_features=[]
				amp_aac = obj7.seqinp2(seq)
				X7 = amp_aac.drop(columns='Sequence')
				pred2 = model23.predict(X7)
		
			if pred2.toarray()[0][0] == 1:
				output = output + "Antibacterial \n"
			if pred2.toarray()[0][1] == 1:
				output = output + "Antiviral \n"
			if pred2.toarray()[0][2] == 1:
				output = output + "Antifungal \n"
			if pred2.toarray()[0][0] == 1:
				if len(seq) > 10:
					obj7 = uniquePAAC_10.PAAC_10()
					obj7.all_features=[]
					amp_aac = obj7.seqinp2(seq)
					X7 = amp_aac.drop(columns='Sequence')
					pred3 = model31.predict(X7)
				else:
					obj7 = uniqueAAC.AAC()
					obj7.all_features=[]
					amp_aac = obj7.seqinp2(seq)
					X7 = amp_aac.drop(columns='Sequence')
					pred3 = model32.predict(X7)
		
				if pred3.toarray()[0][0] == 1:
					output = output + "Anti Gram Positive \n"
				if pred3.toarray()[0][1] == 1:
					output = output + "Anti Gram Negative \n"
	else:
		output = "Invalid Sequence"

	return render_template('index.html', prediction_text = "Sequence is {}".format(output))

if __name__ == "__main__":
	app.run(debug=True)
		