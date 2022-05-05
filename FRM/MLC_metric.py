import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import sklearn
import sklearn.metrics as metrics
import math
def accuracy(a,b):
    count=0
    for i in range(len(a)):
        if a[i] == b[i]:
            count = count + 1
    acc = count/len(a)
    return acc, count

def lab_acc(y_test,pred):
	label = []
	predl = []
	sum_acc = {}
	col = 0
	count = 0
	acc_c = 0
	#calclating average accuracy by collecting values for each value separately
	#and comparing them
	for i in y_test.columns:
	    row = 0
	    for j in y_test.index:
	        label.append(y_test.loc[j,i])
	        predl.append(pred.toarray()[row][col])
	        row = row + 1
	    col = col + 1
	    sum_acc[i], c = accuracy(label,predl)
	    count = count + c
	    acc_c = acc_c + sum_acc[i]
	sum_acc['micro_av'] = acc_c/5
	sum_acc_pd=pd.DataFrame(sum_acc,index=['Accuracy'])
	return sum_acc_pd
def auc_roc(clf,X_test,y_test):
	yhat = clf.predict_proba(X_test)
	auc={}
	plt.plot([0,1], [0,1], linestyle='--', label='No Skill')
	#yhat = yhat[:,1]
	c = ['r','b','g','y','m']
	label = []
	predl = []
	col = 0
	for i in y_test.columns:
	    row = 0
	    for j in y_test.index:
	        label.append(y_test.loc[j,i])
	        predl.append(yhat.toarray()[row][col])
	        row = row + 1
	    col = col + 1
	# keep probabilities for the positive outcome only
	    #predl = predl[predl == 1]
	# calculate roc curves
	    fpr, tpr, thresholds = metrics.roc_curve(label, predl)
	    auc[i] = metrics.roc_auc_score(label, predl)
	# plot the roc curve for the model
	    plt.plot(fpr, tpr, marker=',',color = c[col-1], label= i)
	    label = []
	    predl = []
	# axis labels
	plt.xlabel('False Positive Rate')
	plt.ylabel('True Positive Rate')
	plt.legend()
	# show the plot
	plt.show()
	auc_pd = pd.DataFrame(auc,index=['AUC'])
	return auc_pd
def mcc_lab(y_test,pred):
	mcc = {}
	f1={}
	tp = [0,0]
	fp = [0,0]
	tn = [0,0]
	fn = [0,0]
	label = []
	predl = []
	col = 0
	for i in y_test.columns:
	    row = 0
	    for j in y_test.index:
	        label.append(y_test.loc[j,i])
	        predl.append(pred.toarray()[row][col])
	        row = row + 1
	    col = col + 1
	    mat = metrics.confusion_matrix(label,predl)
	    tn[0], fp[0], fn[0], tp[0] = mat[0,0], mat[0,1], mat[1,0], mat[1,1]
	    tp[1] = tp[1]+tp[0]
	    fn[1] = fn[1]+fn[0]
	    fp[1] = fp[1]+fp[0]
	    tn[1] = tn[1]+tn[0] 
	    den = math.sqrt((tp[0]+fp[0])*(tp[0]+fn[0])*(tn[0]+fp[0])*(tn[0]+fn[0]))
	    num = (tp[0]*tn[0])-(fp[0]*fn[0])
	    f = tp[0]/(tp[0]+(0.5*(fp[0]+fn[0])))
	    mcc[i] = num/den
	    f1[i] = f
	den_mic = math.sqrt((tp[1]+fp[1])*(tp[1]+fn[1])*(tn[1]+fp[1])*(tn[1]+fn[1]))
	num_mic = (tp[1]*tn[1])-(fp[1]*fn[1])
	f1_mic = tp[1]/(tp[1]+(0.5*(fp[1]+fn[1])))
	mcc['micro_av'] = num_mic/den_mic
	f1['micro_av'] = f1_mic
	mcc_pd = pd.DataFrame(mcc,index=['MCC'])
	f1_pd = pd.DataFrame(f1,index=['F1'])
	return mcc_pd, f1_pd
def cv_test_results(clf,X_test,y_test,pred):
	best = np.where(clf.cv_results_['rank_test_accuracy']==1)[0][0]
	cv = clf.n_splits_
	n = [str(a) for a in range(cv)]
	be = {}
	st = {}
	av = {}
	met = ['accuracy','f1_micro','jaccard','h_score']
	for i in met:
		bess = [clf.cv_results_['split'+j+'_test_'+i][best] for j in n]
		be[i] = max(bess)
	
	for i in met:
		nam = 'std_test_'+i
		st[i]= clf.cv_results_[nam][best]
	for i in met:
		nam = 'mean_test_'+i
		av[i]= clf.cv_results_[nam][best]	
	
	mcc_pd, f1_pd = mcc_lab(y_test,pred)
	sum_acc_pd = lab_acc(y_test,pred)
	
	return be,st,av,mcc_pd,f1_pd,sum_acc_pd
			

def classification_report(clf,X_test,y_test,pred):
	mcc_pd, f1_pd = mcc_lab(y_test,pred)
	sum_acc_pd = lab_acc(y_test,pred)
	mcc_pd = mcc_pd.append(sum_acc_pd)
	mcc_pd = mcc_pd.append(f1_pd)
	#code to be used for max 
	best = np.where(clf.cv_results_['rank_test_accuracy']==1)[0][0]
	cv = clf.n_splits_
	n = [str(a) for a in range(cv)]
	be = {}
	st = {}
	met = ['accuracy','f1_micro','jaccard','h_score']
	for i in met:
		bess = [clf.cv_results_['split'+j+'_test_'+i][best] for j in n]
		be[i] = max(bess)
	
	for i in met:
		nam = 'std_test_'+i
		st[i]= clf.cv_results_[nam][best]		

	res={'accuracy': metrics.accuracy_score(y_test,pred),
	'f1_micro':metrics.f1_score(y_test,pred,average='micro'),
	'jaccard':metrics.jaccard_score(y_test,pred,average='micro'),
	'h_score':metrics.hamming_loss(y_test,pred)
	}
	cv_res= {'accuracy':[be['accuracy'],clf.cv_results_['mean_test_accuracy'].max(), st['accuracy']],
         'f1_micro':[be['f1_micro'],clf.cv_results_['mean_test_f1_micro'].max(),st['f1_micro']],
         'jaccard':[be['jaccard'],clf.cv_results_['mean_test_jaccard'].max(),st['jaccard']],
         'h_score':[be['h_score'],clf.cv_results_['mean_test_h_score'].max(),st['h_score']]}
	cv_pd=pd.DataFrame(cv_res,index=["Max","Mean","Std"])
	test_pd=pd.DataFrame(res,index=["Values"]) 
	print("Metrics tested:\n1: Label based accuracy\n2: Matthew's Correlation Coefficient")
	print("3: AUC-ROC\n 4: Example based accuracy\n5: Micro averaged F1 score\n6: Jaccard index")
	print("7: Hamming Loss")
	print()
	print("Label wise metric results:")
	print(mcc_pd)
	auc_pd = auc_roc(clf,X_test,y_test)
	print(auc_pd)
	print()
	print("Sample based metric resuits:")
	print()
	print("Cross Validation")
	print()
	print(cv_pd)
	print()
	print("Testing")
	print()
	print(test_pd)
	return mcc_pd,auc_pd,cv_pd,test_pd
	
