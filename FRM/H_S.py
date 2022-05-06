def hamming_score(y_test,pred):
    return 1-metrics.hamming_loss(y_test,pred)