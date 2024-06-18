#!/bin/python3
# Written by Roberto Lins

import pickle, joblib
import numpy as np
    
def train_base_models(x, y, models):
    trained_models = []
    for model in models:
        model.fit(x, y)
        trained_models.append((model.__class__.__name__, model))
    return trained_models

def super_learner_predictions(X_test, models, meta_model):
    model_predictions = {}
    for name, model in models:
        yhat = model.predict(X_test)
        model_predictions[name] = yhat
    meta_X = np.column_stack(list(model_predictions.values()))
    super_learner_preds = meta_model.predict(meta_X)
    model_predictions["SuperLearner"] = super_learner_preds
    return model_predictions

def predictor(basemodels, superlearner, x, y, rosetta_features):
    models = [joblib.load(filename) for filename in basemodels]
    base_models = train_base_models(x, y, models)
    with open(superlearner, 'rb') as f:
        meta_model = pickle.load(f)
    model_predictions = super_learner_predictions(rosetta_features, base_models, meta_model)
    return model_predictions['SuperLearner'][0]
