# -*- coding: utf-8 -*-
"""
Created on Sun Mar  3 22:42:58 2024

@author: Arjun
"""

#%%
# Basics
import os
import sys
import time
from matplotlib import pyplot as plt
import numpy as np
import arviz as az  # Visualization
import pytensor  # Graph-based tensor library
import jax
import hssm

# import ssms.basic_simulators # Model simulators
import hddm_wfpt
import bambi as bmb

# Setting float precision in pytensor
pytensor.config.floatX = "float32"

from jax.config import config

import pandas as pd
hssm.set_floatX("float32")
config.update("jax_enable_x64", False)



#%%
dataset = pd.read_csv('G:/My Drive/Colab Notebooks/hddm_docker/data/foragingDataNFLW6.csv')
dataset['response'].replace(to_replace = 0, value = -1.0, inplace=True)
dataset['response'].replace(to_replace = 1, value = 1.0, inplace=True)
dataset.rename(columns = {'subj_idx':'participant_id'}, inplace = True)
dataset.head(320)

#%%
# removing trials with RTs less than 350ms and greater than 2s
subset33 = dataset
subsetall = subset33[(subset33["rt"] > 0.35)  & (subset33["rt"] < 2)]
subsetall.head(200)

#%%
# MODEL 1
model_1 = hssm.HSSM(
    data=subsetall,
    model = "ddm",
    )

#%%
hddm_11 = model_1.sample(
    sampler="nuts_numpyro",  # type of sampler to choose, 'nuts_numpyro', 'nuts_blackjax' of default pymc nuts sampler
    cores=2,  # how many cores to use
    chains=2,  # how many chains to run
    draws=1000,  # number of draws from the markov chain
    tune=1000,  # number of burn-in samples
    idata_kwargs=dict(log_likelihood=True),  # return log likelihood
)  # mp_ctx="forkserver")

# #%%
# az.summary(hddm_11)

# #%%
# az.plot_trace(hddm_11);


#%%
# MODEL 2
model_2 = hssm.HSSM(
    data=subsetall,
    model="ddm",
    include=[
        {
            "name": "v",
            "prior": {
                "reward": {"name": "Normal", "mu": 4.0, "sigma": 3.0},
            },
            "formula": "v ~ 1 + (1|participant_id)+ reward",
            "link": "identity",
        },

    ],
)

#%%
hddm_2 = model_2.sample(
    sampler="nuts_numpyro",  # type of sampler to choose, 'nuts_numpyro', 'nuts_blackjax' of default pymc nuts sampler
    cores=4,  # how many cores to use
    chains=4,  # how many chains to run
    draws=8000,  # number of draws from the markov chain
    tune=3000,  # number of burn-in samples
    idata_kwargs=dict(log_likelihood=True),  # return log likelihood
)  # mp_ctx="forkserver")


# #%%
# az.summary(hddm_2)

# #%%
# az.plot_trace(hddm_2);


#%%
# MODEL 3
model_3 = hssm.HSSM(
    data=subsetall,
    model="ddm",
    include=[
        {
            "name": "v",
            "prior": {
                "reward": {"name": "Normal", "mu": 4.0, "sigma": 3.0},
            },
            "formula": "v ~ 1 + (1|participant_id)+ tt + stress2 + reward",
            "link": "identity",
        },

    ],
)

#%%
hddm_3 = model_3.sample(
    sampler="nuts_blackjax",  # type of sampler to choose, 'nuts_numpyro', 'nuts_blackjax' of default pymc nuts sampler
    cores=4,  # how many cores to use
    chains=4,  # how many chains to run
    draws=10000,  # number of draws from the markov chain
    tune=4000,  # number of burn-in samples
    idata_kwargs=dict(log_likelihood=True),  # return log likelihood
)  # mp_ctx="forkserver")


# #%%
# hddm_3_summary = az.summary(hddm_3)

# #%%
# model_3.plot_trace();

# #%%
# model_3.plot_posterior_predictive()


#%%
# MODEL 3A
model_3a = hssm.HSSM(
    data=subsetall,
    model="ddm",
    include=[
        {
            "name": "v",
            "prior": {
                "reward": {"name": "Normal", "mu": 4.0, "sigma": 3.0},
            },
            "formula": "v ~ 1 + (1|participant_id)+ tt + stress2 + reward + stai_trait",
            "link": "identity",
        },

    ],
)

#%%
hddm_3a = model_3a.sample(
    sampler="nuts_numpyro",  # type of sampler to choose, 'nuts_numpyro', 'nuts_blackjax' of default pymc nuts sampler
    cores=4,  # how many cores to use
    chains=4,  # how many chains to run
    draws=10000,  # number of draws from the markov chain
    tune=4000,  # number of burn-in samples
    idata_kwargs=dict(log_likelihood=True),  # return log likelihood
)  # mp_ctx="forkserver")


# #%%
# hddm_3a_summary = az.summary(hddm_3a)

# #%%
# model_3a.plot_posterior_predictive()


# #%%
# compare_data = az.compare(
#     {
#         "reward": model_2.traces,
#         "reward + tt + stress": model_3.traces,
#         "reward + tt + stress + tanx": model_3a.traces
#     }
# )

# compare_data

#%%
# MODEL 4
model_4 = hssm.HSSM(
    data=subsetall,
    model="ddm",
    include=[
        {
            "name": "v",
            "prior": {
                "reward": {"name": "Normal", "mu": 2.0, "sigma": 3.0},
            },
            "formula": "v ~ 1 + (1|participant_id)+ tt + stress2 + reward + stai_trait",
            "link": "identity",
        },
        {
           "name": "a",
           "prior": {
               "reward": {"name": "Normal", "mu": 2.0, "sigma": 3.0},
           },
           "formula": "a ~ 1 + tt + stress2 + reward + stai_trait",
           "link": "identity", 
            
            
            },

    ],
)



#%%

# MODEL 4
model_4 = hssm.HSSM(
    data=subsetall,
    model="ddm",
    include=[
        {
            "name": "v",
            "formula": "v ~ 1 + (1|participant_id) + tt + stress2 + reward + stai_trait",
            "link": "identity",
            "prior": {
                "reward": {"name": "Normal", "mu": 2.0, "sigma": 3.0},
                "stai_trait": {"name": "Normal", "mu": 0.0, "sigma": 3.0},
                "tt": {"name": "Uniform", "lower": -2.0, "upper": 2.0},
                "stress2": {"name": "Uniform", "lower": -2.0, "upper": 2.0},
                "Intercept": {"name": "Uniform", "lower": -2.0, "upper": 2.0},
                "1|participant_id": {"name": "Normal", "mu": 0, "sigma": {"name": "HalfNormal", "sigma": 2.0}}
            },
            "bounds": (-10.0, 10.0),
            
        },
        {
            "name": "z",
            "formula": "z ~ 1 + (1|participant_id)",
            "link": "identity",
            "prior": {
                "Intercept": {"name": "Uniform", "lower": 0.0, "upper": 1.0},
                "1|participant_id": {"name": "Normal", "mu": 0, "sigma": {"name": "HalfNormal", "sigma": 2.0}},
            },
            "bounds": (0.0, 1.0),
            
            },
      
    ],
)
model_4

#%%
hddm_4 = model_4.sample(
    sampler="nuts_numpyro",  # type of sampler to choose, 'nuts_numpyro', 'nuts_blackjax' of default pymc nuts sampler
    cores=4,  # how many cores to use
    chains=4,  # how many chains to run
    draws=12000,  # number of draws from the markov chain
    tune=4000,  # number of burn-in samples
    idata_kwargs=dict(log_likelihood=True),  #  return log likelihood
)  # mp_ctx="forkserver")



# #%%
# hddm_4_summary = az.summary(hddm_4)

# #%%
# model_4.plot_posterior_predictive()




#%%

# MODEL 5
model_5 = hssm.HSSM(
    data=subsetall,
    model="ddm",
    include=[
        {
            "name": "v",
            "formula": "v ~ 1 + (1|participant_id) + tt + stress2 + reward + stai_trait",
            "link": "identity",
            "prior": {
                "reward": {"name": "Normal", "mu": 2.0, "sigma": 3.0},
                "stai_trait": {"name": "Normal", "mu": 0.0, "sigma": 3.0},
                "tt": {"name": "Uniform", "lower": -2.0, "upper": 2.0},
                "stress2": {"name": "Uniform", "lower": -2.0, "upper": 2.0},
                "Intercept": {"name": "Uniform", "lower": -2.0, "upper": 2.0},
                "1|participant_id": {"name": "Normal", "mu": 0, "sigma": {"name": "HalfNormal", "sigma": 2.0}}
            },
            "bounds": (-np.inf, np.inf),
            
        },
        {
            "name": "a",
            "formula": "a ~ 1 + (1|participant_id)",
            "link": "identity",
            "prior": {
                "Intercept": {"name": "Gamma", "mu": 1.5, "sigma": 0.75},
                "1|participant_id": {"name": "Normal", "mu": 0.0, "sigma": {"name": "Weibull", "alpha": 1.5, "beta": 0.3}},
            },
            "bounds": (0.0, np.inf),
            
            },
      
    ],
)
model_5

#%%
hddm_5 = model_5.sample(
    sampler="nuts_numpyro",  # type of sampler to choose, 'nuts_numpyro', 'nuts_blackjax' of default pymc nuts sampler
    cores=4,  # how many cores to use
    chains=4,  # how many chains to run
    draws=10000,  # number of draws from the markov chain
    tune=4000,  # number of burn-in samples
    idata_kwargs=dict(log_likelihood=True),  #  return log likelihood
)  # mp_ctx="forkserver")


#%%
hddm_5_summary = az.summary(hddm_5)

#%%

# MODEL 6
model_6 = hssm.HSSM(
    data=subsetall,
    model="ddm",
    include=[
        {
            "name": "v",
            "formula": "v ~ 1 + (1|participant_id) + tt + stress2 + reward + stai_trait",
            "link": "identity",
            "prior": {
                "reward": {"name": "Normal", "mu": 2.0, "sigma": 3.0},
                "stai_trait": {"name": "Normal", "mu": 0.0, "sigma": 3.0},
                "tt": {"name": "Uniform", "lower": -2.0, "upper": 2.0},
                "stress2": {"name": "Uniform", "lower": -2.0, "upper": 2.0},
                "Intercept": {"name": "Uniform", "lower": -2.0, "upper": 2.0},
                "1|participant_id": {"name": "Normal", "mu": 0, "sigma": {"name": "HalfNormal", "sigma": 2.0}}
            },
            "bounds": (-np.inf, np.inf),
            
        },
        {
            "name": "a",
            "formula": "a ~ 1 + (1|participant_id) + reward + tt + stress2 + stai_trait",
            "link": "identity",
            "prior": {
                "reward": {"name": "Uniform", "lower": 0.01, "upper": 1.0},
                "stai_trait": {"name": "Uniform", "lower": 0.01, "upper": 1.0},
                "tt": {"name": "Uniform", "lower": 0.01, "upper": 0.5},
                "stress2": {"name": "Uniform", "lower": 0.01, "upper": 0.5},
                "Intercept": {"name": "Gamma", "mu": 1.5, "sigma": 0.75},
                "1|participant_id": {"name": "Normal", "mu": 0.0, "sigma": {"name": "Weibull", "alpha": 1.5, "beta": 0.3}},
            },
            "bounds": (0.0, np.inf),
            
            },
      
    ],
)
model_6

#%%
hddm_6 = model_6.sample(
    sampler="nuts_numpyro",  # type of sampler to choose, 'nuts_numpyro', 'nuts_blackjax' of default pymc nuts sampler
    cores=4,  # how many cores to use
    chains=4,  # how many chains to run
    draws=10000,  # number of draws from the markov chain
    tune=4000,  # number of burn-in samples
    idata_kwargs=dict(log_likelihood=True),  #  return log likelihood
)  # mp_ctx="forkserver")


#%%
az.summary(hddm_6)


#%%
# PICKLE HSSM FILES
import pickle

with open('HSSM_model_3a.pkl', 'wb') as f:
    pickle.dump({'model': model_3a, 'trace': model_3a.traces, 'posterior': model_3a.traces.posterior}, f)

#%%
summary_6 = model_6.summary()

#%%
compare_data = az.compare(
    {
        "reward": model_2.traces,
        "reward + tt + stress": model_3.traces,
        "reward + tt + stress + tanx": model_3a.traces,
        "rtst+a all": model_6.traces
    }
)

compare_data


#%%
# MODEL 7
model_7 = hssm.HSSM(
    data=subsetall,
    model="ddm",
    include=[
        {
            "name": "v",
            "formula": "v ~ 1 + (1|participant_id) + tt + stress2 + reward + stai_trait",
            "link": "identity",
            "prior": {
                "reward": {"name": "Normal", "mu": 2.0, "sigma": 3.0},
                "stai_trait": {"name": "Normal", "mu": 0.0, "sigma": 3.0},
                "tt": {"name": "Uniform", "lower": -2.0, "upper": 2.0},
                "stress2": {"name": "Uniform", "lower": -2.0, "upper": 2.0},
                "Intercept": {"name": "Uniform", "lower": -2.0, "upper": 2.0},
                "1|participant_id": {"name": "Normal", "mu": 0, "sigma": {"name": "HalfNormal", "sigma": 2.0}}
            },
            "bounds": (-np.inf, np.inf),
            
        },
        {
            "name": "a",
            "formula": "a ~ 1 + (1|participant_id) + reward + tt + stress2 + stai_trait",
            "link": "identity",
            "prior": {
                "reward": {"name": "Uniform", "lower": -0.4, "upper": 0.4},
                "stai_trait": {"name": "Uniform", "lower": -0.2, "upper": 0.2},
                "tt": {"name": "Uniform", "lower": -0.3, "upper": 0.3},
                "stress2": {"name": "Uniform", "lower": -0.3, "upper": 0.3},
                "Intercept": {"name": "Gamma", "mu": 1.5, "sigma": 0.75},
                "1|participant_id": {"name": "Normal", "mu": 0.0, "sigma": {"name": "Weibull", "alpha": 1.5, "beta": 0.3}},
            },
            "bounds": (0.0, np.inf),
            
            },
      
    ],
)
model_7

#%%
hddm_7 = model_7.sample(
    sampler="nuts_numpyro",  # type of sampler to choose, 'nuts_numpyro', 'nuts_blackjax' of default pymc nuts sampler
    cores=4,  # how many cores to use
    chains=4,  # how many chains to run
    draws=10000,  # number of draws from the markov chain
    tune=4000,  # number of burn-in samples
    idata_kwargs=dict(log_likelihood=True),  #  return log likelihood
)  # mp_ctx="forkserver")


#%%
summary_7 = model_7.summary()

# #%%
model_7.plot_posterior_predictive()
