import hssm
import arviz as az
import matplotlib.pyplot as plt
import pandas as pd

# Set float type to float32 to avoid a current bug in PyMC
hssm.set_floatX("float32")

# Load a package-supplied dataset
##cav_data = hssm.load_data('cavanagh_theta')
##cav_data['rt'] = cav_data['rt'][cav_data['response']==0].apply(lambda x:x*-1)
frg_data = pd.read_csv("frg_norm_hssm.csv")

# Define a basic hierarchical model with trial-level covariates
model = hssm.HSSM(
    model="ddm",
    data=frg_data,
    include=[
        {
            "name": "v",
            "prior": {
                "Intercept": {"name": "Normal", "mu": 0.0, "sigma": 1.0},
                "reward": {"name": "Normal", "mu": 0.0, "sigma": 1.0},
            },
            "formula": "v ~ (1|subj_idx) + reward",
            "link": "identity",
        },
    ],
)

g = model.graph()
g.view()

# Sample from the posterior for this model
model.sample(
    sampler="mcmc",  # type of sampler to choose, 'nuts_numpyro', 'nuts_blackjax' of default pymc nuts sampler "mcmc"
    cores=10,  # how many cores to use
    chains=5,  # how many chains to run
    draws=500,  # number of draws from the markov chain
    tune=500,  # number of burn-in samples
    idata_kwargs=dict(log_likelihood=True)
    )

model.sample_posterior_predictive()
az.plot_ppc(model.traces,show=True)
