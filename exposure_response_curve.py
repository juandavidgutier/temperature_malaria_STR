import numpy as np
import pandas as pd
from scipy.stats import expon
import scipy.stats as stats
from causal_curve import TMLE_Regressor
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from plotnine import ggplot, aes, geom_line, geom_ribbon, ggtitle, labs, ylim, xlim

# Data
top100 = pd.read_csv("D:/data07-23.csv", encoding='latin-1') 

# Confounders as binary
def transform_binary(df):
    for col in df.columns[8:20]:
        median = df[col].median()  
        df[col] = (df[col] > median).astype(int)  
    return df

top100 = transform_binary(top100)

# Temperature between 15-30°C
data_temp = top100[top100['Temperature'].between(15, 30)]

# Eliminate NaNs
top100_clean = data_temp.dropna()

# TMLE
tmle_temp = TMLE_Regressor(n_estimators=2500, random_seed=123, bandwidth=999, verbose=True, 
                           max_depth=4, learning_rate=0.001)  
X_data = top100_clean.iloc[:, 6:20]
T_data = top100_clean['Temperature']
y_data = top100_clean['sir']
     
# Fit
tmle_temp.fit(T=T_data, X=X_data, y=y_data)
    
# CI
tmle_results_temp = tmle_temp.calculate_CDRC(0.95)
    
# Figure 2
plot = (
ggplot(aes(x=tmle_results_temp.Treatment, y=tmle_results_temp.Causal_Dose_Response)) 
+ geom_line() 
+ geom_ribbon(aes(ymin=tmle_results_temp.Lower_CI, ymax=tmle_results_temp.Upper_CI), alpha=0.1)
+ labs(x='Temperature', y='SIR') 
    ) 
print(plot)