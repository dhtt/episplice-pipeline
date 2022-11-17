# %%
import pandas as pd
import os.path
import matplotlib.pylab as plt
import numpy as np
import seaborn as sb
import scipy.stats as st

path_to_DEXSeq = "/home/dhthutrang/Krebs/Reddy/data/mRNA_seq/DEXSeq_output"
df_exp = pd.read_csv(os.path.join(path_to_DEXSeq, "csv", "MCF7_DMSO_MCF7_50nM_res.csv"), sep='\t', header=0)
print(df_exp)
# %%

def get_best_distribution(data):
    dist_names = [ "exponweib", "weibull_max", "weibull_min", "pareto", "genextreme"]
    dist_results = []
    params = {}
    for dist_name in dist_names:
        print(dist_name)
        dist = getattr(st, dist_name)
        print(dist)
        param = dist.fit(data)

        # params[dist_name] = param
        # Applying the Kolmogorov-Smirnov test
    #     D, p = st.kstest(data, dist_name, args=param)
    #     print("p value for " + dist_name + " = " + str(p))
    #     dist_results.append((dist_name, p))

    # # select the best fitted distribution
    # best_dist, best_p = (max(dist_results, key=lambda item: item[1]))
    # # store the name of the best fit and its p value

    # print("Best fitting distribution: "+str(best_dist))
    # print("Best p value: "+ str(best_p))
    # print("Parameters for the best fit: "+ str(params[best_dist]))

    # return best_dist, best_p, params[best_dist]
# %%
np.isfinite(df_exp['stat'].values).all()
df_exp['stat'].replace([np.isfinite, -np.isfinite], NaN, )
# get_best_distribution(df_exp['stat'].values)

# %%
plt.figure()
sb.histplot(df_exp['stat'], bins=100, kde=True, log_scale=True)
plt.show()

# %%
# df_exp['stat'].describe()
kstest(df_exp['stat'], 'norm')

# %%
