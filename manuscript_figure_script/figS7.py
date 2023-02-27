import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator, FixedLocator





hey = "/Users/yeongseon/Dropbox/PhD_Emory_university/1_projects/Project_seg_site/simulated_data/simpleSEIR_R0=16e-1_mu2e-1/seed1234_prop500_win4/seed230201_showparticles_R0,timestart/seed202302.tsv"
df = pd.read_csv(hey, sep = "\t")

plt.rcParams.update({'font.size': 10})
fig, axes = plt.subplots(1, 1, figsize=[13.5, 5], dpi=300)

color = {'X': 'grey', 'V': 'k'}
alpha = {'X': 0.1, 'V': 0.3}
for i in range(len(df)):#range(800):#range(len(df)):
    #if df.iloc[i]['selected'] == "V":
    axes.plot([df.iloc[i]['t_before'], df.iloc[i]['t_after']], [df.iloc[i]['s_before'], df.iloc[i]['s_after']],
              c=color[df.iloc[i]['selected']],
              alpha=alpha[df.iloc[i]['selected']],
              linewidth = 1)

    #
    # if i < 600:
    #     if df.iloc[i]['selected'] == "V":
    #         axes.plot([df.iloc[i]['t_before'], df.iloc[i]['t_after']], [df.iloc[i]['s_before'], df.iloc[i]['s_after']],
    #                   c=color[df.iloc[i]['selected']],
    #                   alpha=alpha[df.iloc[i]['selected']])
    #
    # else:
    #     axes.plot([df.iloc[i]['t_before'], df.iloc[i]['t_after']],[df.iloc[i]['s_before'], df.iloc[i]['s_after']],
    #               c = color[df.iloc[i]['selected']],
    #               alpha = alpha[df.iloc[i]['selected']])




fig.show()