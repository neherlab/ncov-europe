
# Run this from inside of `ncov` repo.
# First run the clades algorithm on whole datset:
#     python assign_clades.py --alignment results/aligned.fasta --output clades_2002-06-26.tsv --chunk-size 100

import pandas as pd 
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick

# We have no sequences from Uzbekistan currently

metadata = pd.read_csv("data/metadata.tsv", delimiter='\t')

# Rename Russia
metadata.loc[metadata.country=='Russia', "country"] = 'Russian Federation'

# Reassign those who should be in Europe but are in Asia (or wherever)
metadata.loc[metadata.country=='Israel', "region"] = 'Europe'
metadata.loc[metadata.country=='Kazakhstan', "region"] = 'Europe'
metadata.loc[metadata.country=='Georgia', "region"] = 'Europe'
metadata.loc[metadata.country=='Uzbekistan', "region"] = 'Europe' #there aren't any yet

# Get just Europe
metadata_europe = metadata[metadata.region == 'Europe']

clades_file = pd.read_csv("clades_2002-06-14a.tsv", delimiter='\t')

meta_clades = metadata_europe.merge(clades_file, left_on='strain', right_on='name')

# Assign new clade based on crossover of GISAID & Nextstrain
meta_clades['new_clades']=''
# Those that are Nextstrain 20A *and* GISAID GH are a distinct set, not perfectly matching up
meta_clades.loc[(meta_clades.GISAID_clade=='GH') & (meta_clades.clade=='20A'), "new_clades"] = '20A & GH (Discordant)'
# Those remaining in 20A are pretty much G
meta_clades.loc[(meta_clades.clade=='20A') & (meta_clades.new_clades==''), "new_clades"] = "20A | G"
# Those in 20C are GH (though GH is larger than just these)
meta_clades.loc[(meta_clades.clade=='20C') & (meta_clades.new_clades==''), "new_clades"] = "20C | GH"

#The last three are direct match-ups
meta_clades.loc[(meta_clades.clade=='20B') & (meta_clades.new_clades==''), "new_clades"] = "20B | GR"
meta_clades.loc[(meta_clades.clade=='19A') & (meta_clades.new_clades==''), "new_clades"] = "19A | L/V/O"
meta_clades.loc[(meta_clades.clade=='19B') & (meta_clades.new_clades==''), "new_clades"] = "19B | S"

# Make a table giving counts by country and clade
#from this: https://stackoverflow.com/questions/44298779/pivot-table-with-group-and-without-value-field

#clade_table = meta_clades.groupby(['country','clade']).size().unstack(fill_value=0)
clade_table = meta_clades.groupby(['country','new_clades']).size().unstack(fill_value=0)

n_country = clade_table.sum(axis=1)

#clades = ['19A', '19B', '20A', '20B', '20C']

#clades = ['19A | L/V/O', '19B | S', '20A | G', '20B | GR', '20C | GH', 'Discordant: 20A & GH']
clades = ['19A | L/V/O', '19B | S', '20A | G', '20A & GH (Discordant)', '20C | GH', '20B | GR']


#set to percents
#from this: https://stackoverflow.com/questions/42006346/pandas-convert-columns-to-percentages-of-the-totals/42006745

clade_table = clade_table[clades].div(clade_table[clades].sum(axis=1), axis=0).multiply(100)

bar_colors = ['C0', 'C1', 'C2', 'C7', 'C3', 'C4']

# Plot
# from this: https://pstblog.com/2016/10/04/stacked-charts
fs=12
rects = clade_table.loc[:,clades].plot.bar(stacked=True, figsize=(10,7),
    color=bar_colors)
ax = plt.gca()
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=fs)
for label in plt.gca().get_xticklabels():
    label.set_horizontalalignment('right')
plt.tick_params('x', rotation=30, labelsize=0.9*fs)
ax.yaxis.set_major_formatter(mtick.PercentFormatter())
plt.ylim(0,100)
plt.xlabel('')

for i,c in enumerate(clade_table.index):
    ax.text(i-0.5,100, f"n={n_country[i]}", rotation=30, fontsize=0.9*fs)
# for p in rects.patches:
#     rects.annotate(str(p.get_height()), xy=(p.get_x(), p.get_height()))

plt.tight_layout()
plt.show()



###### Alternative version with matching coloring to Fig 1

clades = ['19B | S', 
          '19A | L/V/O',
          '20A | G', 
          '20A & GH (Discordant)', 
          '20C | GH', 
          '20B | GR']

clades = ['20B | GR',
          '20C | GH', 
          '20A & GH (Discordant)', 
          '20A | G', 
          '19A | L/V/O',
          '19B | S'
          ]

#set to percents
#from this: https://stackoverflow.com/questions/42006346/pandas-convert-columns-to-percentages-of-the-totals/42006745

clade_table = clade_table[clades].div(clade_table[clades].sum(axis=1), axis=0).multiply(100)

bar_colors = ['#ff8080', '#ffb247', '#cdde87', '#999999', '#aaeeff', '#c6afe9']

# Plot
# from this: https://pstblog.com/2016/10/04/stacked-charts
fs=12
rects = clade_table.loc[:,clades].plot.bar(stacked=True, figsize=(10,7),
    color=bar_colors)
ax = plt.gca()
handles, labels = ax.get_legend_handles_labels()
plt.legend(handles[::-1], labels[::-1], loc='center left', bbox_to_anchor=(1, 0.5), fontsize=fs)
for label in plt.gca().get_xticklabels():
    label.set_horizontalalignment('right')
plt.tick_params('x', rotation=30, labelsize=0.9*fs)
ax.yaxis.set_major_formatter(mtick.PercentFormatter())
plt.ylim(0,100)
plt.xlabel('')

for i,c in enumerate(clade_table.index):
    ax.text(i-0.5,100, f"n={n_country[i]}", rotation=30, fontsize=0.9*fs)
# for p in rects.patches:
#     rects.annotate(str(p.get_height()), xy=(p.get_x(), p.get_height()))

plt.tight_layout()
plt.show()



#########################################################
#########################################################
#########################################################
#########################################################
#########################################################
#########################################################
#########################################################
#########################################################
#########################################################

## Random stuff I tried - bad

countries = list(meta_clades.country.unique())

pd.DataFrame(index=countries, columns=clades)

for cntry in countries:
    mini_meta = meta_clades[meta_clades.country == cntry]



meta_clades.pivot(index='country', columns='clade')

meta_clades.clades.unique()