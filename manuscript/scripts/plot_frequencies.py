import json
import numpy as np
from datetime import datetime
from collections import defaultdict
from augur.utils import numeric_date
import matplotlib.pyplot as plt
from scipy.stats import scoreatpercentile

colors = {
"19B": {"color": "#FF8080", "ls":"-"},
"A": {"color": "#FF8080", "ls":"-"},
#
"19A": {"color": "#ffb247", "ls":"-"},
"B": {"color": "#ffb247", "ls":"-"},
"L": {"color": "#ffb247", "ls":"-"},
"O": {"color": "#FFD938", "ls":"-"},
"V": {"color": "#FFEEAA", "ls":"-"},
#
"20A": {"color": "#CDDE87", "ls":"-"},
"G": {"color": "#CDDE87", "ls":"-"},
"B.1": {"color": "#CDDE87", "ls":"-"},
"B.1.5": {"color": "#CDDE87", "ls":"--"},
"B.1.22": {"color": "#CCCCCC", "ls":"--"},
#
"20B": {"color": "#C6AFE9", "ls":"-"},
"B.1.1": {"color": "#C6AFE9", "ls":"-"},
"B.1.1.1": {"color": "#C6AFE9", "ls":"--"},
"GR": {"color": "#C6AFE9", "ls":"-"},
#
"20C": {"color": "#AAEEFF", "ls":"-"},
"GH": {"color": "#AAEEFF", "ls":"-"},
}


with open("data/ncov_europe.json", 'r') as fh:
    tree_json = json.load(fh)

def get_terminals(T):
    def _recurse(node, d):
        if ('children' in node) and len(node['children']):
            for c in node['children']:
                d = _recurse(c, d)
        else:
            d.append(node)

        return d

    return _recurse(T,[])

def get_frequencies(nodes, attribute, pivots, width_days=14):
    width_years = width_days/365
    freqs = defaultdict(lambda: np.zeros_like(pivots))
    for n in nodes:
        val = np.exp(-0.5*(pivots-n['node_attrs']['num_date']['value'])**2/width_years**2)
        freqs[n["node_attrs"][attribute]['value']] += val

    norm = np.sum(list(freqs.values()), axis=0)
    return {k:v/norm for k,v in freqs.items()}


def resample_countries(nodes_by_country):
    countries = list(nodes_by_country.keys())
    new_samples = []
    for ci in np.random.randint(0,len(countries), size=len(countries)):
        c = countries[ci]
        new_samples.extend(nodes_by_country[c])
    return new_samples


nodes = get_terminals(tree_json["tree"])

nodes_europe = [n for n in nodes if n['node_attrs']["region"]["value"]=="Europe"]

nodes_by_country = defaultdict(list)
for n in nodes_europe:
    nodes_by_country[n['node_attrs']["country"]["value"]].append(n)

start_date = datetime(2020,2,1).toordinal()
today = datetime.today().toordinal()
date_points = [datetime.fromordinal(x) for x in np.arange(start_date, today, 7)]
num_date_points = np.array([numeric_date(x) for x in date_points])


line_order = {'clade_membership':['19B', '19A', '20A', '20C', '20B'],
              'GISAID_clade':['S', 'L', 'O', 'V', 'G', 'GH', 'GR'],
              'pangolin_lineage':['A', 'B', 'B.1', 'B.1.5', 'B.1.22', 'B.1.1', 'B.1.1.1']}

panel_labels = {'clade_membership': 'Nextstrain',
              'GISAID_clade': 'GISAID',
              'pangolin_lineage': 'cov-lineages.org'}

n_samples=100
fs=12
fig, axs = plt.subplots(3, 1, figsize = (10,10), sharex=True)
for ax, attr in zip(axs, ["clade_membership", "GISAID_clade", "pangolin_lineage"]):
    frequencies = get_frequencies(nodes_europe, attr,
                                  num_date_points, width_days=10)

    resampled_frequencies = [get_frequencies(resample_countries(nodes_by_country), attr,
                                             num_date_points, width_days=10)
                             for i in range(n_samples)]


    for ci, clade in enumerate(line_order[attr]):
    # for ci, clade in enumerate(sorted(frequencies.keys())):
        if frequencies[clade].max() > 0.05:
            col = colors[clade]["color"] if clade in colors else f'C{ci%10}'
            ls = colors[clade]["ls"] if clade in colors else f'-'
            tmp_samples = np.array([(f[clade] if clade in f else np.zeros_like(num_date_points))
                                    for f in resampled_frequencies])
            upper = scoreatpercentile(tmp_samples, 75, axis=0)
            lower = scoreatpercentile(tmp_samples, 25, axis=0)
            ax.plot(date_points, frequencies[clade], label=clade, lw=3, c=col, ls=ls)
            ax.fill_between(date_points, lower, upper, color=col, alpha=0.3)

    ax.text(date_points[1], 0.7, panel_labels[attr], fontsize=fs*1.4)
    ax.legend(fontsize=fs, ncol=3, loc=9)
    ax.set_ylabel("frequency", fontsize=fs)
    ax.set_ylim([0,.8])
    ax.tick_params(labelsize=fs)
    ax.tick_params('x', rotation=30)

plt.tight_layout()
plt.savefig(f"figures/clade_frequencies.pdf")
