import json
import numpy as np
from datetime import datetime
from collections import defaultdict
from augur.utils import numeric_date
import matplotlib.pyplot as plt

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

nodes = get_terminals(tree_json["tree"])

nodes_europe = [n for n in nodes if n['node_attrs']["region"]["value"]=="Europe"]


start_date = datetime(2020,2,1).toordinal()
today = datetime.today().toordinal()
date_points = [datetime.fromordinal(x) for x in np.arange(start_date, today, 7)]
num_date_points = np.array([numeric_date(x) for x in date_points])



for attr in ["clade_membership", "gisaid_clade_membership"]:
    frequencies = get_frequencies(nodes_europe, attr,
                                  num_date_points, width_days=10)

    fs=16
    plt.figure(figsize = (12,6))
    for clade in sorted(frequencies.keys()):
        plt.plot(date_points, frequencies[clade], label=clade, lw=3)

    plt.legend(fontsize=fs)
    plt.ylabel("frequency in Europe", fontsize=fs)
    plt.ylim([0,1])
    plt.tick_params(labelsize=fs)
    plt.tick_params('x', rotation=30)
    plt.tight_layout()

    plt.savefig(f"figures/{attr}.png")