import json
import matplotlib.pyplot as plt
from sys import argv

if len(argv) > 1:
    filename = argv[1]
    with open(filename, "r") as f:
        data = json.load(f)
else:
    data = json.loads(input())

hist = {}
for d in data:
    w, p, total = d["w"], d["p"], d["total"]
    if p not in hist:
        hist[p] = {}
    hist[p][w] = [c / total for c in d["hist"]]

for p in hist:
    nw = len(hist[p])
    alpha = 0.8 if nw > 1 else 1
    for i, w in enumerate(hist[p]):
        h = hist[p][w]
        x = [k + i / nw for k in range(len(h))]
        plt.bar(x, h, width=1 / nw, alpha=alpha, label=f"$w={w}$")
    x = list(range(1, len(h)))
    y = [1 / p * (1 - 1 / p) ** (k - 1) for k in x]
    plt.plot(x, y, "r", label="geometric")
    plt.ylim(top=2 / p)
    plt.title(f"PFP phrase size distribution ($p={p}$)")
    plt.legend()
    plt.show()
