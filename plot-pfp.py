import json
import matplotlib.pyplot as plt
from sys import argv

data = json.loads(input())
hist = {}
for d in data:
    w, p, total = d["w"], d["p"], d["total"]
    if p not in hist:
        hist[p] = {}
    hist[p][w] = [c / total for c in d["hist"]]

np = len(hist)
fig, axes = plt.subplots(np, 1, layout="constrained", figsize=(7, 3 + np))
plt.suptitle("PFP phrase size distribution")
for i, p in enumerate(hist):
    ax = axes[i] if np > 1 else axes
    nw = len(hist[p])
    alpha = 0.8 if nw > 1 else 1
    for j, w in enumerate(hist[p]):
        h = hist[p][w]
        x = [k + j / nw for k in range(len(h))]
        ax.bar(x, h, width=1 / nw, alpha=alpha, label=f"$w={w}$")
    x = list(range(1, len(h)))
    y = [1 / p * (1 - 1 / p) ** (k - 1) for k in x]
    ax.plot(x, y, "r", label="geometric")
    ax.set_title(f"$p={p}$", y=1.0, pad=-15)
    ax.set_ylim(top=1.25 / p)
    if i == 0:
        ax.legend()
if len(argv) > 1:
    plt.savefig(argv[1], bbox_inches="tight", dpi=300)
else:
    plt.show()
