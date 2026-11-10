import csv, itertools
def center(f, pose):
    p = 0; xs, ys, zs = [], [], []
    for line in open(f):
        if line.startswith("MODEL"):
            p = int(line.split()[1]); xs, ys, zs = [], [], []
        elif line.startswith(("ATOM", "HETATM")) and p == pose:
            if line.split()[-1] in ("H", "HD"):
                continue
            xs.append(float(line[30:38])); ys.append(float(line[38:46])); zs.append(float(line[46:54]))
        elif line.startswith("ENDMDL") and p == pose:
            return (sum(xs)/len(xs), sum(ys)/len(ys), sum(zs)/len(zs))
items = [("pfoa", 1), ("dbp", 1), ("drometrizole", 1), ("thiabendazole", 1)]
for r in csv.DictReader(open("results_all_poses.csv")):
    if r["compound"] == "preq1" and r["in_pocket"] == "no":
        items.append(("preq1", int(r["pose"])))
cs = {f"{n}_pose{p}": center(f"{n}_dock.pdbqt", p) for n, p in items}
for k, c in cs.items():
    print(k, tuple(round(x, 1) for x in c))
print("\nPairwise center distances (A):")
for (a, ca), (b, cb) in itertools.combinations(cs.items(), 2):
    d = sum((x-y)**2 for x, y in zip(ca, cb)) ** 0.5
    print(f"{a} vs {b}: {d:.1f}")
