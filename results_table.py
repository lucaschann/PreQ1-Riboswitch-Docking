import glob, csv
ref = (32.05, 29.32, 107.90)
cutoff = 3.0
roles = {"preq1": "reference", "pfoa": "control", "dbp": "control"}
header = ["compound", "role", "pose", "score_kcal_mol", "rmsd_lb", "rmsd_ub",
          "dist_from_PreQ1_site_A", "in_pocket", "heavy_atoms", "ligand_efficiency"]
rows = []
for f in sorted(glob.glob("*_dock.pdbqt")):
    name = f.replace("_dock.pdbqt", "")
    pose = 0
    xs, ys, zs = [], [], []
    for line in open(f):
        if line.startswith("MODEL"):
            pose += 1
            xs, ys, zs = [], [], []
        elif line.startswith("REMARK VINA RESULT"):
            p = line.split()
            score, lb, ub = float(p[3]), float(p[4]), float(p[5])
        elif line.startswith(("ATOM", "HETATM")):
            if line.split()[-1] in ("H", "HD"):
                continue
            xs.append(float(line[30:38]))
            ys.append(float(line[38:46]))
            zs.append(float(line[46:54]))
        elif line.startswith("ENDMDL"):
            c = (sum(xs)/len(xs), sum(ys)/len(ys), sum(zs)/len(zs))
            d = sum((a-b)**2 for a, b in zip(c, ref)) ** 0.5
            heavy = len(xs)
            rows.append([name, roles.get(name, "candidate"), pose, score, lb, ub,
                         round(d, 1), "yes" if d <= cutoff else "no",
                         heavy, round(-score/heavy, 2)])
with open("results_all_poses.csv", "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(header)
    w.writerows(rows)
print("compound, role, top_score, top_pose_dist, poses_in_pocket, best_in_pocket_score")
seen = []
for r in rows:
    if r[0] not in seen:
        seen.append(r[0])
for n in seen:
    rs = [r for r in rows if r[0] == n]
    inp = [r for r in rs if r[7] == "yes"]
    best = min(r[3] for r in inp) if inp else "none"
    print(f"{n}, {rs[0][1]}, {rs[0][3]}, {rs[0][6]}, {len(inp)}/{len(rs)}, {best}")
