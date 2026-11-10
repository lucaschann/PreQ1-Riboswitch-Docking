import csv, os, re
ref_atoms = {("U",17,"N3"), ("C",7,"O2"), ("A",85,"OP2"),
             ("C",7,"N3"), ("U",17,"O4"), ("U",17,"O2")}
ref_res = {(r, n) for r, n, a in ref_atoms}
pat = re.compile(r"/A\s+(\S+)\s+(\d+)\s+(\S+)")
skip = {"PRF", "HOH", "UNL"}

def hbonds(name):
    f = f"hbonds_{name}.txt"
    if not os.path.exists(f):
        return None
    contacts, nb = set(), 0
    for line in open(f):
        found = [(r, int(n), a) for r, n, a in pat.findall(line)
                 if r not in skip and not a.startswith("H")]
        if found:
            nb += 1
            contacts.update(found)
    res = {(r, n) for r, n, a in contacts}
    atoms = "; ".join(f"{r}{n} {a}" for r, n, a in sorted(contacts, key=lambda x: x[1]))
    return nb, len(contacts & ref_atoms), len(res & ref_res), atoms

poses = {}
for row in csv.DictReader(open("results_all_poses.csv")):
    poses.setdefault(row["compound"], []).append(row)

order = {"reference": 0, "candidate": 1, "control": 2}
out = []
for name, rows in poses.items():
    top = rows[0]
    inp = [r for r in rows if r["in_pocket"] == "yes"]
    best_in = min(float(r["score_kcal_mol"]) for r in inp) if inp else "none"
    hb = hbonds(name)
    out.append([name, top["role"], float(top["score_kcal_mol"]),
                top["ligand_efficiency"], float(top["dist_from_PreQ1_site_A"]),
                f"{len(inp)}/{len(rows)}", best_in] +
               (list(hb) if hb else ["n/a", "n/a", "n/a", "n/a"]))
out.sort(key=lambda r: (order.get(r[1], 3), r[2]))

header = ["compound", "role", "top_score", "ligand_efficiency", "top_pose_dist_A",
          "poses_in_pocket", "best_in_pocket_score", "hbonds_to_RNA",
          "shared_atoms_of_6", "shared_nt_of_3", "RNA_atoms_contacted"]
with open("master_results.csv", "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(header)
    w.writerows(out)
for r in out:
    print(r[:10])
print("\nWrote master_results.csv")
