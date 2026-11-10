# PreQ1-Riboswitch-Docking
Script to reproduce docking of PreQ1 and plastic additives to the PreQ1-III riboswitch

These scripts dock PreQ1, ten plastic additives, and two negative controls (PFOA and dibutyl phthalate) into the PreQ1-III riboswitch from Faecalibacterium prausnitzii (PDB 4RZD) and compare the poses with the crystallographic PreQ1.

## Files
- run_all.sh: converts each ligand SDF to PDBQT with Meeko and docks it with AutoDock Vina
- dock_one.sh: the same steps for a single compound
- results_table.py: reads every pose, calculates its distance to the crystal PreQ1 center, flags poses within 3 A as in the pocket, and calculates ligand efficiency; writes results_all_poses.csv
- hbond_one.py: compares ChimeraX hydrogen-bond output for one compound with PreQ1's six RNA contact atoms
- master_table.py: combines the pose table and hydrogen-bond files into master_results.csv
- offsite_compare.py: compares where off-pocket poses sit
- chimerax_commands.txt: the ChimeraX commands used for box center, receptor preparation, and hydrogen bonds

## Requirements
AutoDock Vina 1.2.7 (standalone executable, Apple Silicon build, in the working folder), Meeko 0.8.0, RDKit, Python 3.11, ChimeraX 1.12. The Python scripts use only the standard library.

## Inputs (not included)
- 4rzd_receptor.pdbqt: PDB 4RZD prepared in ChimeraX (Dock Prep) and converted with Meeko
- One PubChem 3D conformer SDF per compound, renamed to a short name: preq1, melamine, benzoguanamine, cyanuric_acid, mbi, bit, mbt, drometrizole, thiabendazole, hobt, benzothiazole, pfoa, dbp

## Run order (all in one folder)
    bash run_all.sh
    python3 results_table.py

Then, in ChimeraX, run the hbonds commands in chimerax_commands.txt for each compound, saving files named hbonds_<compound>.txt, then:

    python3 hbond_one.py <compound>
    python3 master_table.py
    python3 offsite_compare.py

## Settings
Search box centered at (32.05, 29.32, 107.90) A, size 20 x 20 x 20 A, exhaustiveness 32, random seed 42, 9 poses per ligand. Each compound was docked once.

## Notes
- The six PreQ1 contact atoms in hbond_one.py and master_table.py (C7 N3, C7 O2, U17 N3, U17 O4, U17 O2, A85 OP2) are hardcoded from the crystal structure.
- A pose counts as "in the pocket" if the center of its non-hydrogen atoms is within 3 A of the crystal PreQ1 center (a threshold chosen by the authors).
- Docking scores are predictions from a rigid-receptor model and are not binding free energies.
2. run_all.sh
#!/bin/bash
# Convert each ligand SDF to PDBQT and dock it with identical settings.
for name in preq1 melamine benzoguanamine cyanuric_acid mbi bit mbt drometrizole thiabendazole hobt benzothiazole pfoa dbp; do
  mk_prepare_ligand.py -i ${name}.sdf -o ${name}.pdbqt
  ./vina_1.2.7_mac_aarch64 --receptor 4rzd_receptor.pdbqt --ligand ${name}.pdbqt \
    --center_x 32.05 --center_y 29.32 --center_z 107.90 \
    --size_x 20 --size_y 20 --size_z 20 \
    --seed 42 --exhaustiveness 32 \
    --out ${name}_dock.pdbqt > ${name}_log.txt
done
3. dock_one.sh
#!/bin/bash
# Usage: bash dock_one.sh melamine
name=$1
mk_prepare_ligand.py -i ${name}.sdf -o ${name}.pdbqt
./vina_1.2.7_mac_aarch64 --receptor 4rzd_receptor.pdbqt --ligand ${name}.pdbqt \
  --center_x 32.05 --center_y 29.32 --center_z 107.90 \
  --size_x 20 --size_y 20 --size_z 20 \
  --seed 42 --exhaustiveness 32 \
  --out ${name}_dock.pdbqt | tee ${name}_log.txt
4. results_table.py
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
5. hbond_one.py
import sys, re
name = sys.argv[1]
ref_atoms = {("U",17,"N3"), ("C",7,"O2"), ("A",85,"OP2"),
             ("C",7,"N3"), ("U",17,"O4"), ("U",17,"O2")}
ref_res = {(r, n) for r, n, a in ref_atoms}
pat = re.compile(r"/A\s+(\S+)\s+(\d+)\s+(\S+)")
skip = {"PRF", "HOH", "UNL"}
contacts, nb = set(), 0
for line in open(f"hbonds_{name}.txt"):
    found = [(r, int(n), a) for r, n, a in pat.findall(line)
             if r not in skip and not a.startswith("H")]
    if found:
        nb += 1
        contacts.update(found)
res_set = {(r, n) for r, n, a in contacts}
print(name)
print(f"  H-bonds to RNA: {nb}")
print("  RNA atoms: " + ", ".join(f"{r}{n} {a}" for r, n, a in sorted(contacts, key=lambda x: x[1])))
print("  nucleotides: " + ", ".join(f"{r}{n}" for r, n in sorted(res_set, key=lambda x: x[1])))
print(f"  shared with PreQ1: {len(res_set & ref_res)}/3 nucleotides, {len(contacts & ref_atoms)}/6 atoms")
6. master_table.py
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
7. offsite_compare.py
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
8. chimerax_commands.txt
# ChimeraX commands used in this study (type into the ChimeraX command bar)

# --- Box center: centroid of crystallographic PreQ1 (residue PRF) ---
open 4rzd
select :PRF
define centroid sel radius 0.5 color orange
# Log reported 32.0547, 29.3151, 107.8961 (rounded to 32.05, 29.32, 107.90)

# --- Receptor: run Dock Prep from the Tools menu (Structure Editing > Dock Prep)
# with default options, then export (Meeko later dropped PRF) ---
select ~:PRF
save 4rzd_receptor.pdb format pdb models sel

# --- Hydrogen bonds of the crystal PreQ1 (reference contacts) ---
close
open 4rzd
hbonds :PRF restrict cross reveal true log true saveFile hbonds_crystal_preq1.txt

# --- Hydrogen bonds of one docked pose (example: melamine; repeat per compound) ---
close
open 4rzd
delete :PRF
delete solvent
open melamine_dock.pdbqt
close #2.2-9
hbonds #2.1 restrict cross reveal true log true saveFile hbonds_melamine.txt

# --- Residues within 5 A of crystallographic PreQ1 ---
close
open 4rzd
select :PRF :<5
info residues sel
9. The Meeko receptor conversion (also in chimerax_commands.txt, or add as its own file)

This step runs in Terminal, not ChimeraX, so it’s cleaner as a separate file named terminal_commands.txt:

# Terminal commands (in the Python 3.11 environment with Meeko installed)

# Receptor: convert to PDBQT. --allow_bad_res lets Meeko drop residues with no chemical
# template (this removed the crystal PreQ1, PRF; the 5'-terminal GTP was kept).
mk_prepare_receptor.py -i 4rzd_receptor.pdb -o 4rzd_receptor -p --allow_bad_res

# Ligand example (PubChem 3D conformer SDF renamed to a short name)
mk_prepare_ligand.py -i melamine.sdf -o melamine.pdbqt
