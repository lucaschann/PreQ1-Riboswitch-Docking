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
