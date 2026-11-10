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
