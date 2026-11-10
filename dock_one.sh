#!/bin/bash
# Usage: bash dock_one.sh melamine
name=$1
mk_prepare_ligand.py -i ${name}.sdf -o ${name}.pdbqt
./vina_1.2.7_mac_aarch64 --receptor 4rzd_receptor.pdbqt --ligand ${name}.pdbqt \
  --center_x 32.05 --center_y 29.32 --center_z 107.90 \
  --size_x 20 --size_y 20 --size_z 20 \
  --seed 42 --exhaustiveness 32 \
  --out ${name}_dock.pdbqt | tee ${name}_log.txt
