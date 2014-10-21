from MDAnalysis import Universe, Writer
u = Universe("protein.gro", "protein-short.xtc")
w = Writer("protein-short.dcd", u.trajectory.numatoms)
for ts in u.trajectory:
    w.write(ts)
w.close_trajectory()
