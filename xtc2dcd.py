from MDAnalysis import Universe, Writer
u = Universe("md.part0001.gro", "md.part0001.xtc")
w = Writer("part1.dcd", u.trajectory.numatoms)
for ts in u.trajectory:
    w.write(ts)
w.close_trajectory()
