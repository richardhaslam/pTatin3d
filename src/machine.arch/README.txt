To add a new configuration here, note the following naming scheme.

For files aimed at operating systems, use
  machine.<os>.<compiler>.[<modifier>.]arch
for example
  machine.linux.gcc.arch
  machine.osx.gcc.externals.arch

For files aimed at specific machines, use
  <machine name>.<os>.<compiler>.[<modifier>.]arch
for example
  daint.crayxc30.craycc.arch
