sequence = "daefrhdsgyevhhqklvffaedvgsnkgaiiglmvggvvia" 
for aa in sequence: cmd._alt(aa)
alter (all),resi=str(int(resi)-1)
alter (all), chain="A"
save surf.pdb
quit 
