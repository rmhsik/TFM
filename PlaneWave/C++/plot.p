set terminal png size 1000,1000
set output "initial.png"
plot "InitialPhi.dat" matrix w image
set output "final.png"
plot "FinalPhi.dat" matrix w image
