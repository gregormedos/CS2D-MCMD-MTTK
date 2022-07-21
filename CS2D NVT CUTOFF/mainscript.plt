path='C:\Users\medos\simulacije\CS2DMD NVT CUTOFF izoterme'
N=100
nseries=20

temp_points=2
temp_start=0.500
temp_step=0.500
dens_points=2
dens_start=0.400
dens_step=0.100
do for [j=1:temp_points] {
  temp=temp_start+(j-1)*temp_step
  do for [k=1:dens_points] {
    dens=dens_start+(k-1)*dens_step
    load 'min.plt'
    load 'ekv.plt'
    load 'sam.plt'
    load 'snapshot.plt'
    load 'slike.plt'
    load 'movie.plt'
    load 'corr.plt'
    load 'msd.plt'
  }
}
