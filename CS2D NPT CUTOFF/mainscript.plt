path='C:\Users\medos\programiranje\urbic\programi\2022\CS2D NPT CUTOFF'
nseries=20

pres_points=1
pres_start=1.000
pres_step=0.
temp_points=1
temp_start=0.500
temp_step=0.
do for [j=1:pres_points] {
  pres=pres_start+(j-1)*pres_step
  do for [k=1:temp_points] {
    temp=temp_start+(k-1)*temp_step
    load 'snapshot.plt'
    load 'min.plt'
    load 'ekv.plt'
    load 'sam.plt'
    load 'corr.plt'
    load 'slike.plt'
    load 'movie.plt'
  }
}
