# Gnuplot script file for plotting data
set encoding utf8
set autoscale                                                       # scale axes automatically
unset logscale                                                      # remove any log-scaling
unset xlabel                                                        # remove any previous xlabel
unset ylabel                                                        # remove any previous ylabel
set xtics auto                                                      # set xtics automatically
set ytics auto                                                      # set ytics automatically
unset key
unset size

cd sprintf('%s\%.3f\%.3f', path, pres, temp)
set terminal postscript eps enhanced color
set size 0.5,0.5
set size square
unset xtics
unset ytics
halflbox=0.5
set xrange [-halflbox:halflbox]
set yrange [-halflbox:halflbox]
set output 'slika_sta.eps'
plot 'snapshot_sta' with circles lc rgb "black"
set output 'slika_min.eps'
plot 'snapshot_min' with circles lc rgb "black"
set output 'slika_ekc.eps'
plot 'snapshot_ekc' with circles lc rgb "black"
set output 'slika_int.eps'
plot 'snapshot_int' with circles lc rgb "black"
set output 'slika_ekd.eps'
plot 'snapshot_ekd' with circles lc rgb "black"

set size nosquare
cd path