#set term gif animate delay 4
#set output "test.gif"

set size square
set xrange[0:30e-6]
set yrange[-15e-6:15e-6]
#set zrange[-10e-6:10e-6]
#set view 90,0

head_e="./outputs/electrons"
head_h="./outputs/holes"

do for [i = 0:2000 ] {
    fname_e=head_e.i.".txt"
    fname_h=head_h.i.".txt"
    set label 1 at graph 0.1,0.1 "time (fs)=".i*10
#    splot fname_e with points pt 7 lt 3,fname_h with points pt 6 lt 1
    plot fname_e with points pt 7 lt 3,fname_h with points pt 6 lt 1

pause 0.1
}

