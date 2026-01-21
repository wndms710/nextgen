#!/bin/bash

if [ "$#" -lt 1 ]; then
  echo "Usage: ./plot_gain.sh <resonance-log-file>"
  exit 1
fi

gnuplot -persist <<-EOFMarker
  set logscale y
  set xlabel "Freq"
  set ylabel "Gain [nd]"
  plot "< cat $1 | grep XXXX" using 2:3 notitle w lp
EOFMarker
