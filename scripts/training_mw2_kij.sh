
for m in 1E0 1E-1 1E-2 1E-3 1E-4 1E-5 ; do
  for w in 1 10 50 100 ; do
    for k in 1E3 1E4 2E4; do
#      for n in 1 2 3 4 5 6 7 8 9 0; do
        python ./src/activesolid/harmonic_train_bench.py --lr 1E9 --V_seed 0 --m $m --k $k --omega $w --xi 0.5
        python ./src/activesolid/harmonic_train_bench.py --lr 1E9 --V_seed 0 --m $m --k $k --omega $w --xi 1.0
        python ./src/activesolid/harmonic_train_bench.py --lr 1E9 --V_seed 0 --m $m --k $k --omega $w --xi 1.5
#      done
    done
  done
done