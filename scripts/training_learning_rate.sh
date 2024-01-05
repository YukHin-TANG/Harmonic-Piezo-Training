
for n in 1 2 3 4 5 6 7 8 9 0; do
  for rate in 1E7 1E8 1E9 1E10; do
      python ./src/activesolid/harmonic_train_bench.py --lr $rate --V_seed $n
  done
done
