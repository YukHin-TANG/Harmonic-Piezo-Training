#for n in 1 2 3 4 5 6 7 8 9 0; do
#    for angle in 0 15 30 45 60 75 90 105 120 135 150 165; do
#      for A in 1E-8 1E-7 1E-6 1E-5 1E-4 1E-3 1E-2 1E-1; do
#        python ./src/activesolid/harmonic_train_bench.py --lr 1E9 --V_seed $n --A $A --angle $(( 180 + angle))
#        python ./src/activesolid/harmonic_train_bench.py --lr 1E9 --V_seed $n --A $A --angle $angle
#    done
#  done
#done
# 100 runs
for i in {1..100}; do
  echo $i
  python ./src/activesolid/harmonic_train_bench.py --lr 1E9 --V_seed 1 --angle $((RANDOM % 360)) --A $((RANDOM % 10))E$((-RANDOM % 7))
done