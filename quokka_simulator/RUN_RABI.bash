#!/bin/bash

# p (physical error rate) values:
p_values=(0 0.0005 0.001 0.002 0.003 0.004)


for p in "${p_values[@]}"; do


	echo "Running p = $p rabi simulations..."

	for i in $(seq 0 10); do
   	 ./rabi_simulator_modules/d3rabi_v5.out 3 0.25 0 $p 2 100000 &
	done

	wait

	echo "Completed"

done

echo "Done"
