#!/bin/bash

# p (physical error rate) values:
p_values=(0)


for p in "${p_values[@]}"; do


	echo "Running p = $p rabi simulations..."

	for i in $(seq 0 10); do
   	 ./rabi_simulator_modules/d3rabi_v4.out 3 0.25 0 $p 2 1000000 &
	done

	wait

	echo "Completed"

done

echo "Done"
