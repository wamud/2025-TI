echo "Running p = 0 rabi simulations..."

for i in $(seq 0 3); do
    ./rabi_simulator_modules/d3rabi_v4 3 0.25 0 0 2 100 > /dev/null &
done

wait

echo "Completed"

echo "Performing p = 0.001 rabi simulations..."

for i in $(seq 0 3); do
    ./rabi_simulator_modules/d3rabi_v4 3 0.25 0 0.001 2 100 > /dev/null &
done

wait

echo "Completed"

echo "Performing p = 0.002 rabi simulations..."

for i in $(seq 0 3); do
    ./rabi_simulator_modules/d3rabi_v4 3 0.25 0 0.002 2 100 > /dev/null &
done

wait


echo "Done"