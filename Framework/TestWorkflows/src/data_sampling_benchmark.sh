#!/usr/bin/env bash

#set -e ;# exit on error
set -u ;# exit when using undeclared variable
#set -x ;# debugging


PAYLOAD_SIZE=(256 4096 65536 1048576 16777216 268435456);
NB_PRODUCERS=(1 4 16 64);
USLEEP_TIME=(1000000 10000 100 0);
#NB_DISPATCHERS=(1 2 4 8);
NB_DISPATCHERS=(1);

echo "payload size        , nb producers        , usleep time         , nb dispatchers      , messages per second" > data-sampling-benchmark

for payload_size in ${PAYLOAD_SIZE[@]}; do
  for nb_producers in ${NB_PRODUCERS[@]}; do
    for usleep_time in ${USLEEP_TIME[@]}; do
      for nb_dispatchers in ${NB_DISPATCHERS[@]}; do
        for run in {1..5}
        do
          echo "***************************
          Launching test for payload size $payload_size bytes, $nb_producers producers, $usleep_time us of their sleeptime, $nb_dispatchers dispatchers"

          printf "%20s," "$payload_size" >> data-sampling-benchmark
          printf "%21s," "$nb_producers" >> data-sampling-benchmark
          printf "%21s," "$usleep_time" >> data-sampling-benchmark
          printf "%21s," "$nb_dispatchers" >> data-sampling-benchmark

          timeout -k 60s 8m o2-testworkflows-datasampling-benchmark -q -b --payload-size $payload_size --producers $nb_producers --dispatchers $nb_dispatchers --usleep $usleep_time --shm-segment-size 300000000
          pkill -9 -f o2-testworkflows-datasampling-benchmark
          printf "\n" >> data-sampling-benchmark
        done
      done
    done
  done
done
