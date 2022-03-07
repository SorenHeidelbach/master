#!/bin/bash -e

mpstat | awk '$3 ~ /CPU/ {print}' > "cpu_log.txt"
while true; do
    mpstat | awk '$3 ~ /all/ {print}' >> "cpu_log.txt"
    sleep 30
done
