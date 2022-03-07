#!/bin/bash -e

echo "      date     time $(free -hm | grep total | sed -E 's/^    (.*)/\1/g')" > "memory_log.txt"
while true; do
    echo "$(date '+%Y-%m-%d %H:%M:%S') $(free -hm | grep Mem: | sed 's/Mem://g')" >> "memory_log.txt"
    sleep 30
done
