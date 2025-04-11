#! /bin/bash

# 1. run disc
java -Xmx100g \
    -jar TIGER.jar \
    -app FastCall2 \
    -mod disc \
    -a chr001.fa \
    -b taxaBamMap.txt \
    -c 0 \
    -d 30 \
    -e 20 \
    -f 2 \
    -g 0.2 \
    -h 3 \
    -i 0.8 \
    -j 0.35 \
    -k 0.2 \
    -l 1 \
    -m 32 \
    -n /ing \
    -o /usr/local/bin/samtools > log.txt &


java -Xmx100g -jar TIGER.jar -app FastCall2 -mod disc -a chr001.fa -b taxaBamMap.txt -c 0 -d 30 -e 20 -f 2 -g 0.2 -h 3 -i 0.8 -j 0.35 -k 0.2 -l 1 -m 32 -n /ing -o /usr/local/bin/samtools > log.txt &

# 2. run blib
java -Xmx100g -jar TIGER.jar -app FastCall2 -mod blib -a chr001.fa -b 1 -c 2 -d 32 -e /ing -f /vLib > log.txt &