#!/bin/bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )"

# for p in 32 #16
# do
#     for i in 1 5 10
#     do
#         $DIR/cmake-build-debug/pairhmm $p 8 8 $i
#         $DIR/cmake-build-debug/pairhmm $p 16 16 $i
#         $DIR/cmake-build-debug/pairhmm $p 24 24 $i
#         $DIR/cmake-build-debug/pairhmm $p 32 32 $i
#         $DIR/cmake-build-debug/pairhmm $p 40 40 $i
#
#         $DIR/cmake-build-debug/pairhmm $p 8 16 $i
#
#         $DIR/cmake-build-debug/pairhmm $p 8 24 $i
#         $DIR/cmake-build-debug/pairhmm $p 16 24 $i
#
#         $DIR/cmake-build-debug/pairhmm $p 8 32 $i
#         $DIR/cmake-build-debug/pairhmm $p 16 32 $i
#         $DIR/cmake-build-debug/pairhmm $p 24 32 $i
#
#         $DIR/cmake-build-debug/pairhmm $p 8 40 $i
#         $DIR/cmake-build-debug/pairhmm $p 16 40 $i
#         $DIR/cmake-build-debug/pairhmm $p 24 40 $i
#         $DIR/cmake-build-debug/pairhmm $p 32 40 $i
#
#         $DIR/cmake-build-debug/pairhmm $p 8 48 $i
#         $DIR/cmake-build-debug/pairhmm $p 16 48 $i
#         $DIR/cmake-build-debug/pairhmm $p 24 48 $i
#         $DIR/cmake-build-debug/pairhmm $p 32 48 $i
#         $DIR/cmake-build-debug/pairhmm $p 40 48 $i
#
#         $DIR/cmake-build-debug/pairhmm $p 8 56 $i
#         $DIR/cmake-build-debug/pairhmm $p 16 56 $i
#         $DIR/cmake-build-debug/pairhmm $p 24 56 $i
#         $DIR/cmake-build-debug/pairhmm $p 32 56 $i
#         $DIR/cmake-build-debug/pairhmm $p 40 56 $i
#     done
# done

for p in 32
do
    for i in $(seq 1 100)
    do
        $DIR/cmake-build-debug/pairhmm $p 20 20 $i
    done
done
