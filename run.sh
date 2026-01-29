#!/bin/bash

mkdir -p build/Release
cmake -DCMAKE_POLICY_VERSION_MINIMUM=3.5 -DRTAUDIO_API_JACK=OFF -DRTMIDI_API_JACK=OFF -DCMAKE_BUILD_TYPE=Release -Wno-deprecated -DBUILD_EXAMPLES=0 -B build/Release -S .

(
  # utilizing cmake's parallel build options
  # for cmake >= 3.12: -j <number of processor cores + 1>
  # for older cmake: -- -j5
  cmake --build build/Release --config Release -j 9
)

result=$?
if [ ${result} == 0 ]; then
  cd bin
  ./RAGECAGE
fi
