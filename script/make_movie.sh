#!/bin/sh

IN=../build/png
OUT=../build/png

ffmpeg -framerate 5 -i "${IN}/m%03d.png" -c:v libx264 -r 24 -pix_fmt yuv420p "${OUT}/movie.mp4"
