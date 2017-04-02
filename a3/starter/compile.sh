#!/bin/sh
g++ -O4 -g -fopenmp -std=c++11 svdDynamic.c RayTracer.c utils.c -lm -o RayTracer
