#!/bin/bash

kubectl cp watch-evaluator-lw-pod:/evaluator-storage/ evaluator-storage

du -hsc evaluator-storage

tail -n 1 evaluator-storage/working-directory/SimulationLayer/*/*/openmm_statistics.csv


