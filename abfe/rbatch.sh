#!/bin/sh
until sbatch $1; do echo "Submission failed - retrying..." ; done
