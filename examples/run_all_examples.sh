#!/bin/bash
# Run all example scripts

echo "Running Example 1: Single Q-point analysis"
python example_01_single_q.py

echo ""
echo "Running Example 2: Multiple Q-points"
python example_02_multiple_q.py

echo ""
echo "Running Example 3: Low-level API"
python example_03_lowlevel_api.py

echo ""
echo "All examples completed!"
