import sys
sys.path.insert(0, 'sixcircle_1p85')

missing = []
try:
    from sixcircle_rqd import *
    print("SUCCESS - all dependencies installed!")
except ModuleNotFoundError as e:
    module = str(e).split("'")[1]
    print(f"Missing: {module}")
    print(f"\nInstall with: pip install {module}")
