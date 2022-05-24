#!/bin/bash

# 1. Quickstart ---------------------------------------------------------------
echo "* [1/9] Quickstart"
python Demo_Quickstart.py

# 2,3. Landscape Creation and Update ------------------------------------------
echo "* [2,3/9] Landscape Creation and Update"
python Demo_XY.py

# 4. Sites and Trap Types -----------------------------------------------------
echo "* [4/9] Landscape Creation"
python Demo_Types.py

# 5. GA Optimization ----------------------------------------------------------
echo "* [5/9] GA Optimization"
python Demo_GA.py
python Demo_GA-Simple.py

# 6. GA with Immovable Traps --------------------------------------------------
echo "* [6/9] GA with Immovable Traps"
python Demo_GACustom.py
python Demo_GACustom-Simple.py

# 7. GA Custom with Multi-Point Type ------------------------------------------
echo "* [7/9] GA Custom with Multi-Point Type"
python Demo_pointTypes.py
python Demo_pointTypes-Simple.py

# 8. GA with Sex Kernels ------------------------------------------------------
echo "* [8/9] GA with Sex-Specific Traps"
python Demo_GASex.py
python Demo_GASex-Simple.py

# 9. GA with Immovable Traps---------------------------------------------------
echo "* [9/9] Stage-Specific Traps"
python Demo_TrapsComplex.py