# Calibrating The Meuse and Rhine Wflow Hourly Models

This renewed effort in calibrating the two largest operational hydrologic models for the Netherlands is intended to address:

1. Primarily timing lead and lag times in event onset and arrival on a per catchment basis.
2. Improving high and low flow simulation accuracy. 
3. Providing uncertainty estimates that stem from calibration. 
4. Reproducibility. 

A workflow is built for the Meuse, which is then deployed for the Rhine. 

This workflow is borrowing and expanding on the work done for:
1. Interreg (calibration, Laurene Bouaziz and Anais Cousanon)
2. Puget (Joost Buitink and Brendan Dalmijn)

# TODO: names and characters, licensing, projects any additional names

The innovations and contributions of this new effort in model calibration are:
1. Event timing metrics as calibration objectives (mean timing offset and std of lead/lag)
2. Cascading calibration (dependency tree solver)
3. 