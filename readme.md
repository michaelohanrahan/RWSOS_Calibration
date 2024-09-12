# Wflow Full Automatic Calibration Workflow

## Context

In comparison to lumped models, Wflow has proven to be difficult to calibrate. Two reasons are evident for this; firstly, the model is fully distributed and secondly, lumped or semi-distributed models can be run much faster unlocking the monte-carlo, 'brute force' -type calibration approaches for parameter estimation. Wflow does have the benefit, however of deriving parameters from oben data (e.g. earth observations) that allow for surprising performance unadjusted. We calibrate the models with this more physically-meaningful basis in mind and assume that we can allow for slightly runs to get to the same or better performance, bringing a workflow that resembles an adapted brute force approach within the reach of realistic computation times on a high performance computer(HPC). 

The distribution is still a problem and historically this has been dealt with in a multitude of ways. Whole basin, or per subcatchment adjustments in a manual or semi-automated way has typically been the approach, adjusting the most sensitive to the least sensitive parameters first and assessing individually. This approach carries a high risk of equifinality, a risk we aim to reduce. We accomplish here a calibration that solves the graph of dependent subcatchments, brute forcing runs within a parameter space per grouped catchment level to output a single best or ensemble of Wflow models assessed on a per subcatchment basis. 

### Aim:

We aim to produce a workflow that can be utilised in any project where the user with a base wflow model and complementary observations can calibrate to a varying intensity over their chosen parameter space.

### Case: Improving The Meuse and Rhine Wflow Hourly Model Performance

This renewed effort in calibrating the two largest operational hydrologic models for the Netherlands is intended to address:

    1. Primarily timing lead and lag times in event onset and arrival on a per catchment basis.
    2. Improving high and low flow simulation accuracy. 
    3. Providing uncertainty estimates that stem from calibration. 
    4. Reproducibility. 

A workflow is built and tested with the Meuse, which is then deployed for the Rhine and is essentially transferrable to any other wflow model.

This workflow is borrowing from and expanding on the work done for:

    1. GRADE & Rijkswaterstaat Operational (Joost & Laurene)
    2. EU Interreg (calibration, Laurene Bouaziz and Anais Cousanon)
    3. USGS Puget Sound (Joost Buitink and Brendan Dalmijn)

The innovations and contributions of this new effort in model calibration are:
    1. Event timing metrics as calibration objectives (mean timing offset and std of lead/lag)
    2. Solving the dependency (levels) of each catchment in a subcatch file
    3. Establishing uncertainty bounds by producing a variable amount of output models. 
   
The folder as it looks now contains three key folders:

    1. meuse:
       - The basic calibration. Using a minimal calibration recipe, we tested this accomplishing 1350 runs per level, outputting the best params per subcatchment. 
    2. meuse_random:
       - 