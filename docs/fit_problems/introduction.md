# Introduction

This section serves as a tutorial on how to treat problems that frequently appear when running profile likelihood (PL) fits with a focus on TRExFitter.
This tutorial does not intent to focus on the basics of the PL fit nor on the basic configuration of TRExFitter.
It is assumed that the people doing this tutorial are familiar with the TRExFitter framework and they know how to create histograms from ntuples/read histograms/run the fits, etc.
Part of this tutorial is based on the twiki from exotics group, available [here](https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/Profilelikelihood).

Useful links:

- If you want to read about PL, see [slides by Nicolas Berger](https://indico.cern.ch/event/664144/contributions/2728127/attachments/1558588/2452300/HBSMWS-20171115-expanded.pdf)
- Basic documentation about RooStats (the framework that is doing the heavy lifting in terms of technicalities) is available on cds [here](https://cds.cern.ch/record/1456844)
- Minuit documentation (minimiser) [pdf](https://root.cern.ch/download/minuit.pdf)
- TRExFitter documentation is available on this [twiki](https://twiki.cern.ch/twiki/bin/view/AtlasProtected/TtHFitter)
    - For the list of available commands and general documentation of TRExFitter, see the [settings](../settings.md)
    - Some basic tutorials for TRExFitter are available [15 min TRExFitter introduction](https://indico.cern.ch/event/684782/contributions/2957566/attachments/1634431/2607059/20180417_TRExFitter_overview.pdf) from April 17, 2018 and [90 min TRExFitter tutorial](https://indico.cern.ch/event/700646/contributions/2936548/attachments/1627872/2595678/TRExFitterTutorial.pdf) from April 5, 2018
- Effects of smoothing on constraints: [presentation](https://indico.cern.ch/event/761804/contributions/3160985/attachments/1733339/2802398/Defranchis_template_constraints.pdf)
- Smoothing and pruning [presentation](https://indico.cern.ch/event/691683/contributions/2873279/attachments/1593521/2522846/PruningSmoothing.pdf)
- SUSY group tips and tricks [twiki](https://twiki.cern.ch/twiki/bin/view/AtlasProtected/SUSYMultiBinRecommendations)

- Contact:
    - Mailing list: [atlas-phys-stat-tthfitter@cern.ch](https://e-groups.cern.ch/e-groups/EgroupsSubscription.do?egroupName=atlas-phys-stat-tthfitter)
    - JIRA: [TRExFitter JIRA](https://its.cern.ch/jira/browse/TTHFITTER/)
