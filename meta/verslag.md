 - Uitvoermethode werd aangepast
 - Gebruik van Rayon voor parallelisatie
 - Opkuisen van ongebruikte format parameter 
 - Gebruik van library whiteread & modulair maken van het inlezen van de training files


Benchmarks
==========

 10 runs whole genome (example/NC_000913.fna) met train/complete
    32.809s = 3.28s per run (geoptimaliseerd)
    39.947s = 3.99s per run (FragGeneScan1.31)

 10 runs non whole genome (example/NC_000913-454.fna) met train/454_10
    28.886s = 2.88s per run (geoptimaliseerd)
    1m43.831s = 10.3831s per run (FragGeneScan1.31)

UMGAP Resultaten
================
Gemiddeld over 10 runs.

    FGSRS:
    Sensitiviteit: 9291
    Precisie: 4596.5

    FGSPP:
    Sensitiviteit: 9260.6
    Precisie: 4460.8

    FGS1.31:
    Sensitiviteit: 9284
    Precisie: 4605.3

FragGeneScanRs scoort dus ongeveer even goed gemiddeld genomen over 10 runs als de originele FragGeneScan1.31 en ongeveer 160 punten beter dan FragGeneScan++.