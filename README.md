# Neonatal EEG Analysis

This repository contains Matlab code for the analysis of neonatal EEG. 

## Algorithms

Burst Detection - Reference [1]

EEG Maturational Age Estimator (version 1) - Reference [2]

EEG Maturational Age Estimator (version 2) - Reference [3]

Seizure Detection (CNN) - Reference [4]

Preterm Maturational Features - from the literature

Stand-alone Feature Set for Preterm EEG Analysis

##


## EEG file format

EDF format (European Data Format)

see https://www.edfplus.info/specs/edf.html

## Planned Posts

| HIE grading algorithms

## Prerequisites

Matlab 2017a, 2018b, 2020b

## Example 

An example for the burst detection algorithm is demo_burst_detection.m

An example for the EMA estimator is estimate_ema.m

An example for seizure detection with CNNs is demo.m

An example for implementing the preterm feature set is demo_preterm.m

## Built With

Matlab 2017a, 2018b, 2020b

## Authors

Nathan Stevenson

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## References

[1] K Palmu, N Stevenson, S Wikström, L Hellström-Westas, S Vanhatalo, JM Palva. Optimization of an NLEO-based algorithm for automated detection of spontaneous activity transients in early preterm EEG. Physiol Meas, 2010; 31: N85-N93 

[2] O’Toole JM, Boylan GB, Vanhatalo S, Stevenson NJ. Estimating functional brain maturity in very and extremely preterm neonates using automated analysis of the electroencephalogram. Clinical Neurophysiology. 2016; 127(8): 2910-8.

[3] Stevenson NJ, Oberdorfer L, Koolen N, O’Toole JM, Werther T, Klebermass-Schrehof K, Vanhatalo S. Functional maturation in preterm infants measured by serial recording of cortical activity. Scientific Reports. 2017; 7(1): 12969.

[4] Stevenson NJ, Tapani K, Vanhatalo S. Hybrid neonatal EEG seizure detection algorithms achieve the benchmark of visual interpretation of the human expert. In Proc. IEEE EMBC, Berlin, Germany 2019.

## Contact

Nathan Stevenson

QIMR Berghofer, Previously at University of Helsinki

email: nathan.stevenson@qimrberghofer.edu.au

