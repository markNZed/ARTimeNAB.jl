# ARTimeNAB

This is an open source (AGPL3) and simplified version of ARTime, an anomaly detection algorithm.
It supports the Numenta anomaly benchmark ([NAB](https://github.com/numenta/NAB)) ARTime detector.

A fork of [NAB](https://github.com/markNZed/NAB/tree/ARTime) includes the ARTime detector, please start there to see & reproduce the ARTime results. The NAB Python package in that fork uses [PythonCall](https://github.com/cjdoris/PythonCall.jl) to install ARTime from this repository.

ARTime was developed by Mark Hampton.

# Acknowledgements

Stephen Grossberg and Gail Carpenter developed adaptive resonance theory (ART). Grossberg's 2021 book *Conscious Mind, Resonant Brain: How Each Brain Makes a Mind* was the major inspiration for ARTime.

Numenta provided NAB to inspire innovation in anomaly detection. It was very valuable in testing ARTime. The paper introducing NAB is from Ahmad, S., Lavin, A., Purdy, S., & Agha, Z. (2017). Unsupervised real-time anomaly detection for streaming data. Neurocomputing, Available online 2 June 2017, ISSN 0925-2312, https://doi.org/10.1016/j.neucom.2017.04.070

The excellent Julia package [AdaptiveResonance.jl](https://github.com/AP6YC/AdaptiveResonance.jl) was extremely useful in getting ARTime off the ground. Modifications in the DVFA implementation of AdaptiveResonance.jl led to a compact version of AdaptiveResonance being included in ARTime. 

[@isentropic](https://github.com/isentropic) was a great help in introducing me to Julia and improving the quality of the code.

# Where to from here

Unfortunately deep learning catastrophically forgot why it was not a panacea in the 1980s. A talented team of computational neuroscientists could push ART much further in machine learning... For an introduction to Grossbergian Neuroscience look no further than Yohan John's [Neurologos](https://www.youtube.com/watch?v=HrOxj-3hBiw&list=PLTEtXsHFKZTsxmKyVn69ZmLBxghBH1tNR) channel on YouTube.
