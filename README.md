# Cloud-Algorithms
PEFT (OCT Rank) Algorithm 
Efficient application scheduling algorithms are important for obtaining high performance in heterogeneous computing
systems. In this paper, we present a novel list-based scheduling algorithm called Predict Earliest Finish Time (PEFT) for
heterogeneous computing systems. The algorithm has the same time complexity as the state-of-the-art algorithm for the same
purpose, that is, O(n^2logn):pÞ for v tasks and p processors, but offers significant makespan improvements by introducing a look-ahead
feature without increasing the time complexity associated with computation of an optimistic cost table (OCT). The calculated value is
an optimistic cost because processor availability is not considered in the computation. Our algorithm is only based on an OCT that is
used to rank tasks and for processor selection. The analysis and experiments based on randomly generated graphs with various
characteristics and graphs of real-world applications show that the PEFT algorithm outperforms the state-of-the-art list-based
algorithms for heterogeneous systems in terms of schedule length ratio, efficiency, and frequency of best results.
