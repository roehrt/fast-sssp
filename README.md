# Fast Negative-Weight SSSP

This repository provides a standalone C++ implementation of the fast single-source shortest paths (SSSP) algorithm by Bringmann, Cassis, and Fischer [^2]. The implementation is based on the foundational algorithm developed by Bernstein, Nanongkai, and Wulff-Nilsen [^1], with additional insights drawn from Li and Mowry [^3].

The runtime complexity described in [^2] is $\mathcal{O}((m + n \log \log n) \log^2 n \log(nW))$, where $-W$ is the smallest edge weight. This complexity is achieved using Thorup's integer priority queue [^4], but implementing such a priority queue would be a substantial project on its own.

Our implementation offers two alternatives for priority queues:
- **Thin Heap:** If `-DTHIN_HEAP` is defined during compilation, the implementation uses [GNU's policy-based thin heap](https://gcc.gnu.org/onlinedocs/libstdc++/ext/pb_ds/pq_performance_tests.html#thin_heap_note). This results in a runtime complexity of $\mathcal{O}((m + n \log n) \log^2 n \log(nW))$.
- **Pairing Heap (Default):** Without `-DTHIN_HEAP`, the implementation defaults to using a pairing heap, which achieves a runtime complexity of $\mathcal{O}(m \log^3 n \log(nW))$. Despite having a worse theoretical complexity, the pairing heap is faster in practice.


_This is primarily a fun project, and it's important to note that a non-naive implementation of classic algorithms like Bellman-Ford will generally be faster in practice. Nonetheless, unlike many other recent theoretical breakthrough algorithms for classic problems, this one seems to be actually somewhat usable in practice._

## Building

To build the project, you will need `g++`.

### Using CMake

Follow these steps to build using CMake:
```bash
cmake -S . -B build && cmake --build build
```
The compiled binary will be available at `build/sssp`.

### Using `g++`

Alternatively, you can build the project directly with `g++`:
```bash
g++ -O3 -march=native -std=c++20 main.cpp -o sssp
```

## References

[^1]: Bernstein, Aaron, Danupon Nanongkai, and Christian Wulff-Nilsen. 2022.  
“Negative-Weight Single-Source Shortest Paths in Near-Linear Time.” In  
*63rd IEEE Annual Symposium on Foundations of Computer Science, FOCS 2022, Denver, CO, USA, October 31 - November 3, 2022*, 600–611. IEEE.  
<https://doi.org/10.1109/FOCS54457.2022.00063>.

[^2]: Bringmann, Karl, Alejandro Cassis, and Nick Fischer. 2023.  
“Negative-Weight Single-Source Shortest Paths in Near-Linear Time: Now Faster!” In  
*64th IEEE Annual Symposium on Foundations of Computer Science, FOCS 2023, Santa Cruz, CA, USA, November 6-9, 2023*, 515–538. IEEE.  
<https://doi.org/10.1109/FOCS57990.2023.00038>.

[^3]: Li, Jason, and Connor Mowry. 2024.  
“A Bottom-up Algorithm for Negative-Weight SSSP with Integrated Negative Cycle Finding.”  
*CoRR* abs/2411.19449.  
<https://doi.org/10.48550/ARXIV.2411.19449>.

[^4]: Thorup, Mikkel. 2003. “Integer Priority Queues with Decrease Key in
Constant Time and the Single Source Shortest Paths Problem.” In
*Proceedings of the 35th Annual ACM Symposium on Theory of Computing,
June 9-11, 2003, San Diego, CA, USA*, edited by Lawrence L. Larmore and
Michel X. Goemans, 149–58. ACM. <https://doi.org/10.1145/780542.780566>.