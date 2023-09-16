# LSTest

# Compile and Run

```
g++ rain2r.cpp -o rain -O3
./rain

g++ brute_force.cpp -o bf -O3
./bf
```

# Complexity

We conducted experiments on the 128, 192 and 256-bit versions.  Both of our approach and brute force follow the Guess-and-Determine framework, which allowed us to estimate the complexity by multiplying the duration of each trial by the total number of trials. 

Experimental results revealed that our method requires 3.1ms, 13.3ms, and 26.4ms per trial, with the number of trials being 2^120, 2^176, and 2^240, respectively. For the brute force method, it requires 0.35ms, 0.9ms, and 1.8ms per trial, with the number of trials being 2^128, 2^192, and 2^256 respectively.

Thus the recalculated actual complexities are 2^123.1, 2^179.9, and 2^243.8. Meanwhile, our theoretical estimates of complexity were 2^120.3, 2^180.4, and 2^243.1.  The experimental outcomes align well with the theoretical predictions (with the 192-bit version even doing better).
