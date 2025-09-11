# PFP distribution

Visualize the distribution of PFP phrase sizes

```sh
# plot PFP distribution on <input.fa> for w=10 and p=100
cargo r -r -- -w 10 -p 100 -i <input.fa> | python plot-pfp.py

# skipping the input will use a random sequence instead
cargo r -r -- -w 10 -p 100 | python plot-pfp.py

# plot multiple values of w at once
cargo r -r -- -w 10 21 31 -p 100 -i <input.fa> | python plot-pfp.py

# use the original version of NtHash (biased)
cargo r -r -F nthash -- -w 10 -p 100 -i <input.fa> | python plot-pfp.py

# use the multiplicative version of NtHash (MulHash)
cargo r -r -F mulhash -- -w 10 -p 100 -i <input.fa> | python plot-pfp.py
```

<img width="640" height="480" alt="plot_human" src="https://github.com/user-attachments/assets/8065067e-e265-4f0e-82d4-076a52162f90" />
