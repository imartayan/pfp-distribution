# PFP distribution

Visualize the distribution of PFP phrase sizes

```sh
# plot PFP distribution on <input.fa> for w=10 and p=32
cargo r -r -- -w 10 -p 32 -i <input.fa> | python plot-pfp.py

# skipping the input will use a random sequence instead
cargo r -r -- -w 10 -p 32 | python plot-pfp.py

# plot multiple values of w at once
cargo r -r -- -w 10 21 31 -p 32 -i <input.fa> | python plot-pfp.py

# use an alternative version of NtHash (rotate by 7)
cargo r -r -F alt -- -w 10 -p 32 -i <input.fa> | python plot-pfp.py

# use the multiplicative version of NtHash (MulHash)
cargo r -r -F mul -- -w 10 -p 32 -i <input.fa> | python plot-pfp.py
```
