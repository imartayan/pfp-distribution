use clap::Parser;
use core::array::from_fn;
use needletail::parse_fastx_file;
use rayon::prelude::*;
use rayon::{ThreadPoolBuilder, current_num_threads};
use serde::Serialize;
use std::time::Instant;

use simd_minimizers::collect::CollectAndDedup;
use simd_minimizers::packed_seq::{PackedSeqVec, SeqVec};
use simd_minimizers::private::S;
use simd_minimizers::seq_hash::KmerHasher;

#[cfg(feature = "mulhash")]
type SIMDHasher = simd_minimizers::seq_hash::MulHasher<false>;
#[cfg(all(feature = "nthash", not(feature = "mulhash")))]
type SIMDHasher = simd_minimizers::seq_hash::NtHasher<false, 1>;
#[cfg(not(any(feature = "mulhash", feature = "nthash")))]
type SIMDHasher = simd_minimizers::seq_hash::NtHasher<false>;

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Input file (FASTA, possibly compressed) [default: random seq]
    #[arg(short, long)]
    input: Option<String>,
    /// Output file for histogram (JSON) [default: stdout]
    #[arg(short, long)]
    output: Option<String>,
    /// Window size
    #[arg(short, num_args = 1..)]
    w: Vec<usize>,
    /// Sampling factor of the windows
    #[arg(short, num_args = 1..)]
    p: Vec<usize>,
    /// Hash seed [default: none]
    #[arg(short, long)]
    seed: Option<u32>,
    /// Number of threads [default: all]
    #[arg(short, long)]
    threads: Option<usize>,
}

fn main() {
    const DEFAULT_LEN: usize = 100_000_000;

    let args = Args::parse();
    let threads = if let Some(t) = args.threads {
        ThreadPoolBuilder::new()
            .num_threads(t)
            .build_global()
            .unwrap();
        t
    } else {
        current_num_threads()
    };

    let seqs = if let Some(path) = args.input {
        let mut seqs = vec![PackedSeqVec::default(); threads];
        let mut i = 0;
        let mut reader = parse_fastx_file(path).unwrap();
        while let Some(record) = reader.next() {
            seqs[i].push_ascii(&record.unwrap().seq());
            i += 1;
            if i == threads {
                i = 0;
            }
        }
        seqs
    } else {
        (0..threads)
            .map(|_| PackedSeqVec::random(DEFAULT_LEN))
            .collect()
    };

    let input_size = seqs.iter().map(|seq| seq.len()).sum::<usize>();
    let start = Instant::now();

    let mut hists = Vec::with_capacity(args.w.len() * args.p.len());
    args.p.iter().for_each(|&p| {
        args.w.iter().for_each(|&w| {
            hists.push(compute_hist(&seqs, w, p, args.seed));
        })
    });

    let time = start.elapsed().as_secs_f64();
    eprintln!(
        "computed PFP distribution in {time:.2}s - {:.2} GB/s",
        (input_size * args.w.len() * args.p.len()) as f64 / 1e9 / time
    );

    println!("{}", serde_json::to_string(&hists).unwrap());
}

#[derive(Serialize)]
struct Hist {
    w: usize,
    p: usize,
    total: usize,
    hist: Vec<u32>,
}

fn compute_hist(seqs: &[PackedSeqVec], w: usize, p: usize, seed: Option<u32>) -> Hist {
    let mut dists: Vec<Vec<u32>> = Vec::with_capacity(seqs.len());
    let mut hist = vec![0; p * 3];

    seqs.par_iter()
        .map(|seq| {
            let hasher = match seed {
                Some(seed) => SIMDHasher::new_with_seed(w, seed),
                None => SIMDHasher::new(w),
            };
            let hash_iter = hasher.hash_kmers_simd(seq.as_slice(), 1);
            let threshold = S::new([u32::MAX / p as u32 + 1; 8]);
            let offset = hash_iter.it.len();
            let offsets = from_fn(|i| (i * offset) as u32);
            let mut positions = Vec::with_capacity(offset * 9 / p);

            let mut pos = S::new(offsets);
            let mut selected_pos = pos;
            let selected_pos_iter = hash_iter.map(|hash| {
                let is_selected = hash.cmp_lt(threshold);
                selected_pos = is_selected.blend(pos, selected_pos);
                pos += S::ONE;
                selected_pos
            });
            selected_pos_iter.collect_and_dedup_into::<false>(&mut positions);
            for i in 0..(positions.len() - 1) {
                positions[i] = positions[i + 1] - positions[i];
            }
            positions.pop();
            positions
        })
        .collect_into_vec(&mut dists);

    let total = dists.iter().map(|v| v.len()).sum::<usize>();
    dists.iter().for_each(|dist| {
        dist.iter().map(|&d| d as usize).for_each(|d| {
            if d < hist.len() {
                hist[d] += 1;
            }
        })
    });
    Hist { w, p, total, hist }
}
