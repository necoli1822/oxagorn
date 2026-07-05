// Quick cross-check: print detected tRNA/tmRNA genes via the safe API.
use oxagorn::api::{detect, SearchOptions};

fn main() {
    let path = std::env::args().nth(1).expect("usage: detect <genome.fna>");
    let opts = SearchOptions::default(); // tRNA+tmRNA, both strands, circular, standard code
    let genes = detect(&path, &opts, 1).expect("detect");
    let mut trna = 0;
    let mut tmrna = 0;
    for g in &genes {
        let kind = if g.genetype == 1 {
            tmrna += 1;
            "tmRNA"
        } else {
            trna += 1;
            "tRNA"
        };
        println!(
            "{}\t{}\t{}..{}\t{}\tE={:.1}\t{}",
            g.seqname,
            kind,
            g.start,
            g.stop,
            if g.comp == 1 { "-" } else { "+" },
            g.energy,
            g.anticodon
        );
    }
    eprintln!("API: {} genes ({} tRNA, {} tmRNA)", genes.len(), trna, tmrna);
}
