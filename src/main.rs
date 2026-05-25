use anyhow::Result;
use itertools::Itertools;

use crate::space::{ProbabilityMap, ThresholdedSpace};

mod grid;
mod space;

fn main() -> Result<()> {
    let mut pmap = ProbabilityMap::empty(5, 10);

    // Obstacle
    pmap.update((1, 4), 1.0)?;
    pmap.update((1, 3), 1.0)?;

    // Obstacle
    pmap.update((3, 7), 1.0)?;
    pmap.update((3, 6), 1.0)?;

    // Obstacle
    pmap.update((3, 1), 1.0)?;

    let tspace = ThresholdedSpace::from_probability_map(&pmap, 0.5);
    println!("Thresholded Space:");
    println!("{}", tspace);

    for sig in (0..3).map(|_| 0..2).multi_cartesian_product() {
        println!(
            "Path with sig [{}]:",
            sig.iter().map(|c| c.to_string()).join(",")
        );
        if let Some(path) = tspace.a_star((4, 0), (0, 9), sig) {
            tspace.print_with_path(&path);
        } else {
            println!("no path found");
        }
    }

    Ok(())
}
