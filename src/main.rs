use anyhow::Result;
use itertools::Itertools;

use crate::{
    persistence::Solver,
    space::{ProbabilityMap, ThresholdedSpace},
};

mod grid;
mod persistence;
mod space;

fn main() -> Result<()> {
    let mut pmap = ProbabilityMap::empty(5, 5);
    pmap.update((1, 1), 0.3)?;
    pmap.update((1, 2), 0.5)?;
    pmap.update((3, 3), 0.9)?;

    let solver = Solver::builder()
        .pmap(pmap)
        .delta(0.2)
        .start((0, 0))
        .goal((4, 4))
        .build();
    let bc = solver.barcode();
    println!("{}", bc);

    // TODO
    // - Better visualization for the barcode
    // - Test the barcode generation
    // - Test hausdorff distance

    Ok(())
}
