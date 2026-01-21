use crate::ObservatoryConfigParams;

pub fn print_observatory_config<R>(params: &ObservatoryConfigParams<R>)
where
    R: std::fmt::Display + std::fmt::LowerExp,
{
    println!("┌──────────────────────────────────────────────────────────────────────────┐");
    println!("│  OBSERVATORY CONFIGURATION                                               │");
    println!("├──────────────────────────────────────────────────────────────────────────┤");
    println!("│  Temporal Coordinates:                                                   │");
    println!("│    T1: E12 (circular orbit) - Reference baseline                         │");
    println!("│    T2: E18 (eccentric orbit) - Cross-correlation baseline                │");
    println!("│    T3: E14 at t₀ (eccentric) - Dynamic clock                             │");
    println!("│    T4: E14 at t₀+1h - Temporal offset 1                                  │");
    println!("│    T5: E14 at t₀+2h - Temporal offset 2                                  │");
    println!("├──────────────────────────────────────────────────────────────────────────┤");
    println!("│  Physical Constants:                                                     │");
    println!(
        "│    G (NIST):     {:<12.5e} m³/(kg·s²)                                 │",
        params.g_ref
    );
    println!(
        "│    c:            {:<12.0} m/s                                        │",
        params.c
    );
    println!(
        "│    GM_Earth:     {:<12.9e} m³/s²                                    │",
        params.gm_earth
    );
    println!(
        "│    Time delay:   {:<4} seconds ({:<1} hours)                                  │",
        params.time_delay_s,
        params.time_delay_s / 3600
    );
    println!("└──────────────────────────────────────────────────────────────────────────┘");
}
