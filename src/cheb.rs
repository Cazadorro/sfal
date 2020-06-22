use std::f64;

///https://www.embeddedrelated.com/showarticle/152.php
/// x is your normal x, u is the normalized coordinate
fn u_from_x(x: f64, a: f64, b: f64) -> f64 {
    ((2.0 * x) - a - b) / (b - a)
}

fn x_from_u(u: f64, a: f64, b: f64) -> f64 {
    ((b - a) / 2.0) * u + ((a + b) / 2.0)
}

fn chebyshev_node(i: u64, n: u64) -> f64 {
    return (((2.0 * i as f64 - 1.0) * std::f64::consts::PI) / (2.0 * n as f64)).cos();
}

fn chebyshev_nodes(node_count:u64) -> Vec<f64>{
    (1..=node_count).map(|i| chebyshev_node(i, node_count)).collect::<Vec<f64>>()
}

fn craft_approximator(x_range_min: f64, x_range_max: f64, c_values :Vec<f64>) -> impl Fn(f64) -> f64{
    let approximator = move |x: f64| -> f64{
        // let c_values = c_values;
        // let x_range_min= x_range_min;
        // let x_range_max= x_range_max;
        let u = u_from_x(x, x_range_min, x_range_max);
        let mut t_prev = 1.0;
        let mut t = u;
        let mut y = c_values[0];
        for i in 1..c_values.len() {
            y += t * c_values[i];
            let t_next = 2.0 * u * t - t_prev;
            t_prev = t;
            t = t_next;
        }
        y
    };
    approximator
}

pub fn gen_chebyshev_approximator(x_range_min: f64, x_range_max: f64, node_count:u64, fx: fn(f64) -> f64) -> impl Fn(f64) -> f64{
    let nodes = chebyshev_nodes(node_count);
    let x_values = nodes.iter().map(|u| x_from_u(*u, x_range_min, x_range_max)).collect::<Vec<f64>>();
    let y_values = x_values.iter().map(|x| (fx)(*x)).collect::<Vec<f64>>();
    println!("{:#?}",nodes);
    println!("{:#?}",x_values);
    println!("{:#?}",y_values);
    let mut c_values = vec![0.0; node_count as usize];
    let mut t_values_0 = vec![1.0; node_count as usize];
    let mut t_values_1 = nodes.clone();
    if node_count > 0 {
        c_values[0] = y_values.iter().sum::<f64>() / node_count as f64;
        if node_count > 1 {
            c_values[1] = 2.0 * t_values_1.iter().zip(&y_values).map(|(t, y)| (*t) * (*y)).sum::<f64>() / node_count as f64;
            for i in 2..node_count as usize {
                //tn+1 = 2*u*tn - tn-1;
                let new_t_values = if i & 1 == 0 {
                    for j in 0..node_count as usize {
                        let t_value_prev = t_values_0[j];
                        let t_value = t_values_1[j];
                        t_values_0[j] = (2.0 * nodes[j] * t_value) - t_value_prev;
                    }
                    &t_values_0
                }else{
                    for j in 0..node_count as usize {
                        let t_value_prev = t_values_1[j];
                        let t_value = t_values_0[j];
                        t_values_1[j] = (2.0 * nodes[j] * t_value) - t_value_prev;
                    }
                    &t_values_1
                };
                c_values[i] = 2.0 * new_t_values.iter().zip(&y_values).map(|(t, y)| (*t) * (*y)).sum::<f64>() / node_count as f64;
            }
        }
    }
    craft_approximator(x_range_min, x_range_max, c_values)
}