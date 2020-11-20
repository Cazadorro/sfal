fn partial_partition_array(partition: &mut Vec<u64>, vec_partition: &mut Vec<Vec<u64>>, max_i: u64, k: u64) {
    // for i in 1..max_i{
    //     partition[i] = 0;
    // }
    // parition[0] = k-max_i;
    // vec_partition.push(partition);
    for i in 2..max_i {
        let index = (i - 1) as usize;
        while partition[0] >= i {
            partition[0] -= i;
            let prev = partition[0];
            partition[index] += 1;
            vec_partition.push(partition.to_vec());
            partial_partition_array(partition, vec_partition, i, partition[0]);
            partition[0] = prev;
        }
        partition[index] = 0;
        partition[0] = k;
    }
    //resetting
    for i in 0..max_i as usize {
        partition[i] = 0;
    }
}

pub fn integer_partition_array(k: u64) -> Vec<Vec<u64>> {
    let mut partition = vec![0_u64; k as usize];
    let mut vec_partition = Vec::<Vec::<u64>>::new();
    let last_index = (k - 1) as usize;
    partition[last_index] = 1;
    vec_partition.push(partition.to_vec());
    partition[last_index] = 0;
    partition[0] = k;
    vec_partition.push(partition.to_vec());
    partial_partition_array(&mut partition, &mut vec_partition, k, k);
    vec_partition
}

pub fn integer_partition(n: u64) -> Vec<Vec<u64>> {
    let mut a = vec![0_u64; n as usize];
    let mut k = 1_u64;
    let mut y = n - 1;
    let mut partition_vec = Vec::<Vec::<u64>>::new();
    while k != 0 {
        let mut x = a[(k - 1) as usize] + 1;
        k -= 1;
        while 2 * x <= y {
            a[k as usize] = x;
            y -= x;
            k += 1;
        }
        let mut l = k + 1;
        while x <= y {
            a[k as usize] = x;
            a[l as usize] = y;
            partition_vec.push((a[0..(k + 2) as usize]).to_vec());
            x += 1;
            y -= 1;
        }
        a[k as usize] = x + y;
        y = x + y - 1;
        partition_vec.push((a[0..(k + 1) as usize]).to_vec());
    }
    partition_vec
}

fn calc_p0(partitions: &Vec<u64>) -> u64 {
    let k = partitions.len() as u64;
    return k - partitions.iter().sum::<u64>();
}

pub fn factorial(value: u64) -> u64{
    let mut fact = 1_u64;
    for i in 2..=value {
        fact *= i;
    }
    fact
}


pub fn calc_r_k(phi_x0: f64, phi_dx_x0: &Vec<f64>, u: u64, v: u64, k: u64) -> f64 {
    if k == 0{
        1.0
    }else {
        let partitions = integer_partition_array(k);
        let mut sum = 0.0;
        for partition in (&partitions).iter() {
            let p0 = calc_p0(partition);
            let sign = if (k - p0) & 1 == 0 {
                1.0
            } else {
                -1.0
            };
            let mut outer_prod = 1.0;
            for r in 0..=(k - p0 - 1) {
                outer_prod *= (u as f64 / v as f64) + r as f64;
            }
            let mut inner_prod = 1.0;
            for s in 1..=k {
                let s_idx = (s - 1) as usize;
                let ps = partition[s_idx];
                let psfact = factorial(partition[s_idx]) as f64;
                inner_prod /= psfact;
                let sfact = factorial(s) as f64;
                inner_prod *= (phi_dx_x0[s_idx] / (sfact * phi_x0)).powi(ps as i32);
            }
            sum += sign * outer_prod * inner_prod;
        }
        sum
    }
}

pub fn calc_burmann_coefficent(phi_x0: f64, phi_dx_x0: &Vec<f64>, n: u64, v: u64) -> f64 {
    let mut sum = 0.0;
    for r in 0..=(n - 1) {
        let r_idx = r as usize;
        let r_fact = factorial(r) as f64;
        sum += (phi_dx_x0[r_idx] * calc_r_k(phi_x0, phi_dx_x0, n, v, n-r-1))/(r_fact * n as f64);
    }
    sum
}


